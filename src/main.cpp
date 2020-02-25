#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <regex>
#include <stdexcept>
#include <cstdio>

// C fileno
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

// for timing
#include <chrono>

// local header
#include "utils.hpp"
#include "split_iterator.hpp"
#include "compress.hpp"
#include "sparse.hpp"
#include "string_t.h"

int usage() {
    std::cerr << "./main [compress|decompress|sparsify] <input_file> <output_file>" << std::endl;
    return 1;
}

class VcfCoordinateQuery {
public:
    VcfCoordinateQuery(
            const std::string& reference_name,
            uint64_t start_position,
            uint64_t end_position) {
        this->reference_name = reference_name;
        this->start_position = start_position;
        this->has_start_position = true;
        this->end_position = end_position;
        this->has_end_position = true;
    }

    VcfCoordinateQuery() {
        this->reference_name = "";
        this->start_position = -1;
        this->end_position = -1;
        this->has_start_position = false;
        this->has_end_position = false;
    }

    bool matches(const std::string& reference_name, uint64_t position) {
        if (this->reference_name.size() > 0 && this->reference_name != reference_name) {
            return false;
        }
        if (this->has_start_position && this->start_position > position) {
            return false;
        }
        if (this->has_end_position && this->end_position < position) {
            return false;
        }
        return true;
    }

    bool has_criteria() {
        return this->reference_name.size() > 0
            || this->start_position >= 0
            || this->end_position >= 0;
    }

    const std::string& get_reference_name() {
        return this->reference_name;
    }

    uint64_t get_start_position() {
        return this->start_position;
    }

    uint64_t get_end_position() {
        return this->end_position;
    }

private:
    std::string reference_name;
    uint64_t start_position;
    bool has_start_position;
    uint64_t end_position;
    bool has_end_position;
};

class VcfLineStateMachine {
public:
    enum State {
        UNINITIALIZED, META, HEADER, VARIANT
    };

    VcfLineStateMachine() {
        current_state = UNINITIALIZED;
    }

    void to_meta() {
        if (current_state == META) {
            return;
        }
        else if (current_state == HEADER || current_state == VARIANT) {
            throw std::runtime_error("Cannot move to line state META");
        }
        current_state = META;
    }

    void to_header() {
        if (current_state == HEADER) {
            return;
        }
        // can't go backwards, and can't repeat the header line
        else if (current_state == VARIANT || current_state == HEADER) {
            throw std::runtime_error("Cannot move to line state VARIANT");
        }
        current_state = HEADER;
    }

    void to_variant() {
        if (current_state == VARIANT) {
            return;
        }
        else if (current_state == UNINITIALIZED || current_state == META) {
            throw std::runtime_error("Cannot move to line state VARIANT");
        }
        current_state = VARIANT;
    }
private:
    State current_state;
};

byte_array byte_vector_to_bytearray(const std::vector<byte_t>& v) {
    byte_array ba;
    ba.bytes = new byte_t[v.size()];
    ba.len = v.size();
    for (size_t i = 0; i < v.size(); i++) {
        ba.bytes[i] = v[i];
    }
    return ba;
}


void query_sparse_file_fd(const std::string& input_filename, VcfCoordinateQuery query) {
    int input_fd = open(input_filename.c_str(), O_RDONLY);
    if (input_fd < 0) {
        perror("open");
        throw std::runtime_error("Failed to open file: " + input_filename);
    }

    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);

    #ifdef DEBUG
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    start = std::chrono::steady_clock::now();
    #endif
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);
    #ifdef DEBUG
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    debugf("decompress2_metadata_headers time: %lu\n", duration.count());
    #endif

    // Leave default sparse config
    SparsificationConfiguration sparse_config;

    long off = lseek(input_fd, 0, SEEK_CUR);
    if (off < 0) {
        throw std::runtime_error("ftell failed: " + std::to_string(off));
    }
    long data_start_offset = off + 8;
    uint64_t first_line_offset = 0;
    if (read(input_fd, &first_line_offset, sizeof(uint64_t) == 0)) {
        throw std::runtime_error("Failed to read first_line_offset value from file");
    }

    debugf("data_start_offset = %lu\n", data_start_offset);
    debugf("first_line_offset = %lu\n", first_line_offset);

    // Single variant lookup
    if (query.has_criteria() && (query.get_start_position() == query.get_end_position())) {
        debugf("Single variant lookup\n");
        size_t variant_offset = sparse_config.compute_sparse_offset(query.get_reference_name(), query.get_start_position());

        debugf("variant_offset = %lu, file_offset = %lu\n", variant_offset, data_start_offset + variant_offset);

        long initial_lookup_offset = lseek(input_fd, data_start_offset + variant_offset, SEEK_SET);

        union {
            struct {
                uint64_t distance_to_previous;
                uint64_t distance_to_next;
            } lengths;
            uint8_t bytes[16];
        } _length_headers;
        if (read(input_fd, &_length_headers.bytes, 16) == 0) {
            throw std::runtime_error("Reached end of file unexpectedly");
        }
        debugf("distance_to_previous = %lu, distance_to_next = %lu\n",
            _length_headers.lengths.distance_to_previous, _length_headers.lengths.distance_to_next);

        // IF prev offset value is zero and this is not the first line, must be an invalid location
        if (_length_headers.lengths.distance_to_previous == 0 &&
                initial_lookup_offset != (long)(first_line_offset + data_start_offset)) {
            // seek ahead to next viable line start
            long seek_distance = sparse_config.multiplication_factor * sparse_config.block_size;
            seek_distance -= 16; // Already read this many bytes
            lseek(input_fd, seek_distance, SEEK_CUR);
            debugf("Offset %ld was not a data line for single variant lookup, output no data\n",
                tellfd(input_fd) + 16 - seek_distance);

        } else {
            debugf("Found requested single variant line\n");
            // std::string linebuf;
            // linebuf.reserve(4 * 1024); // 4 KiB
            string_t linebuf;
            string_reserve(&linebuf, 4 * 1024);
            size_t linelength;
            decompress2_data_line_fd2(input_fd, schema, &linebuf, &linelength);
            for (size_t i = 0; i < linebuf.size; i++) {
                printf("%c", linebuf.buf[i]);
            }
        }
    }
    // Multi-variant lookup
    else if (query.has_criteria() && (query.get_start_position() != query.get_end_position())) {
        debugf("Multiple variant lookup\n");
        size_t start_variant_offset = sparse_config.compute_sparse_offset(query.get_reference_name(), query.get_start_position());

        debugf("start of range: variant_offset = %lu, file_offset = %lu\n", start_variant_offset, data_start_offset + start_variant_offset);

        // seek to start of range
        long initial_lookup_offset = lseek(input_fd, data_start_offset + start_variant_offset, SEEK_SET);

        long initial_seek_data = lseek(input_fd, initial_lookup_offset, SEEK_DATA);

        debugf("initial_lookup_offset = %ld, initial_seek_data = %ld\n",
            initial_lookup_offset, initial_seek_data);

        // Landed in a sparse hole
        if (initial_lookup_offset != initial_seek_data) {
            debugf("SEEK_DATA moved from initially requested offset\n");
            // Update to correspond to a viable start-of-line offset
            // Below is the viable offset modulo for offsets minus data_start_offset
            long viable_line_offset_modulo = sparse_config.multiplication_factor * sparse_config.block_size;
            if ((initial_seek_data - data_start_offset) % viable_line_offset_modulo != 0) {
                // go backwards to previous viable start offset
                //long previous_viable_line_distance = (initial_seek_data - data_start_offset) % viable_line_offset_modulo;
                long next_viable_line_distance = viable_line_offset_modulo - ((initial_seek_data - data_start_offset) % viable_line_offset_modulo);
                //previous_viable_line_distance += data_start_offset;

                // lseek(input_fd, -previous_viable_line_distance, SEEK_CUR);
                // debugf("Seeked backwards previous_viable_line_distance = %ld\n", previous_viable_line_distance);

                lseek(input_fd, next_viable_line_distance, SEEK_CUR);
                debugf("Seeked forwards next_viable_line_distance = %ld\n", next_viable_line_distance);
            }
        }

        // Determine if the current offset corresponds to a data line, if not, seek to the next one
        while (true) {
            union {
                struct {
                    uint64_t distance_to_previous;
                    uint64_t distance_to_next;
                } lengths;
                uint8_t bytes[16];
            } _length_headers;
            if (read(input_fd, &_length_headers.bytes, 16) == 0) {
                throw std::runtime_error("Reached end of file unexpectedly");
            }
            debugf("distance_to_previous = %lu, distance_to_next = %lu\n",
                _length_headers.lengths.distance_to_previous, _length_headers.lengths.distance_to_next);

            // IF prev offset value is zero and this is not the first line, must be an invalid location
            if (_length_headers.lengths.distance_to_previous == 0 &&
                    initial_lookup_offset != (long)(first_line_offset + data_start_offset)) {
                // seek ahead to next viable line start
                long seek_distance = sparse_config.multiplication_factor * sparse_config.block_size;
                seek_distance -= 16; // Already read this many bytes
                lseek(input_fd, seek_distance, SEEK_CUR);
                debugf("Offset %ld was not a data line, seeked to next viable offset %ld\n",
                    tellfd(input_fd) + 16 - seek_distance, tellfd(input_fd));
            } else {
                debugf("Offset was a data location, begin linear traversal\n");
                lseek(input_fd, -16, SEEK_CUR);
                break;
            }
        }

        debugf("Determined actual start offset for data in the query range: %ld\n", tellfd(input_fd));

        // max offset of the beginning of the last matching variant line
        // size_t end_variant_offset = sparse_config.compute_sparse_offset(query.get_reference_name(), query.get_end_position());

        string_t linebuf;
        string_reserve(&linebuf, 4 * 1024);
        size_t linelength;
        debugf("Starting linear variant enumeration from reference = %s %lu to %lu\n",
                query.get_reference_name().c_str(),
                query.get_start_position(),
                query.get_end_position());

        while (true) {
            // Important
            string_clear(&linebuf);

            uint64_t distance_to_previous, distance_to_next;
            size_t line_start_offset = lseek(input_fd, 0, SEEK_CUR);
            if (read(input_fd, &distance_to_previous, sizeof(uint64_t)) <= 0) {
                throw std::runtime_error("couldn't read from file");
            }
            if (read(input_fd, &distance_to_next, sizeof(uint64_t)) <= 0) {
                throw std::runtime_error("couldn't read from file");
            }

            debugf("distance_to_previous = %lu, distance_to_next = %lu\n", distance_to_previous, distance_to_next);

            // Not positioned at a data line
            if (distance_to_previous == 0 && distance_to_next == 0) {
                throw std::runtime_error("No previous or next distance values");
            }
            bool end_of_reference = false;
            if (distance_to_next == 0) {
                end_of_reference = true;
            }

            #ifdef DEBUG
            start = std::chrono::steady_clock::now();
            #endif
            debugf("current offset: %ld\n", lseek(input_fd, 0, SEEK_CUR));
            decompress2_data_line_fd2(input_fd, schema, &linebuf, &linelength);
            #ifdef DEBUG
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            debugf("decompress2_data_line time: %lu\n", duration.count());
            #endif

            debugf("compressed bytes read: %lu\n", linelength);
            // Update distance_to_next based on bytes already read from input stream
            //distance_to_next += linelength + 16; // line length plus 2 uint64s at start
            long bytes_read_so_far = (tellfd(input_fd) - line_start_offset);
            distance_to_next -= bytes_read_so_far;
            debugf("bytes_read_so_far = %lu, new distance_to_next = %lu\n", bytes_read_so_far, distance_to_next);

            // time copy of linebuf
            #ifdef DEBUG
            start = std::chrono::steady_clock::now();
            #endif
            SplitIterator spi(linebuf.buf, "\t");
            #ifdef DEBUG
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            debugf("SplitIterator.constructor time: %lu\n", duration.count());
            #endif

            #ifdef DEBUG
            start = std::chrono::steady_clock::now();
            #endif
            std::string reference_name = spi.next();
            std::string pos_str = spi.next();
            #ifdef DEBUG
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            debugf("SplitIterator.next time: %lu\n", duration.count());
            #endif

            #ifdef DEBUG
            start = std::chrono::steady_clock::now();
            #endif
            char *endptr = NULL;
            size_t pos = std::strtoul(pos_str.c_str(), &endptr, 10);
            if (endptr != pos_str.data() + pos_str.size()) {
                throw new std::runtime_error("Couldn't parse pos column: " + pos_str);
            }
            #ifdef DEBUG
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            debugf("strtoul time: %lu\n", duration.count());
            #endif

            debugf("line reference_name = %s, pos = %lu; query reference_name = %s, end_position = %lu\n",
                    reference_name.c_str(), pos,
                    query.get_reference_name().c_str(), query.get_end_position());

            if (reference_name == query.get_reference_name() && pos <= query.get_end_position()) {
                // Meets filter criteria, print the line
                fwrite(linebuf.buf, sizeof(char), linebuf.size, stdout); // newline included already

                if (end_of_reference) {
                    debugf("Reached end of reference %s\n", query.get_reference_name().c_str());
                    break;
                } else if (pos >= query.get_end_position()) {
                    debugf("Reached end of query range %lu\n", query.get_end_position());
                    break;
                } else {
                    debugf("Seeking ahead to next line\n");
                    #ifdef DEBUG
                    start = std::chrono::steady_clock::now();
                    long current_offset = tellfd(input_fd);
                    #endif
                    lseek(input_fd, distance_to_next, SEEK_CUR);
                    #ifdef DEBUG
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    debugf("lseek time: %lu\n", duration.count());
                    #endif

                    debugf("Previously at address: %ld, now at address: %ld\n", current_offset, tellfd(input_fd));
                }
            } else {
                break;
            }
        }

    }
    // No filter
    else {
        debugf("No filter criteria\n");
        throw std::runtime_error("sparse query with no filter is not yet implemented\n");
        uint8_t first_skip_bytes[8];
        if (read(input_fd, &first_skip_bytes, 8*sizeof(uint8_t)) <= 0) {
            throw std::runtime_error("Couldn't read from file");
        }
        uint64_t first_skip_count;
        uint8_array_to_uint64(first_skip_bytes, &first_skip_count);
        debugf("first_skip_count: %lu\n", first_skip_count);
    }
    close(input_fd);
}





class VcfPackedBinningIndexConfiguration {
public:
    VcfPackedBinningIndexConfiguration(int entries_per_bin):
            entries_per_bin(entries_per_bin) {

    }

    int entries_per_bin;

};

void create_binning_index(const std::string& compressed_input_filename, const std::string& index_filename) {
    int input_fd = open(compressed_input_filename.c_str(), O_RDONLY);
    if (input_fd < 0) {
        perror("open");
        throw std::runtime_error("Failed to open file: " + compressed_input_filename);
    }
    int status;
    VcfCompressionSchema schema;
}


void query_compressed_file(const std::string& input_filename, VcfCoordinateQuery query) {
    debugf("Querying %s for %s:%lu-%lu\n", input_filename.c_str(),
        query.get_reference_name().c_str(),
        query.get_start_position(),
        query.get_end_position());
    int input_fd = open(input_filename.c_str(), O_RDONLY);
    int status;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);

    size_t matched_line_count = 0;
    // std::string variant_line;
    // variant_line.reserve(1024 * 1024); // 1MiB
    string_t variant_line;
    string_reserve(&variant_line, 1024 * 1024); // 1 MiB

    while (true) {
        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                tellfd(input_fd),
                tellfd(input_fd));

        // right now all length headers are 4 bytes
        // TODO update to interpret variable-length length headers
        uint8_t line_length_header_bytes[4] = {0,0,0,0};
        status = read(input_fd, &line_length_header_bytes, 4);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished querying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("line_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                line_length_header_bytes[0],
                line_length_header_bytes[1],
                line_length_header_bytes[2],
                line_length_header_bytes[3]);

        uint8_t required_columns_length_header_bytes[4] = {0,0,0,0};
        status = read(input_fd, &required_columns_length_header_bytes, 4);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                required_columns_length_header_bytes[0],
                required_columns_length_header_bytes[1],
                required_columns_length_header_bytes[2],
                required_columns_length_header_bytes[3]);

        int64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                tellfd(input_fd),
                tellfd(input_fd));

        char c;
        std::string ref;
        while (true) {
            status = read(input_fd, &c, 1);
            if (status < 0) {
                throw VcfValidationError("Failed to read from file");
            } else if (status == 0) {
                throw VcfValidationError("Unexpected EOF while reading for reference name!");
            }

            read_bytes++;
            if (c == '\t') {
                break;
            }
            ref.push_back(c);
        }

        std::string pos_str;
        while (true) {
            status = read(input_fd, &c, 1);
            if (status < 0) {
                throw VcfValidationError("Failed to read from file");
            } else if (status == 0) {
                throw VcfValidationError("Unexpected EOF while reading for position!");
            }

            read_bytes++;
            if (c == '\t') {
                break;
            }
            pos_str.push_back(c);
        }

        bool conversion_success = false;
        uint64_t pos = str_to_uint64(pos_str, conversion_success);
        if (!conversion_success) {
            throw VcfValidationError(string_format("Malformed line, pos column must be int but got: %s", pos_str.c_str()).c_str());
        }
        debugf("read reference_name = %s, pos = %s from compressed line\n", ref.c_str(), pos_str.c_str());
        bool matches = query.matches(ref, pos);

        if (matches) {
            // TODO seekg backwards by appropriate number of bytes
            long seek_bytes = -1 * read_bytes;
            debugf("Line matches, so seeking %ld bytes\n", seek_bytes);
            //input_fstream.seekg(seek_bytes, input_fstream.cur);
            lseek(input_fd, seek_bytes, SEEK_CUR);

            debugf("Now positioned so next byte is at position %ld (0x%08lx)\n",
                    tellfd(input_fd), tellfd(input_fd));
            // variant_line.clear();
            string_clear(&variant_line);
            size_t compressed_line_length = 0;

            int status = decompress2_data_line_fd2(input_fd, schema, &variant_line, &compressed_line_length);
            if (status < 0) {
                break;
            }

            matched_line_count++;
            std::cout << variant_line.buf; // newline is included in decompress2_data_line

        } else {
            debugf("Line reference_name = %s, pos = %lu did not match\n", ref.c_str(), pos);
            // TODO interpret line_length_bytes to skip to next line
            LineLengthHeader line_length_header;
            line_length_header.set_extension_count(3); // TODO interpret
            line_length_header.deserialize(line_length_header_bytes);
            uint32_t line_length = line_length_header.length;
            // length header is not included in the "line length" header value itself
            // so subtract 4 from read_bytes
            // TODO update to dynamically account for size of header
            uint32_t skip_count = line_length - (read_bytes - 4);
            debugf("line length = %u, already read = %ld, so moving %u bytes from position 0x%08lx to next line\n",
                line_length, read_bytes - 4, skip_count, tellfd(input_fd));
            lseek(input_fd, skip_count, SEEK_CUR);
        }


    } // end line loop
    debugf("matched_line_count: %lu\n", matched_line_count);
    close(input_fd);
}

void gap_analysis(const std::string& input_filename) {
    int input_fd = open(input_filename.c_str(), O_RDONLY);
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    //const char newline = '\n';
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);

    size_t variant_line_count = 0;
    string_t variant_line;
    string_reserve(&variant_line, 1024 * 1024); // 1 MiB

    std::string start_position_filename("start-positions.txt");
    std::ofstream start_position_fstream(start_position_filename);

    unsigned char c;
    while (true) {
        if (peekfd(input_fd, &c) == 0) {
            // done
            debugf("Finished decompressing lines");
            break;
        }
        variant_line_count++;
        string_clear(&variant_line);
        size_t compressed_line_length = 0;
        int status = decompress2_data_line_fd2(input_fd, schema, &variant_line, &compressed_line_length);
        if (status < 0) {
            break;
        }
        // split the line
        std::vector<std::string> line_terms = split_string(std::string(variant_line.buf), "\t");
        start_position_fstream << line_terms[1];
        start_position_fstream << " ";
        start_position_fstream << std::to_string(variant_line.size);
        start_position_fstream << " ";
        start_position_fstream << std::to_string(compressed_line_length);
        start_position_fstream << "\n";

    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    start_position_fstream.flush();
    close(input_fd);
}


int main(int argc, char **argv) {
    // mostly just using posix fd functions, so disable sync
    std::ios_base::sync_with_stdio(false);

    int status;
    std::string action(argv[1]);

    if (action == "gap-analysis") {
        std::string input_filename(argv[2]);
        gap_analysis(input_filename);
    } else if (action == "compress" || action == "decompress") {
        std::string input_filename(argv[2]);
        std::string output_filename(argv[3]);
        if (input_filename == output_filename) {
            throw std::runtime_error("input and output file are the same");
        }
        if (action == "compress") {
            status = compress(input_filename, output_filename);
        } else {
            status = decompress2_fd(input_filename, output_filename);
        }

        if (status < 0) {
            std::cerr << "Error in compression of file" << std::endl;
        }
    } else if (action == "query") {
        std::string input_filename(argv[2]);
        std::string reference_name(argv[3]);
        std::string start_position_str(argv[4]);
        std::string end_position_str(argv[5]);
        bool conversion_success = false;
        uint64_t start_position = str_to_uint64(start_position_str, conversion_success);
        if (!conversion_success) {
            std::cerr << "Start position must be an integer: " << start_position_str << std::endl;
            return 1;
        }
        uint64_t end_position = str_to_uint64(end_position_str, conversion_success);
        if (!conversion_success) {
            std::cerr << "End position must be an integer: " << end_position_str << std::endl;
            return 1;
        }
        VcfCoordinateQuery query(reference_name, start_position, end_position);
        query_compressed_file(input_filename, query);
    } else if (action == "gap-analysis") {
        std::string input_filename(argv[2]);
        gap_analysis(input_filename);
    } else if (action == "sparsify") {
        if (argc < 4) {
            return usage();
        }
        std::string input_filename(argv[2]);
        std::string output_filename(argv[3]);
        sparsify_file_fd(input_filename, output_filename);
    } else if (action == "sparse-query") {
        std::string input_filename(argv[2]);
        std::string reference_name(argv[3]);
        std::string start_position_str(argv[4]);
        std::string end_position_str(argv[5]);
        bool conversion_success = false;
        uint64_t start_position = str_to_uint64(start_position_str, conversion_success);
        if (!conversion_success) {
            std::cerr << "Start position must be an integer: " << start_position_str << std::endl;
            return 1;
        }
        uint64_t end_position = str_to_uint64(end_position_str, conversion_success);
        if (!conversion_success) {
            std::cerr << "End position must be an integer: " << end_position_str << std::endl;
            return 1;
        }
        VcfCoordinateQuery query(reference_name, start_position, end_position);
        query_sparse_file_fd(input_filename, query);

    }
    else {
        std::cout << "Unknown action name: " << action << std::endl;
    }
}