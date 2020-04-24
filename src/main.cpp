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

#define SPARSE_EXTERNAL_INDEX_BLOCK_SIZE 256


int usage() {
    std::cerr << "./main [compress|decompress|sparsify] <input_file> <output_file>" << std::endl;
    return 1;
}

class VcfCoordinateQuery {
public:
    VcfCoordinateQuery(
        const std::string& reference_name,
        uint64_t start_position,
        uint64_t end_position) :
            reference_name(reference_name),
            start_position(start_position),
            end_position(end_position),
            has_start_position(true),
            has_end_position(true) {}

    VcfCoordinateQuery(
        const std::string& reference_name) :
            reference_name(reference_name),
            start_position(0),
            end_position(0),
            has_start_position(false),
            has_end_position(false) {}

    VcfCoordinateQuery() :
            reference_name(""),
            start_position(0),
            end_position(0),
            has_start_position(false),
            has_end_position(false) {}

    VcfCoordinateQuery& operator=(const VcfCoordinateQuery& rhs) {
        this->reference_name = rhs.reference_name;
        this->start_position = rhs.start_position;
        this->end_position = rhs.end_position;
        this->has_start_position = rhs.has_start_position;
        this->has_end_position = rhs.has_end_position;
        // this->ref_name_map = rhs.ref_name_map;
        return *this;
    }

    /**
     * Returns true if the input values are within the range of this query
     */
    bool matches(const std::string& reference_name, uint64_t position) {
        if (this->reference_name.size() > 0 && this->reference_name != reference_name) {
            return false;
        }
        if (this->has_start_position && position < this->start_position) {
            return false;
        }
        if (this->has_end_position && position > this->end_position) {
            return false;
        }
        return true;
    }

    int compare_to(const std::string& reference_name, uint64_t position) {
        uint8_t input_reference_name_idx = this->ref_name_map.reference_to_int(reference_name);
        uint8_t this_reference_name_idx = this->ref_name_map.reference_to_int(this->reference_name);

        // This is greater than arguments
        if (input_reference_name_idx < this_reference_name_idx
                || (input_reference_name_idx == this_reference_name_idx
                        && position < this->start_position)) {
            return 1;
        }
        // This is less than arguments
        else if (input_reference_name_idx > this_reference_name_idx
                || (input_reference_name_idx == this_reference_name_idx
                        && position > this->end_position)) {
            return -1;
        }
        // This is equal to (within range of) arguments
        else {
            return 0;
        }
    }

    int compare_to_range(const std::string& reference_name, uint64_t start_input, uint64_t end_input) {
        uint8_t input_reference_name_idx = this->ref_name_map.reference_to_int(reference_name);
        uint8_t this_reference_name_idx = this->ref_name_map.reference_to_int(this->reference_name);
        int ret = 0;

        // This is greater than arguments
        if (input_reference_name_idx < this_reference_name_idx
                || (input_reference_name_idx == this_reference_name_idx
                        && end_input < this->start_position)) {
            ret = 1;
        }
        // This is less than arguments
        else if (input_reference_name_idx > this_reference_name_idx
                || (input_reference_name_idx == this_reference_name_idx
                        && start_input > this->end_position)) {
            ret = -1;
        }
        // This is equal to (within range of) arguments
        else {
            ret = 0;
        }
        debugf("%s query(%s, %lu, %lu) - input(%s, %lu, %lu): %d\n",
            __FUNCTION__,
            this->reference_name.c_str(), this->start_position, this->end_position,
            reference_name.c_str(), start_input, end_input,
            ret);
        return ret;
    }


    // TODO no has_*_position checks here
    // bool after(const std::string& reference_name, uint64_t position) {
    //     uint8_t input_reference_name_idx = this->ref_name_map.reference_to_int(reference_name);
    //     uint8_t this_reference_name_idx = this->ref_name_map.reference_to_int(reference_name);

    //     if (input_reference_name_idx > this_reference_name_idx
    //             || (input_reference_name_idx == this_reference_name_idx
    //                     && position >= this->end_position)) {
    //         return true;
    //     }
    //     return false;
    // }

    bool has_criteria() {
        return this->reference_name.size() > 0
            || this->has_start_position
            || this->has_end_position;
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
    uint64_t end_position;
    bool has_start_position;
    bool has_end_position;
    reference_name_map ref_name_map;
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
    off_t lseek_ret = 0;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);

    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    start = std::chrono::steady_clock::now();
    #endif
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING decompress2_metadata_headers: %lu\n", duration.count());
    #endif

    // Leave default sparse config
    SparsificationConfiguration sparse_config;

    long off = tellfd(input_fd);
    if (off < 0) {
        throw std::runtime_error("ftell failed: " + std::to_string(off));
    }
    long data_start_offset = off + 8;
    uint64_t first_line_offset = 0;
    debugf("Reading first line offset value from file offset %ld\n", tellfd(input_fd));
    if (read(input_fd, &first_line_offset, sizeof(uint64_t)) < (int)sizeof(uint64_t)) {
        throw std::runtime_error("Failed to read first_line_offset value from file");
    }

    debugf("data_start_offset = %lu\n", data_start_offset);
    debugf("first_line_offset = %lu\n", first_line_offset);

    // Single variant lookup
    if (query.has_criteria() && (query.get_start_position() == query.get_end_position())) {
        debugf("Single variant lookup\n");
        size_t variant_offset = sparse_config.compute_sparse_offset(query.get_reference_name(), query.get_start_position());

        off64_t new_offset = data_start_offset + variant_offset;
        debugf("variant_offset = %lu, file_offset = %lu\n", variant_offset, new_offset);

        loff_t initial_lookup_offset = lseek64(input_fd, new_offset, SEEK_SET);
        debugf("initial_lookup_offset = %ld\n", initial_lookup_offset);
        if (initial_lookup_offset != new_offset) {
            perror("lseek");
            debugf("Failed to seek to line in file\n");
            return;
        }

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
            std::string linebuf;
            linebuf.reserve(4 * 1024); // 4 KiB
            // string_t linebuf;
            // string_reserve(&linebuf, 4 * 1024);
            size_t linelength;
            int status = decompress2_data_line_FILEwrapper(input_fd, schema, linebuf, &linelength);
            if (status == 0) {
                throw std::runtime_error("Unexpected EOF\n");
            } else if (status < 0) {
                throw std::runtime_error("Failed to decompress data line\n");
            }
            for (size_t i = 0; i < linebuf.size(); i++) {
                printf("%c", linebuf[i]);
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
                off_t current_offset = tellfd(input_fd);
                lseek_ret = lseek(input_fd, next_viable_line_distance, SEEK_CUR);
                if (lseek_ret != next_viable_line_distance + current_offset) {
                    perror("lseek");
                    debugf("Could not seek to next viable line address %ld, got %ld\n",
                        next_viable_line_distance + current_offset, lseek_ret);
                    return;
                }
                debugf("Seeked forwards next_viable_line_distance = %ld to %ld\n", next_viable_line_distance, tellfd(input_fd));
            }
        }

        // Determine if the current offset corresponds to a data line, if not, seek to the next one
        while (true) {
            // union {
            //     struct {
            //         uint64_t distance_to_previous;
            //         uint64_t distance_to_next;
            //     } lengths;
            //     uint8_t bytes[16];
            // } _length_headers;
            // if (read(input_fd, &_length_headers.bytes, 16) == 0) {
            //     throw std::runtime_error("Reached end of file unexpectedly");
            // }
            // debugf("distance_to_previous = %lu, distance_to_next = %lu\n",
            //     _length_headers.lengths.distance_to_previous, _length_headers.lengths.distance_to_next);
            uint8_t distance_headers[16];
            memset(distance_headers, 0, 16);

            if (read(input_fd, distance_headers, 16) < 16) {
                throw std::runtime_error("Reached end of file unexpectedly when reading distance headers");
            }
            uint64_t distance_to_previous = 0, distance_to_next = 0;
            uint8_array_to_uint64(distance_headers, &distance_to_previous);
            uint8_array_to_uint64(distance_headers + 8, &distance_to_next);
            debugf("distance_to_previous = %lu (0x%08lx), distance_to_next = %lu (0x%08lx)\n",
                distance_to_previous, distance_to_previous, distance_to_next, distance_to_next);

            // IF prev offset value is zero and this is not the first line, must be an invalid location
            if (distance_to_previous == 0 &&
                    initial_lookup_offset != (long)(first_line_offset + data_start_offset)) {
                // seek ahead to next viable line start
                long seek_distance = sparse_config.multiplication_factor * sparse_config.block_size;
                seek_distance -= 16; // Already read this many bytes
                lseek64(input_fd, seek_distance, SEEK_CUR);
                debugf("Offset %ld was not a data line, seeked to next viable offset %ld\n",
                    tellfd(input_fd) + 16 - seek_distance, // same as initial_lookup_offset
                    tellfd(input_fd));
            } else {
                debugf("Offset was a data location, begin linear traversal\n");
                lseek64(input_fd, -16, SEEK_CUR);
                break;
            }
        }

        debugf("Determined actual start offset for data in the query range: %ld\n", tellfd(input_fd));

        // max offset of the beginning of the last matching variant line
        // size_t end_variant_offset = sparse_config.compute_sparse_offset(query.get_reference_name(), query.get_end_position());

        std::string linebuf;
        linebuf.reserve(16 * 1024);
        size_t linelength;
        debugf("Starting linear variant enumeration from reference = %s %lu to %lu\n",
                query.get_reference_name().c_str(),
                query.get_start_position(),
                query.get_end_position());

        while (true) {
            // Important
            linebuf.clear();
            size_t line_start_offset = tellfd(input_fd);
            debugf("line_start_offset = %lu\n", line_start_offset);

            // uint64_t distance_to_previous, distance_to_next;
            // if (read(input_fd, &distance_to_previous, sizeof(uint64_t)) <= 0) {
            //     throw std::runtime_error("couldn't read from file");
            // }
            // if (read(input_fd, &distance_to_next, sizeof(uint64_t)) <= 0) {
            //     throw std::runtime_error("couldn't read from file");
            // }
            // debugf("distance_to_previous = %lu, distance_to_next = %lu\n", distance_to_previous, distance_to_next);

            uint8_t distance_headers[16];
            memset(distance_headers, 0, 16);

            if (read(input_fd, distance_headers, 16) < 16) {
                throw std::runtime_error("Reached end of file unexpectedly when reading distance headers");
            }
            uint64_t distance_to_previous = 0, distance_to_next = 0;
            uint8_array_to_uint64(distance_headers, &distance_to_previous);
            uint8_array_to_uint64(distance_headers + 8, &distance_to_next);
            debugf("distance_to_previous = %lu (0x%08lx), distance_to_next = %lu (0x%08lx)\n",
                distance_to_previous, distance_to_previous, distance_to_next, distance_to_next);

            // Not positioned at a data line
            if (distance_to_previous == 0 && distance_to_next == 0) {
                throw std::runtime_error("No previous or next distance values");
            }
            bool end_of_reference = false;
            if (distance_to_next == 0) {
                end_of_reference = true;
            }

            #ifdef TIMING
            start = std::chrono::steady_clock::now();
            #endif
            debugf("current offset: %ld\n", tellfd(input_fd));
            int status = decompress2_data_line_FILEwrapper(input_fd, schema, linebuf, &linelength);
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("TIMING decompress2_data_line: %lu\n", duration.count());
            #endif
            if (status == 0) {
                throw std::runtime_error("Unexpected EOF");
            } else if (status < 0) {
                throw std::runtime_error("Failed to decompress data line\n");
            }

            debugf("compressed bytes read: %lu\n", linelength);
            // Update distance_to_next based on bytes already read from input stream
            //distance_to_next += linelength + 16; // line length plus 2 uint64s at start
            long bytes_read_so_far = (tellfd(input_fd) - line_start_offset);
            distance_to_next -= bytes_read_so_far;
            debugf("bytes_read_so_far = %lu, new distance_to_next = %lu\n", bytes_read_so_far, distance_to_next);

            // time copy of linebuf
            #ifdef TIMING
            start = std::chrono::steady_clock::now();
            #endif
            SplitIterator spi(linebuf, "\t");
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("TIMING SplitIterator.constructor: %lu\n", duration.count());
            #endif

            #ifdef TIMING
            start = std::chrono::steady_clock::now();
            #endif
            std::string reference_name = spi.next();
            std::string pos_str = spi.next();
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("TIMING SplitIterator.next: %lu\n", duration.count());
            #endif

            #ifdef TIMING
            start = std::chrono::steady_clock::now();
            #endif
            char *endptr = NULL;
            size_t pos = std::strtoul(pos_str.c_str(), &endptr, 10);
            if (endptr != pos_str.data() + pos_str.size()) {
                throw new std::runtime_error("Couldn't parse pos column: " + pos_str);
            }
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("TIMING strtoul: %lu\n", duration.count());
            #endif

            debugf("line reference_name = %s, pos = %lu; query reference_name = %s, end_position = %lu\n",
                    reference_name.c_str(), pos,
                    query.get_reference_name().c_str(), query.get_end_position());

            if (reference_name == query.get_reference_name() && pos <= query.get_end_position()) {
                // Meets filter criteria, print the line
                fwrite(linebuf.c_str(), sizeof(char), linebuf.size(), stdout); // newline included already

                if (end_of_reference) {
                    debugf("Reached end of reference %s\n", query.get_reference_name().c_str());
                    break;
                } else if (pos >= query.get_end_position()) {
                    debugf("Reached end of query range %lu\n", query.get_end_position());
                    break;
                } else {
                    debugf("Seeking ahead to next line\n");
                    #ifdef DEBUG
                    long current_offset = tellfd(input_fd);
                    #endif

                    #ifdef TIMING
                    start = std::chrono::steady_clock::now();
                    #endif
                    lseek64(input_fd, distance_to_next, SEEK_CUR);
                    #ifdef TIMING
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING lseek: %lu\n", duration.count());
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


size_t struct_index_entry_size =
    sizeof(uint8_t) + sizeof(uint32_t) + sizeof(uint64_t);

struct index_entry {
    uint8_t reference_name_idx;
    uint32_t position;
    uint64_t byte_offset;
};

int write_index_entry_fd(int fd, struct index_entry *entry) {
    write(fd, &entry->reference_name_idx, sizeof(entry->reference_name_idx));
    write(fd, &entry->position, sizeof(entry->position));
    write(fd, &entry->byte_offset, sizeof(entry->byte_offset));
    return 0;
}

long write_index_entry(FILE *file, struct index_entry *entry) {
    // write(fd, &entry->reference_name_idx, sizeof(entry->reference_name_idx));
    // write(fd, &entry->position, sizeof(entry->position));
    // write(fd, &entry->byte_offset, sizeof(entry->byte_offset));

    long total_bytes = 0;
    total_bytes += fwrite(&entry->reference_name_idx, sizeof(entry->reference_name_idx), 1, file);
    total_bytes += fwrite(&entry->position, sizeof(entry->position), 1, file);
    total_bytes += fwrite(&entry->byte_offset, sizeof(entry->byte_offset), 1, file);
    return total_bytes;
}

/**
 * Returns 0 on success. Sets `bytes_read` to the number of bytes read. Should be the
 * sum of the size of the fields in a struct index_entry.
 */
int read_index_entry_fd(int fd, struct index_entry *entry, int *bytes_read) {
    int read_count = 0, n = 0;
    n = read(fd, &(entry->reference_name_idx), sizeof(uint8_t));
    if (n == 0) {
        debugf("Unexpected EOF\n");
        return -1;
    } else if (n < 0) {
        perror("read");
        debugf("Failed to read reference_name_idx, status = %d\n", n);
        return -1;
    } else {
        read_count += n;
    }

    n = read(fd, &entry->position, sizeof(uint32_t));
    if (n == 0) {
        debugf("Unexpected EOF\n");
        return -1;
    } else if (n < 0) {
        perror("read");
        debugf("Failed to read position, status = %d\n", n);
        return -1;
    } else {
        read_count += n;
    }

    n = read(fd, &entry->byte_offset, sizeof(uint64_t));
    if (n == 0) {
        debugf("Unexpected EOF\n");
        return -1;
    } else if (n < 0) {
        perror("read");
        debugf("Failed to read byte_offset, status = %d\n", n);
        return -1;
    } else {
        read_count += n;
    }

    *bytes_read = read_count;
    return 0;
}

/**
 * Returns 0 on success. Sets `bytes_read` to the number of bytes read. Should be the
 * sum of the size of the fields in a struct index_entry.
 */
int read_index_entry(FILE *file, struct index_entry *entry, int *bytes_read) {
    int read_count = 0, n = 0;
    // n = read(fd, &(entry->reference_name_idx), sizeof(uint8_t));
    n = fread(&(entry->reference_name_idx), sizeof(uint8_t), 1, file);
    if (n < 1) {
        debugf("Failed to read reference_name_idx, status = %d\n", n);
        debugf("Unexpected EOF\n");
        perror("fread");
        return -1;
    } else {
        read_count += n * sizeof(uint8_t);
    }

    n = fread(&entry->position, sizeof(uint32_t), 1, file);
    if (n < 1) {
        debugf("Failed to read position, status = %d\n", n);
        debugf("Unexpected EOF\n");
        perror("fread");
        return -1;
    } else {
        read_count += n * sizeof(uint32_t);
    }

    n = fread(&entry->byte_offset, sizeof(uint64_t), 1, file);
    if (n < 1) {
        debugf("Failed to read byte_offset, status = %d\n", n);
        debugf("Unexpected EOF\n");
        perror("read");
        return -1;
    } else {
        read_count += n * sizeof(uint64_t);
    }

    *bytes_read = read_count;
    return 0;
}


void create_sparse_binning_index(
        const std::string& compressed_input_filename,
        const std::string& index_filename,
        SparsificationConfiguration& sparse_config) {
    FILE *input_file = fopen(compressed_input_filename.c_str(), "r");
    if (input_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open file: " + compressed_input_filename);
    }
    int output_fd = open(index_filename.c_str(), DEFAULT_FILE_CREATE_FLAGS, DEFAULT_FILE_CREATE_MODE);
    if (output_fd < 0) {
        perror("open");
        throw std::runtime_error("Failed to open output file: " + index_filename);
    }

    // write one byte at beginning of file
    // char c = 'A';
    // fwrite(&c, 1, 1, output_file);

    int status;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(input_file, meta_header_lines, schema);

    reference_name_map ref_name_map;

    std::vector<byte_t> line_bytes;
    line_bytes.reserve(16 * 1024);

    // Counter to handle entries per bin
    size_t line_number = 0;

    // Iterate through lines of compressed file
    while (true) {
        long line_byte_offset = ftell(input_file);

        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                line_byte_offset,
                line_byte_offset);

        compressed_line_length_headers line_length_headers;
        memset(&line_length_headers, 0, sizeof(compressed_line_length_headers));
        status = read_compressed_line_length_headers(input_file, &line_length_headers);
        if (status == 0) {
            debugf("Finished creating index\n");
            break;
        } else if (status != (int) compressed_line_length_headers_size) {
            throw std::runtime_error("Failed to read line length headers");
        }

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                ftell(input_file),
                ftell(input_file));

        uint8_t line_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.line_length, line_length_header_bytes);
        uint8_t required_columns_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.required_columns_length, required_columns_length_header_bytes);

        debugf("Line length: %d\n", line_length_headers.line_length);

        // collect bytes to end of line
        size_t i = 0;
        line_bytes.clear();
        if (line_bytes.capacity() < line_length_headers.line_length + read_bytes) {
            line_bytes.reserve(line_length_headers.line_length + read_bytes);
        }

        // two uint64 placeholders for diffs to previous, next line
        for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
            line_bytes.push_back(0);
        }

        line_bytes.push_back(line_length_header_bytes[0]);
        line_bytes.push_back(line_length_header_bytes[1]);
        line_bytes.push_back(line_length_header_bytes[2]);
        line_bytes.push_back(line_length_header_bytes[3]);
        line_bytes.push_back(required_columns_length_header_bytes[0]);
        line_bytes.push_back(required_columns_length_header_bytes[1]);
        line_bytes.push_back(required_columns_length_header_bytes[2]);
        line_bytes.push_back(required_columns_length_header_bytes[3]);

        debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        bool got_reference_name = false;
        std::string reference_name;
        reference_name.reserve(32);

        bool got_pos = false;
        std::string pos_str;
        pos_str.reserve(32);

        size_t pos = 0;

        // iterate over the rest of the line
        // subtract 4 due to length of required_columns_length_header_bytes, which are already read
        while (i++ < line_length_headers.line_length - 4) {
            unsigned char b;
            if (fread(&b, 1, 1, input_file) < 1) {
                std::string msg = string_format(
                    "Unexpectedly reached end of compressed file, line header said %d, but only read %d bytes from line",
                    line_length_headers.line_length, i);
                throw VcfValidationError(msg.c_str());
            }
            line_bytes.push_back(b);

            if (!got_reference_name) {
                if (b != '\t') {
                    reference_name.push_back(b);
                } else {
                    if (reference_name.size() == 0) {
                        throw std::runtime_error("Line did not contain a reference name");
                    } else {
                        debugf("Got reference name: %s\n", reference_name.c_str());
                        got_reference_name = true;
                    }
                }
            } else if (!got_pos) {
                if (b != '\t') {
                    pos_str.push_back(b);
                } else {
                    if (pos_str.size() == 0) {
                        throw std::runtime_error("Line did not contain a position value");
                    } else {
                        debugf("Got position: %s\n", pos_str.c_str());
                        got_pos = true;
                        char *endptr = NULL;
                        pos = strtoul(pos_str.c_str(), &endptr, 10);
                        if (endptr != pos_str.c_str() + pos_str.size()) {
                            throw std::runtime_error("Failed to parse full position value to long: " + pos_str);
                        }
                    }
                }
            }
        }

        size_t sparse_offset = sparse_config.compute_sparse_offset(reference_name, pos);
        off_t lseek_ret = lseek(output_fd, sparse_offset, SEEK_SET);
        if (lseek_ret != (off_t) sparse_offset) {
            perror("lseek");
            debugf("Failed to seek to offset %ld, got %ld\n", sparse_offset, lseek_ret);
            return;
        }
        struct index_entry entry;
        entry.reference_name_idx = ref_name_map.reference_to_int(reference_name);
        entry.position = pos;
        entry.byte_offset = line_byte_offset;
        debugf("Writing entry to index %d %d %ld\n",
            entry.reference_name_idx, entry.position, entry.byte_offset);
        write_index_entry_fd(output_fd, &entry);

        // Determine whether to write this entry to the index or whether it is inside the current bin
        // if (line_number % index_configuration.entries_per_bin == 0) {
        //     struct index_entry entry;
        //     // char *endptr;
        //     entry.reference_name_idx = ref_name_map.reference_to_int(reference_name.buf);
        //     entry.position = pos;
        //     entry.byte_offset = line_byte_offset;
        //     debugf("Writing entry to index %d %d %ld\n",
        //         entry.reference_name_idx, entry.position, entry.byte_offset);
        //     write_index_entry(output_file, &entry);
        // } else {
        //     debugf("Not writing entry to index\n");
        // }
        line_number++;
    }

    close(output_fd);
    fclose(input_file);
}


void query_sparse_binned_index(
        const std::string& compressed_filename,
        const std::string& index_filename,
        VcfCoordinateQuery& query,
        SparsificationConfiguration& sparse_config) {
    int status;
    off_t lseek_ret;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    debugf("Opening compressed file\n");
    FILE *compressed_file = fopen(compressed_filename.c_str(), "r");
    if (compressed_file == NULL) {
        perror("fopen");
        debugf("Failed to open input file: %s\n", compressed_filename.c_str());
        return;
    }
    debugf("Successfully opened file\n");
    if (!file_exists(index_filename.c_str())) {
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }
    debugf("Opening index file\n");
    int index_fd = open(index_filename.c_str(), O_RDONLY);
    if (index_fd < 0) {
        perror("open");
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }

    // test SEEK_DATA
    off_t test_seek_data_offset = 10000000;
    lseek_ret = lseek(index_fd, test_seek_data_offset, SEEK_DATA);
    if (lseek_ret < 0) {
        perror("fseek");
        debugf("SEEK_DATA failed: %ld\n", lseek_ret);
        return;
    }


    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING file_open: %lu\n", duration.count());
    #endif

    debugf("Parsing metadata lines and header line\n");
    VcfCompressionSchema schema;
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(compressed_file, meta_header_lines, schema);

    struct index_entry start_entry;
    memset(&start_entry, 0, sizeof(start_entry));

    reference_name_map ref_name_map;

    // uint8_t query_reference_name_idx = ref_name_map.reference_to_int(query.get_reference_name());

    size_t sparse_offset_s = sparse_config.compute_sparse_offset(
        query.get_reference_name(), query.get_start_position());
    long sparse_offset = (long) sparse_offset_s;
    long index_size = file_size(index_filename.c_str());
    long entry_count = index_size / struct_index_entry_size;
    debugf("Index of size %ld has %ld entries\n", index_size, entry_count);

    debugf("Seeking to offset %ld\n", sparse_offset);
    lseek_ret = lseek(index_fd, sparse_offset, SEEK_SET);
    if (lseek_ret != sparse_offset) {
        perror("lseek");
        debugf("Failed to seek to offset %ld, returned %ld\n", sparse_offset, lseek_ret);
        return;
    }

    struct index_entry entry;
    memset(&entry, 0, sizeof(struct index_entry));

    // Try to read an entry
    int bytes_read = 0;
    if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
        debugf("Failed to read index entry from file\n");
        fclose(compressed_file);
        close(index_fd);
        return;
    }

    // Not realistically possible given semantic constraints on values
    if (entry.reference_name_idx == 0 && entry.position == 0 && entry.byte_offset == 0) {
        debugf("Entry was empty, looking for next available entry\n");

        // Determine if in a hole
        long current_offset = tellfd(index_fd);
        debugf("Calling SEEK_DATA from offset %ld\n", current_offset);
        lseek_ret = lseek(index_fd, current_offset, SEEK_DATA);
        if (lseek_ret < 0) { // new_offset is same as lseek_ret
            perror("fseek");
            debugf("Failed to fseek SEEK_DATA, status = %ld\n", lseek_ret);
            return;
        }
        long new_offset = tellfd(index_fd);
        debugf("SEEK_DATA moved offset from %ld to %ld\n", current_offset, new_offset);
        if (new_offset != current_offset) {
            // was in a hole
            debugf("Was in a hole, now is not\n");
        }

        // Brute search ahead
        while (true) {
            debugf("Searching ahead for non empty index\n");
            long offset_before = tellfd(index_fd);
            if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
                debugf("Failed to read index entry from file\n");
                fclose(compressed_file);
                close(index_fd);
                return;
            }
            // long offset_after = tellfd(index_fd);

            if (entry.reference_name_idx == 0 && entry.position == 0 && entry.byte_offset == 0) {
                long next_viable_entry_address = offset_before + SPARSE_EXTERNAL_INDEX_BLOCK_SIZE;
                debugf("Entry still empty, seeking ahead %ld bytes to address %ld\n",
                    next_viable_entry_address - tellfd(index_fd), next_viable_entry_address);
                lseek_ret = lseek(index_fd, next_viable_entry_address, SEEK_SET);
                if (lseek_ret != next_viable_entry_address) {
                    perror("lseek");
                    debugf("Could not seek to next viable entry address\n");
                }
                continue;
            } else {
                debugf("Found non empty entry\n");
                break;
            }
        }
    }

    debugf("Found index entry reference_name = %u, position = %u, byte_offset = %lu\n",
        entry.reference_name_idx, entry.position, entry.byte_offset);

    std::string linebuf;
    linebuf.reserve(16 * 4096);


    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    enum query_state {
        BEFORE_QUERY, IN_QUERY, AFTER_QUERY
    };
    enum query_state current_state = query_state::BEFORE_QUERY;
    #endif

    if (entry_count > 0) {
        debugf("entry reference_name_idx = %u, position = %u, byte_offset = %lu\n",
            entry.reference_name_idx, entry.position, entry.byte_offset);

        // Iterate through data file
        fseek(compressed_file, entry.byte_offset, SEEK_SET);

        // Record the number of lines scanned before hitting the desired range
        int before_count = 0;

        while (true) {
            linebuf.clear();
            size_t compressed_line_length;
            status = decompress2_data_line(compressed_file, schema, linebuf, &compressed_line_length);
            if (status == 0) {
                // EOF
                debugf("End of input file\n");
                break;
            } else if (status < 0) {
                throw VcfValidationError("Error, failed to read next line from compressed input file");
            }
            SplitIterator spi(linebuf, "\t");
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string reference_name = spi.next();
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string pos_str = spi.next();
            bool conversion_success = false;
            uint64_t pos = str_to_uint64(pos_str, conversion_success);
            if (!conversion_success) {
                throw VcfValidationError(("Failed to parse integer pos from " + pos_str).c_str());
            }
            debugf("Checking reference_name = %s, pos = %lu against query reference_name = %s, start = %lu, end = %lu\n",
                reference_name.c_str(), pos, query.get_reference_name().c_str(), query.get_start_position(), query.get_end_position());

            int query_comparison = query.compare_to(reference_name, pos);
            if (query_comparison > 0) {
                // still before query
                debugf("Query did not match, state is still BEFORE_QUERY\n");
                before_count++;
                continue;
            } else if (query_comparison == 0) {
                debugf("State is in range of query\n");
                #ifdef TIMING
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }
                current_state = query_state::IN_QUERY;
                #endif
                printf(linebuf.c_str());
            } else if (query_comparison < 0) {
                debugf("State is after query\n");
                // after query
                #ifdef TIMING
                // If previous state was before, now not before, output decompress_seeking timer
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }

                current_state = query_state::AFTER_QUERY;
                #endif
                break;

            } else {
                throw std::runtime_error("Invalid state, query could not be compared");
            }
        }
        debugf("lines decompressed before query = %d\n", before_count);

    } else {
        debugf("Index was empty\n");
    }

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING decompress_iteration: %lu\n", duration.count());
    #endif

    close(index_fd);
    fclose(compressed_file);
}

long read_to(FILE *input_file, unsigned char end, bool remove_end, std::string& out) {
    unsigned char cur = 0;
    long counter = 0;
    while (fread(&cur, 1, 1, input_file) == 1) {
        counter++;
        if (cur != end) {
            out.push_back(reinterpret_cast<char&>(cur));
        } else {
            if (!remove_end) {
                off_t fseek_ret = fseek(input_file, -1, SEEK_CUR);
                if (fseek_ret != 0) {
                    return -1;
                }
                return counter - 1;
            }
            return counter;
        }
    }
    return -1;
}

int parse_kvp(const std::string& input, std::map<std::string,std::string>& output_map) {
    const std::vector<std::string> pairs = split_string(input, ";");
    for (size_t pair_i = 0; pair_i < pairs.size(); pair_i++) {
        const std::string& pair = pairs[pair_i];
        const std::vector<std::string> parts = split_string(pair, "=");
        if (parts.size() == 2) {
            const std::string& key = parts[0];
            const std::string& val = parts[1];
            // No duplicate checking
            output_map[key] = val;
        } else if (parts.size() == 1) {
            const std::string& key = parts[0];
            const std::string& val = "";
            // No duplicate checking
            output_map[key] = val;
        } else {
            throw std::runtime_error("Invalid kvp format: " + input);
        }
    }
    return 0;
}

bool alt_is_structural(const std::string& alt) {
    return alt.find('<') != std::string::npos;
}

long compute_end_position(
        long pos,
        const std::string& reference_name,
        const std::string& alt,
        const std::string& info) {
    long read_to_ret = 0;
    long end_position = 0;
    // debugf("ID=%s\n", id.c_str());
    // debugf("REF=%s\n", ref.c_str());
    debugf("ALT=%s\n", alt.c_str());

    if (alt_is_structural(alt)) {
        debugf("ALT is structural: %s\n", alt.c_str());
        // std::string info;
        debugf("INFO=%s\n", info.c_str());
        std::map<std::string,std::string> info_kvp;
        if (parse_kvp(info, info_kvp) < 0) {
            throw std::runtime_error("Failed to parse info kvp");
        }
        std::string svtype = info_kvp["SVTYPE"];
        debugf("Structural variant type: %s\n", svtype.c_str());

        // SV without END info: SVA, ALU (alt begins with: <INS:)
        if (info_kvp.count("END")) {
            long end;
            std::vector<std::string> end_strings = split_string(info_kvp["END"], ",");
            long max_end = 0;
            for (auto iter = end_strings.begin(); iter != end_strings.end(); iter++) {
                if (str_to_long(*iter, &end) != 0) {
                    throw std::runtime_error("Failed to parse END integer: " + *iter);
                }
                if (end > max_end) {
                    max_end = end;
                }
            }
            end_position = std::abs(max_end);
        } else if (info_kvp.count("SVLEN")) {
            long svlen = 0, max_svlen = 0;
            std::vector<std::string> svlen_strings = split_string(info_kvp["SVLEN"], ",");
            for (size_t i = 0; i < svlen_strings.size(); i++) {
                if (str_to_long(svlen_strings[i], &svlen) != 0) {
                    throw std::runtime_error("Failed to parse SVLEN integer: " + svlen_strings[i]);
                }
                // Using abs to handle INS types, and DEL, DUP together
                if (std::abs(svlen) > max_svlen) {
                    max_svlen = std::abs(svlen);
                }
            }
            end_position = pos + max_svlen - 1;
        } else {
            debugf("Could not find END or SVLEN to determine end position of structural variant, using start_position\n");
            end_position = pos;
        }

    } else {
        debugf("Non structural variant\n");
        std::vector<std::string> alts = split_string(alt, ",");
        size_t max_alt_size = 0;
        for (auto iter = alts.begin(); iter != alts.end(); iter++) {
            if (iter->size() > max_alt_size) {
                max_alt_size = iter->size();
            }
        }
        // POS plus longest region, either insertion or deletion
        // end_position = pos + std::abs((long)(reference_name.size() - max_alt_size));
        if (reference_name.size() >= max_alt_size) {
            end_position = pos + reference_name.size();
        } else {
            end_position = pos + max_alt_size - 1;
        }
    }
    return end_position;
}

void create_binned_index2(
        const std::string& compressed_input_filename,
        const std::string& index_filename,
        VcfPackedBinningIndexConfiguration& index_configuration) {
    FILE *input_file = fopen(compressed_input_filename.c_str(), "r");
    if (input_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open file: " + compressed_input_filename);
    }
    FILE *output_file = fopen(index_filename.c_str(), "w");
    if (output_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open output file: " + index_filename);
    }

    int status;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(input_file, meta_header_lines, schema);

    reference_name_map ref_name_map;

    std::vector<byte_t> line_bytes;
    line_bytes.reserve(16 * 1024);

    // Counter to handle entries per bin
    size_t line_number = 0;

    // Iterate through lines of compressed file
    // long total_write_bytes = 0;

    // keep map of binned position -> end_position, index entry
    std::map<size_t,std::pair<size_t,struct index_entry>> index_map;


    while (true) {
        long line_byte_offset = ftell(input_file);

        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                line_byte_offset,
                line_byte_offset);

        compressed_line_length_headers line_length_headers;
        memset(&line_length_headers, 0, sizeof(compressed_line_length_headers));
        status = read_compressed_line_length_headers(input_file, &line_length_headers);
        if (status == 0) {
            debugf("Finished creating index\n");
            break;
        } else if (status != (int) compressed_line_length_headers_size) {
            throw std::runtime_error("Failed to read line length headers");
        }

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                ftell(input_file),
                ftell(input_file));

        uint8_t line_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.line_length, line_length_header_bytes);
        uint8_t required_columns_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.required_columns_length, required_columns_length_header_bytes);
        debugf("Line length: %d\n", line_length_headers.line_length);

        // collect bytes to end of line
        // size_t i = 0;
        line_bytes.clear();
        // if (line_bytes.capacity() < line_length_headers.line_length + read_bytes) {
        //     line_bytes.reserve(line_length_headers.line_length + read_bytes);
        // }

        // // two uint64 placeholders for diffs to previous, next line
        // for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
        //     line_bytes.push_back(0);
        // }

        // line_bytes.push_back(line_length_header_bytes[0]);
        // line_bytes.push_back(line_length_header_bytes[1]);
        // line_bytes.push_back(line_length_header_bytes[2]);
        // line_bytes.push_back(line_length_header_bytes[3]);
        // line_bytes.push_back(required_columns_length_header_bytes[0]);
        // line_bytes.push_back(required_columns_length_header_bytes[1]);
        // line_bytes.push_back(required_columns_length_header_bytes[2]);
        // line_bytes.push_back(required_columns_length_header_bytes[3]);
        // debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        // bool got_reference_name = false;
        std::string reference_name;
        reference_name.reserve(32);

        // bool got_pos = false;
        std::string pos_str;
        pos_str.reserve(32);
        size_t pos = 0;
        size_t end_position = 0;

        long read_to_ret;
        read_to_ret = read_to(input_file, '\t', true, reference_name);
        if (read_to_ret <= 0) {
            throw std::runtime_error("Failed to read reference name");
        }
        read_to_ret = read_to(input_file, '\t', true, pos_str);
        if (read_to_ret <= 0) {
            throw std::runtime_error("Failed to read position");
        }
        bool success;
        pos = str_to_uint64(pos_str, success);
        if (!success) {
            throw std::runtime_error("Failed to parse pos: " + pos_str);
        }
        std::string id;
        read_to_ret = read_to(input_file, '\t', true, id);
        std::string ref;
        read_to_ret = read_to(input_file, '\t', true, ref);
        std::string alt;
        read_to_ret = read_to(input_file, '\t', true, alt);


        std::string info;
        std::string qual;
        std::string filter;
        if (read_to(input_file, '\t', true, qual) <= 0) {
            throw std::runtime_error("Failed to read qual");
        }
        if (read_to(input_file, '\t', true, filter) <= 0) {
            throw std::runtime_error("Failed to read filter");
        }
        if (read_to(input_file, '\t', true, info) <= 0) {
            throw std::runtime_error("Failed to read info");
        }
        debugf("CHR=%s\n", reference_name.c_str());
        debugf("POS=%s\n", pos_str.c_str());
        debugf("ID=%s\n", id.c_str());
        debugf("REF=%s\n", ref.c_str());
        debugf("ALT=%s\n", alt.c_str());
        debugf("QUAL=%s\n", qual.c_str());
        debugf("FILTER=%s\n", filter.c_str());
        debugf("INFO=%s\n", info.c_str());

        end_position = compute_end_position(pos, reference_name, alt, info);

        // if (alt_is_structural(alt)) {
        //     debugf("ALT is structural: %s\n", alt.c_str());
        //     std::string info;
        //     std::string qual;
        //     std::string filter;
        //     if (read_to(input_file, '\t', true, qual) <= 0) {
        //         throw std::runtime_error("Failed to read qual");
        //     }
        //     if (read_to(input_file, '\t', true, filter) <= 0) {
        //         throw std::runtime_error("Failed to read filter");
        //     }
        //     if (read_to(input_file, '\t', true, info) <= 0) {
        //         throw std::runtime_error("Failed to read info");
        //     }
        //     debugf("QUAL=%s\n", qual.c_str());
        //     debugf("FILTER=%s\n", filter.c_str());
        //     debugf("INFO=%s\n", info.c_str());
        //     std::map<std::string,std::string> info_kvp;
        //     if (parse_kvp(info, info_kvp) < 0) {
        //         throw std::runtime_error("Failed to parse info kvp");
        //     }
        //     std::string svtype = info_kvp["SVTYPE"];
        //     debugf("Structural variant type: %s\n", svtype.c_str());

        //     // SV without END info: SVA, ALU (alt begins with: <INS:)
        //     if (info_kvp.count("END")) {
        //         long end;
        //         std::vector<std::string> end_strings = split_string(info_kvp["END"], ",");
        //         long max_end = 0;
        //         for (auto iter = end_strings.begin(); iter != end_strings.end(); iter++) {
        //             if (str_to_long(*iter, &end) != 0) {
        //                 throw std::runtime_error("Failed to parse END integer: " + *iter);
        //             }
        //             if (end > max_end) {
        //                 max_end = end;
        //             }
        //         }
        //         end_position = std::abs(max_end);
        //     } else if (info_kvp.count("SVLEN")) {
        //         long svlen = 0, max_svlen = 0;
        //         std::vector<std::string> svlen_strings = split_string(info_kvp["SVLEN"], ",");
        //         for (size_t i = 0; i < svlen_strings.size(); i++) {
        //             if (str_to_long(svlen_strings[i], &svlen) != 0) {
        //                 throw std::runtime_error("Failed to parse SVLEN integer: " + svlen_strings[i]);
        //             }
        //             // Using abs to handle INS types, and DEL, DUP together
        //             if (std::abs(svlen) > max_svlen) {
        //                 max_svlen = std::abs(svlen);
        //             }
        //         }
        //         end_position = pos + max_svlen - 1;
        //     } else {
        //         debugf("Could not find END or SVLEN to determine end position of structural variant, using start_position\n");
        //         end_position = pos;
        //     }

        // } else {
        //     debugf("Non structural variant\n");
        //     std::vector<std::string> alts = split_string(alt, ",");
        //     size_t max_alt_size = 0;
        //     for (auto iter = alts.begin(); iter != alts.end(); iter++) {
        //         if (iter->size() > max_alt_size) {
        //             max_alt_size = iter->size();
        //         }
        //     }
        //     // POS plus longest region, either insertion or deletion
        //     // end_position = pos + std::abs((long)(reference_name.size() - max_alt_size));
        //     if (reference_name.size() >= max_alt_size) {
        //         end_position = pos + reference_name.size();
        //     } else {
        //         end_position = pos + max_alt_size;
        //     }
        // }

        // insert index entries for every position between [pos, end_position]
        // keep map of binned position -> index entry

        // insert or update index entries which exist if this for this positional range [pos, end_position]

        uint8_t reference_name_idx = ref_name_map.reference_to_int(reference_name);

        size_t map_max_start_position = 0;
        size_t map_max_end_position = 0;
        if (index_map.size() != 0) {
            const std::pair<size_t,std::pair<size_t,struct index_entry>> max_index_map_entry = *(index_map.rbegin());
            map_max_start_position = max_index_map_entry.first;
            map_max_end_position = max_index_map_entry.second.first;
            debugf("index current max position = %lu\n", map_max_start_position);

        }

        size_t first_index_position_to_write;
        if (pos >= map_max_start_position) {
            first_index_position_to_write = pos;
        } else {
            first_index_position_to_write = map_max_start_position;
        }

        // if (first_index_position_to_write <= end_position) {
            size_t count = (size_t) ((end_position - first_index_position_to_write) / index_configuration.entries_per_bin);
            count += 1; // allow trailing fragment for rounding in above calculation

            for (
                    size_t index_position = first_index_position_to_write;
                    index_position <= (end_position + count*index_configuration.entries_per_bin);
                    index_position += index_configuration.entries_per_bin) {
                struct index_entry new_entry;
                new_entry.reference_name_idx = reference_name_idx;
                new_entry.position = index_position; // starting position for a query
                new_entry.byte_offset = line_byte_offset;
                debugf("Inserting new index ref_idx=%d, entry start=%u, end=%lu, byte_offset=%lu\n",
                    new_entry.reference_name_idx, new_entry.position, end_position, new_entry.byte_offset);
                index_map[index_position] = std::pair<size_t,struct index_entry>(end_position, new_entry);
            }
        // }


        // Determine whether to write this entry to the index or whether it is inside the current bin
        // if (line_number % index_configuration.entries_per_bin == 0) {
        //     struct index_entry entry;
        //     // char *endptr;
        //     entry.reference_name_idx = ref_name_map.reference_to_int(reference_name);
        //     entry.position = pos;
        //     entry.byte_offset = line_byte_offset;
        //     debugf("Writing entry to index %d %d %ld\n",
        //         entry.reference_name_idx, entry.position, entry.byte_offset);

        //     long index_entry_bytes = write_index_entry(output_file, &entry);
        //     total_write_bytes += index_entry_bytes;
        //     debugf("Wrote %ld bytes to index, new size: %ld\n", index_entry_bytes, total_write_bytes);
        // } else {
        //     debugf("Not writing entry to index\n");
        // }
        line_number++;

        // TODO seek to next line
        long next_line_distance = line_length_headers.line_length - (ftell(input_file) - line_byte_offset);
        next_line_distance += 4; // include line_length header bytes itself
        debugf("Line length is %u, currently have read %ld bytes, seeking ahead %ld more\n",
            line_length_headers.line_length, ftell(input_file) - line_byte_offset, next_line_distance);
        off_t fseek_ret = fseek(input_file, next_line_distance, SEEK_CUR);
        if (fseek_ret != 0) {
            throw std::runtime_error(string_format(
                "Failed to seek forward %ld bytes to offset %ld\n",
                next_line_distance, ftell(input_file) + next_line_distance));
        }
    }

    for (auto iter = index_map.begin(); iter != index_map.end(); iter++) {
        struct index_entry entry = iter->second.second;
        write_index_entry(output_file, &entry);
    }

    fclose(output_file);
    fclose(input_file);
}


void create_binned_index(const std::string& compressed_input_filename, const std::string& index_filename, VcfPackedBinningIndexConfiguration& index_configuration) {
    FILE *input_file = fopen(compressed_input_filename.c_str(), "r");
    if (input_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open file: " + compressed_input_filename);
    }
    FILE *output_file = fopen(index_filename.c_str(), "w");
    if (output_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open output file: " + index_filename);
    }

    int status;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(input_file, meta_header_lines, schema);

    reference_name_map ref_name_map;

    std::vector<byte_t> line_bytes;
    line_bytes.reserve(16 * 1024);

    // Counter to handle entries per bin
    size_t line_number = 0;

    // Iterate through lines of compressed file
    long total_write_bytes = 0;
    while (true) {
        // long line_byte_offset = tellfd(input_fd);
        long line_byte_offset = ftell(input_file);

        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                line_byte_offset,
                line_byte_offset);

        compressed_line_length_headers line_length_headers;
        memset(&line_length_headers, 0, sizeof(compressed_line_length_headers));
        status = read_compressed_line_length_headers(input_file, &line_length_headers);
        if (status == 0) {
            debugf("Finished creating index\n");
            break;
        } else if (status != (int) compressed_line_length_headers_size) {
            throw std::runtime_error("Failed to read line length headers");
        }

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                ftell(input_file),
                ftell(input_file));

        uint8_t line_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.line_length, line_length_header_bytes);
        uint8_t required_columns_length_header_bytes[4];
        uint32_to_uint8_array(line_length_headers.required_columns_length, required_columns_length_header_bytes);


        debugf("Line length: %d\n", line_length_headers.line_length);

        // collect bytes to end of line
        size_t i = 0;
        line_bytes.clear();
        if (line_bytes.capacity() < line_length_headers.line_length + read_bytes) {
            line_bytes.reserve(line_length_headers.line_length + read_bytes);
        }

        // two uint64 placeholders for diffs to previous, next line
        for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
            line_bytes.push_back(0);
        }

        line_bytes.push_back(line_length_header_bytes[0]);
        line_bytes.push_back(line_length_header_bytes[1]);
        line_bytes.push_back(line_length_header_bytes[2]);
        line_bytes.push_back(line_length_header_bytes[3]);
        line_bytes.push_back(required_columns_length_header_bytes[0]);
        line_bytes.push_back(required_columns_length_header_bytes[1]);
        line_bytes.push_back(required_columns_length_header_bytes[2]);
        line_bytes.push_back(required_columns_length_header_bytes[3]);

        debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        bool got_reference_name = false;
        string_t reference_name;
        string_reserve(&reference_name, 32);

        bool got_pos = false;
        string_t pos_str;
        string_reserve(&pos_str, 32);

        size_t pos = 0;

        // iterate over the rest of the line
        // subtract 4 due to length of required_columns_length_header_bytes, which are already read
        while (i++ < line_length_headers.line_length - 4) {
            unsigned char b;
            // if (read(input_fd, &b, 1) <= 0) {
            if (fread(&b, 1, 1, input_file) < 1) {
                std::string msg = string_format(
                    "Unexpectedly reached end of compressed file, line header said %d, but only read %d bytes from line",
                    line_length_headers.line_length, i);
                throw VcfValidationError(msg.c_str());
            }
            line_bytes.push_back(b);

            if (!got_reference_name) {
                if (b != '\t') {
                    string_appendc(&reference_name, b);
                } else {
                    if (reference_name.size == 0) {
                        throw std::runtime_error("Line did not contain a reference name");
                    } else {
                        debugf("Got reference name: %s\n", reference_name.buf);
                        got_reference_name = true;
                    }
                }
            } else if (!got_pos) {
                if (b != '\t') {
                    // pos_str += b;
                    string_appendc(&pos_str, b);
                } else {
                    // if (pos_str.size() == 0) {
                    if (pos_str.size == 0) {
                        throw std::runtime_error("Line did not contain a position value");
                    } else {
                        debugf("Got position: %s\n", pos_str.buf);
                        got_pos = true;
                        char *endptr = NULL;
                        pos = strtoul(pos_str.buf, &endptr, 10);
                        if (endptr != pos_str.buf + pos_str.size) {
                            throw std::runtime_error("Failed to parse full position value to long: "
                                + std::string(pos_str.buf));
                        }
                    }
                }
            }
        }


        // Determine whether to write this entry to the index or whether it is inside the current bin
        if (line_number % index_configuration.entries_per_bin == 0) {
            struct index_entry entry;
            // char *endptr;
            entry.reference_name_idx = ref_name_map.reference_to_int(reference_name.buf);
            entry.position = pos;
            entry.byte_offset = line_byte_offset;
            debugf("Writing entry to index %d %d %ld\n",
                entry.reference_name_idx, entry.position, entry.byte_offset);

            long index_entry_bytes = write_index_entry(output_file, &entry);
            total_write_bytes += index_entry_bytes;
            debugf("Wrote %ld bytes to index, new size: %ld\n", index_entry_bytes, total_write_bytes);
        } else {
            debugf("Not writing entry to index\n");
        }
        line_number++;
    }

    fclose(output_file);
    fclose(input_file);
}


void create_binning_index_fd(const std::string& compressed_input_filename, const std::string& index_filename, VcfPackedBinningIndexConfiguration& index_configuration) {
    int input_fd = open(compressed_input_filename.c_str(), O_RDONLY);
    if (input_fd < 0) {
        perror("open");
        throw std::runtime_error("Failed to open file: " + compressed_input_filename);
    }
    int output_fd = open(index_filename.c_str(), DEFAULT_FILE_CREATE_FLAGS, DEFAULT_FILE_CREATE_MODE);
    if (output_fd < 0) {
        perror("open");
        throw std::runtime_error("Failed to open output file: " + index_filename);
    }

    int status;
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);

    reference_name_map ref_name_map;

    std::vector<byte_t> line_bytes;
    line_bytes.reserve(4096);

    // Counter to handle entries per bin
    size_t line_number = 0;

    // Iterate through lines of compressed file
    while (true) {
        long line_byte_offset = tellfd(input_fd);

        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                line_byte_offset,
                line_byte_offset);

        // right now all length headers are 4 bytes
        // TODO update to interpret variable-length length headers
        uint8_t line_length_header_bytes[4] = {0,0,0,0};

        status = read(input_fd, &line_length_header_bytes, 4);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
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
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                required_columns_length_header_bytes[0],
                required_columns_length_header_bytes[1],
                required_columns_length_header_bytes[2],
                required_columns_length_header_bytes[3]);

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                tellfd(input_fd),
                tellfd(input_fd));

        LineLengthHeader line_length_header;
        line_length_header.set_extension_count(3); // TODO interpret
        line_length_header.deserialize(line_length_header_bytes);

        debugf("Line length: %d\n", line_length_header.length);

        // collect bytes to end of line
        size_t i = 0;
        line_bytes.clear();
        if (line_bytes.capacity() < line_length_header.length + read_bytes) {
            line_bytes.reserve(line_length_header.length + read_bytes);
        }

        // two uint64 placeholders for diffs to previous, next line
        for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
            line_bytes.push_back(0);
        }

        line_bytes.push_back(line_length_header_bytes[0]);
        line_bytes.push_back(line_length_header_bytes[1]);
        line_bytes.push_back(line_length_header_bytes[2]);
        line_bytes.push_back(line_length_header_bytes[3]);
        line_bytes.push_back(required_columns_length_header_bytes[0]);
        line_bytes.push_back(required_columns_length_header_bytes[1]);
        line_bytes.push_back(required_columns_length_header_bytes[2]);
        line_bytes.push_back(required_columns_length_header_bytes[3]);

        debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        bool got_reference_name = false;
        string_t reference_name;
        string_reserve(&reference_name, 32);

        bool got_pos = false;
        string_t pos_str;
        string_reserve(&pos_str, 32);

        size_t pos = 0;

        // iterate over the rest of the line
        // subtract 4 due to length of required_columns_length_header_bytes, which are already read
        while (i++ < line_length_header.length - 4) {
            unsigned char b;
            if (read(input_fd, &b, 1) <= 0) {
                std::string msg = string_format(
                    "Unexpectedly reached end of compressed file, line header said %d, but only read %d bytes from line",
                    line_length_header.length, i);
                throw VcfValidationError(msg.c_str());
            }
            line_bytes.push_back(b);

            if (!got_reference_name) {
                if (b != '\t') {
                    string_appendc(&reference_name, b);
                } else {
                    if (reference_name.size == 0) {
                        throw std::runtime_error("Line did not contain a reference name");
                    } else {
                        debugf("Got reference name: %s\n", reference_name.buf);
                        got_reference_name = true;
                    }
                }
            } else if (!got_pos) {
                if (b != '\t') {
                    // pos_str += b;
                    string_appendc(&pos_str, b);
                } else {
                    // if (pos_str.size() == 0) {
                    if (pos_str.size == 0) {
                        throw std::runtime_error("Line did not contain a position value");
                    } else {
                        debugf("Got position: %s\n", pos_str.buf);
                        got_pos = true;
                        char *endptr = NULL;
                        pos = strtoul(pos_str.buf, &endptr, 10);
                        if (endptr != pos_str.buf + pos_str.size) {
                            throw std::runtime_error("Failed to parse full position value to long: "
                                + std::string(pos_str.buf));
                        }
                    }
                }
            }
        }

        // Determine whether to write this entry to the index or whether it is inside the current bin
        if (line_number % index_configuration.entries_per_bin == 0) {
            struct index_entry entry;
            // char *endptr;
            entry.reference_name_idx = ref_name_map.reference_to_int(reference_name.buf);
            entry.position = pos;
            entry.byte_offset = line_byte_offset;
            debugf("Writing entry to index %d %d %ld\n",
                entry.reference_name_idx, entry.position, entry.byte_offset);
            write_index_entry_fd(output_fd, &entry);
        } else {
            debugf("Not writing entry to index\n");
        }
        line_number++;
    }

    close(output_fd);
    close(input_fd);
}



/**
 * Reads a reference name and position column from uncompressed TSV file input.
 *
 * Returns the number of bytes read.
 *
 * Returns negative on error, if two TSV columns are not read, or the pos field is not an integer.
 */
int read_reference_name_and_pos(int input_fd, std::string *reference_name_out, uint64_t *pos_out) {
    int status, read_bytes = 0;
    char c;
    std::string ref;
    while (true) {
        status = read(input_fd, &c, 1);
        if (status < 0) {
            debugf("Failed to read from file\n");
            return -1;
        } else if (status == 0) {
            debugf("Unexpected EOF while reading for reference name!\n");
            return -1;
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
            debugf("Failed to read from file\n");
            return -1;
        } else if (status == 0) {
            debugf("Unexpected EOF while reading for position!");
            return -1;
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
        debugf("Malformed line, pos column must be int but got: %s", pos_str.c_str());
        return -1;
        // throw VcfValidationError(string_format("Malformed line, pos column must be int but got: %s", pos_str.c_str()).c_str());
    }
    *reference_name_out = ref;
    *pos_out = pos;
    return read_bytes;
}


void query_binned_index_binarysearch(const std::string& compressed_filename, VcfCoordinateQuery query) {
    int status;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    std::string index_filename = compressed_filename + VCFC_BINNING_INDEX_EXTENSION;
    debugf("Opening %s\n", compressed_filename.c_str());
    FILE *compressed_file = fopen(compressed_filename.c_str(), "r");
    if (compressed_file == NULL) {
        perror("fopen");
        debugf("Failed to open input file: %s\n", compressed_filename.c_str());
        return;
    }
    debugf("Opening %s\n", index_filename.c_str());
    FILE *index_file = fopen(index_filename.c_str(), "r");
    if (index_file == NULL) {
        perror("fopen");
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }

    debugf("Parsing metadata lines and header line\n");
    VcfCompressionSchema schema;
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(compressed_file, meta_header_lines, schema);

    struct index_entry start_entry;
    memset(&start_entry, 0, sizeof(start_entry));

    reference_name_map ref_name_map;
    uint8_t query_reference_name_idx = ref_name_map.reference_to_int(query.get_reference_name());

    // already checked file existence
    long index_size = file_size(index_filename.c_str());

    if (index_size % struct_index_entry_size != 0) {
        throw std::runtime_error(string_format(
            "Index size %ld was not a multiple of entry size: %ld",
            index_size, struct_index_entry_size));
    }

    long entry_count = index_size / struct_index_entry_size;
    debugf("Index of size %ld has %ld entries\n", index_size, entry_count);

    // Find bin before the first bin start that matches the query
    // long start_entry_address = 0;

    auto get_mid = [](long start, long end) -> int {
        debugf("get_mid %ld, %ld\n", start, end);
        if (start == end) {
            return start;
        } else if (start > end) {
            throw std::runtime_error("start was greater than end");
        } else {
            return (int)((start + end) / 2);
        }
    };

    // Entry index for binary search range
    long search_start = 0;
    long search_end = entry_count - 1;
    long search_mid = get_mid(search_start, search_end);
    struct index_entry entry;

    // int bytes_read;
    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    while (true) {
        search_mid = get_mid(search_start, search_end);
        long mid_offset = search_mid * struct_index_entry_size;
        debugf("Search mid_index = %ld, mid_offset = %ld\n", search_mid, mid_offset);

        fseek(index_file, mid_offset, SEEK_SET);

        int bytes_read = 0;
        if (read_index_entry(index_file, &entry, &bytes_read) != 0) {
            debugf("Failed to read index entry from index file");
            fclose(compressed_file);
            fclose(index_file);
            return;
        }

        if (search_start >= search_end) {
            // base case

            // if entry is greater than query, go back one entry
            if (search_mid > 0 &&
                    (entry.reference_name_idx > query_reference_name_idx
                        || (entry.reference_name_idx == query_reference_name_idx
                            && entry.position >= query.get_start_position()))) {
                search_mid = search_mid - 1;
                long mid_offset = search_mid * struct_index_entry_size;
                debugf("Search mid_index = %ld, mid_offset = %ld\n", search_mid, mid_offset);

                fseek(index_file, mid_offset, SEEK_SET);

                if (read_index_entry(index_file, &entry, &bytes_read) != 0) {
                    debugf("Failed to read index entry from index file");
                    fclose(compressed_file);
                    fclose(index_file);
                    return;
                }
            }
            debugf("Found desired bin %ld\n", search_mid);
            break;
        }

        debugf("entry reference_name_idx=%d, position=%u\n",
            entry.reference_name_idx, entry.position);

        // Read entry, determine if after or before query
        if (entry.reference_name_idx == query_reference_name_idx
                && entry.position == query.get_start_position()) {
            // This bin entry starts exactly at the query start range
            debugf("Entry is exactly at start of query range\n");
            break;
        }
        // This bin entry is after the query start range, so search before
        else if (entry.reference_name_idx > query_reference_name_idx
                || (entry.reference_name_idx == query_reference_name_idx
                        && entry.position > query.get_start_position())) {
            debugf("Entry is after the start of the query range\n");
            if (search_mid == 0) {
                debugf("Query starts before start of index, use first index entry\n");
                break;
            }
            search_end = search_mid - 1;
            debugf("Updated search_end to %ld\n", search_end);
        } else if (entry.reference_name_idx < query_reference_name_idx
                || (entry.reference_name_idx == query_reference_name_idx
                        && entry.position < query.get_start_position())) {
            debugf("Entry is before the start of the query range\n");
            search_start = search_mid + 1;
            debugf("Updated search_start to %ld\n", search_start);
        } else {
            throw std::runtime_error(
                "Unknown state, index wasn't equal, greater, or less");
        }

        // if (search_start >= search_end) {
        //     search_mid = search_end;
        //     break;
        // }
    }

    debugf("Got bin index %ld, reference = %d, pos = %d\n",
        search_mid, entry.reference_name_idx, entry.position);



    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("index_search time: %lu\n", duration.count());
    #endif

    std::string linebuf;
    linebuf.reserve(4 * 4096);

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    enum query_state {
        BEFORE_QUERY, IN_QUERY, AFTER_QUERY
    };
    enum query_state current_state = query_state::BEFORE_QUERY;
    #endif

    if (entry_count > 0) { // unnecessary check, would have failed earlier
        debugf("entry reference_name_idx = %u, position = %u, byte_offset = %lu\n",
            entry.reference_name_idx, entry.position, entry.byte_offset);

        // Iterate through data file from entry.byte_offset
        fseek(compressed_file, entry.byte_offset, SEEK_SET);
        // struct compressed_line_length_headers length_headers;

        // Record the number of lines scanned before hitting the desired range
        int before_count = 0;

        while (true) {
            linebuf.clear();
            size_t compressed_line_length;
            status = decompress2_data_line(compressed_file, schema, linebuf, &compressed_line_length);
            if (status == 0) {
                // EOF
                debugf("End of input file\n");
                break;
            } else if (status < 0) {
                throw VcfValidationError("Error, failed to read next line from compressed input file");
            }
            SplitIterator spi(linebuf, "\t");
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string reference_name = spi.next();
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string pos_str = spi.next();
            bool conversion_success = false;
            uint64_t pos = str_to_uint64(pos_str, conversion_success);
            if (!conversion_success) {
                throw VcfValidationError(("Failed to parse integer pos from " + pos_str).c_str());
            }

            std::string id = spi.next();
            std::string ref = spi.next();
            std::string alt = spi.next();
            // TODO optimize, dont need these if alt is not structural
            std::string qual = spi.next();
            std::string filter = spi.next();
            std::string info = spi.next();

            long end_position = compute_end_position(
                (long) pos,
                reference_name,
                alt,
                info);

            debugf("Checking reference_name = %s, pos = %lu against query reference_name = %s, start = %lu, end = %lu\n",
                reference_name.c_str(), pos, query.get_reference_name().c_str(), query.get_start_position(), query.get_end_position());

            int query_compare_to_line = query.compare_to_range(reference_name, pos, end_position);

            // if (query.matches(reference_name, pos)) {
            if (query_compare_to_line == 0 ) {
                // If previous state was before, now not before, output decompress_seeking timer
                #ifdef TIMING
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }
                current_state = query_state::IN_QUERY;
                #endif

                debugf("Query matched line, outputting\n");
                printf(linebuf.c_str());
            // } else if (pos > query.get_end_position()) {
            } else if (query_compare_to_line < 0) { // line is after query
                #ifdef TIMING
                // If previous state was before, now not before, output decompress_seeking timer
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }

                current_state = query_state::AFTER_QUERY;
                #endif

                debugf("Query did not match, state is now AFTER_QUERY\n");
                // debugf("Reached end of query range\n");
                break;
            // } else if (pos < query.get_start_position()) {
            } else if (query_compare_to_line > 0) { // line is before query
                debugf("Query did not match, state is still BEFORE_QUERY\n");
                before_count++;
                continue;
            }
        }
        debugf("lines decompressed before query = %d\n", before_count);

    } else {
        debugf("Index was empty\n");
    }

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING decompress_iteration: %lu\n", duration.count());
    #endif

    fclose(index_file);
    fclose(compressed_file);
}


void query_binned_index_FILE(const std::string& compressed_filename, VcfCoordinateQuery query) {
    int status;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    std::string index_filename = compressed_filename + VCFC_BINNING_INDEX_EXTENSION;
    // int compressed_fd = open(compressed_filename.c_str(), O_RDONLY);
    // if (compressed_fd < 0) {
    //     perror("open");
    //     debugf("Failed to open input file: %s\n", compressed_filename.c_str());
    //     return;
    // }
    debugf("Opening compressed file\n");
    FILE *compressed_file = fopen(compressed_filename.c_str(), "r");
    if (compressed_file == NULL) {
        perror("fopen");
        debugf("Failed to open input file: %s\n", compressed_filename.c_str());
        return;
    }
    debugf("Successfully opened file\n");
    if (!file_exists(index_filename.c_str())) {
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }
    debugf("Opening index file\n");
    int index_fd = open(index_filename.c_str(), O_RDONLY);
    if (index_fd < 0) {
        perror("open");
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING file_open: %lu\n", duration.count());
    #endif

    debugf("Parsing metadata lines and header line\n");
    VcfCompressionSchema schema;
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    // decompress2_metadata_headers_fd(compressed_fd, meta_header_lines, schema);
    decompress2_metadata_headers(compressed_file, meta_header_lines, schema);

    struct index_entry start_entry;
    memset(&start_entry, 0, sizeof(start_entry));

    reference_name_map ref_name_map;

    uint8_t query_reference_name_idx = ref_name_map.reference_to_int(query.get_reference_name());

    // Find bin before the first bin start that matches the query
    long start_entry_address = 0;
    bool at_least_one_entry = false;
    struct index_entry entry;

    int bytes_read;
    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    while (true) {
        size_t current_entry_address = tellfd(index_fd);
        if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
            debugf("Failed to read a full index_entry from index file");
            fclose(compressed_file);
            close(index_fd);
            return;
            // throw std::runtime_error("Failed to read a full index_entry from index file");
        }

        at_least_one_entry = true;
        // This bin entry starts exactly at the query start range
        if (entry.reference_name_idx == query_reference_name_idx
                && entry.position == query.get_start_position()) {
            debugf("Index bin starts exactly at start of query range\n");
            start_entry_address = current_entry_address;
            // current_state = query_state::IN_QUERY; // starting in the query range already
            break;
        }
        // This bin entry is after the query start range, so should start before this
        else if (entry.reference_name_idx >= query_reference_name_idx
                && entry.position >= query.get_start_position()) {
            debugf("Index bin is after the start of the query range\n");
            start_entry_address = current_entry_address;
            // If there is a preceeding bin, start there
            if (current_entry_address != 0) {
                start_entry_address -= struct_index_entry_size;
                lseek(index_fd, start_entry_address, SEEK_SET);
                // re-read previous entry
                debugf("Re-reading previous index entry\n");
                if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
                    debugf("Failed to read a full index_entry from index file");
                    fclose(compressed_file);
                    close(index_fd);
                    return;
                    // throw std::runtime_error("Failed to read a full index_entry from index file");
                }
            }
            break;
        }
    }
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING index_search: %lu\n", duration.count());
    #endif

    debugf("query start_entry_address = %ld\n", start_entry_address);
    #ifdef DEBUG
    int bin_idx = start_entry_address / struct_index_entry_size;
    debugf("bin_idx = %d\n", bin_idx);
    #endif
    // string_t linebuf;
    // string_reserve(&linebuf, 4 * 4096);
    std::string linebuf;
    linebuf.reserve(16 * 4096);

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    enum query_state {
        BEFORE_QUERY, IN_QUERY, AFTER_QUERY
    };
    enum query_state current_state = query_state::BEFORE_QUERY;
    #endif

    if (at_least_one_entry) {
        debugf("entry reference_name_idx = %u, position = %u, byte_offset = %lu\n",
            entry.reference_name_idx, entry.position, entry.byte_offset);

        // Iterate through data file
        // lseek(compressed_fd, entry.byte_offset, SEEK_SET);
        fseek(compressed_file, entry.byte_offset, SEEK_SET);
        // struct compressed_line_length_headers length_headers;

        // Record the number of lines scanned before hitting the desired range
        int before_count = 0;

        while (true) {
            // string_clear(&linebuf);
            linebuf.clear();
            size_t compressed_line_length;
            // status = decompress2_data_line_fd2_string(compressed_fd, schema, linebuf, &compressed_line_length);
            status = decompress2_data_line(compressed_file, schema, linebuf, &compressed_line_length);
            if (status == 0) {
                // EOF
                debugf("End of input file\n");
                break;
            } else if (status < 0) {
                throw VcfValidationError("Error, failed to read next line from compressed input file");
            }
            SplitIterator spi(linebuf, "\t");
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string reference_name = spi.next();
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string pos_str = spi.next();
            bool conversion_success = false;
            uint64_t pos = str_to_uint64(pos_str, conversion_success);
            if (!conversion_success) {
                throw VcfValidationError(("Failed to parse integer pos from " + pos_str).c_str());
            }
            debugf("Checking reference_name = %s, pos = %lu against query reference_name = %s, start = %lu, end = %lu\n",
                reference_name.c_str(), pos, query.get_reference_name().c_str(), query.get_start_position(), query.get_end_position());

            if (query.matches(reference_name, pos)) {
                // If previous state was before, now not before, output decompress_seeking timer
                #ifdef TIMING
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }
                current_state = query_state::IN_QUERY;
                #endif

                debugf("Query matched line, outputting\n");
                printf(linebuf.c_str());
            } else if (pos > query.get_end_position()) {
                #ifdef TIMING
                // If previous state was before, now not before, output decompress_seeking timer
                if (current_state == query_state::BEFORE_QUERY) {
                    end = std::chrono::steady_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    printf("TIMING decompress_seeking: %lu\n", duration.count());
                    start = std::chrono::steady_clock::now();
                }

                current_state = query_state::AFTER_QUERY;
                #endif

                debugf("Query did not match, state is now AFTER_QUERY\n");
                // debugf("Reached end of query range\n");
                break;
            } else if (pos < query.get_start_position()) {
                debugf("Query did not match, state is still BEFORE_QUERY\n");
                before_count++;
                continue;
            }
        }
        debugf("lines decompressed before query = %d\n", before_count);

    } else {
        debugf("Index was empty\n");
    }

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING decompress_iteration: %lu\n", duration.count());
    #endif

    // free(linebuf.buf);
    close(index_fd);
    fclose(compressed_file);
}

void query_binned_index(const std::string& compressed_filename, VcfCoordinateQuery query) {
    int status;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    std::string index_filename = compressed_filename + VCFC_BINNING_INDEX_EXTENSION;
    int compressed_fd = open(compressed_filename.c_str(), O_RDONLY);
    if (compressed_fd < 0) {
        perror("open");
        debugf("Failed to open input file: %s\n", compressed_filename.c_str());
        return;
    }
    if (!file_exists(index_filename.c_str())) {
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }
    int index_fd = open(index_filename.c_str(), O_RDONLY);
    if (index_fd < 0) {
        perror("open");
        // throw std::runtime_error("Failed to open index file: " + index_filename);
        debugf("Failed to open index file: %s\n", index_filename.c_str());
        return;
    }
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING file_open: %lu\n", duration.count());
    #endif

    debugf("Parsing metadata lines and header line\n");
    VcfCompressionSchema schema;
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(compressed_fd, meta_header_lines, schema);


    struct index_entry start_entry;
    memset(&start_entry, 0, sizeof(start_entry));

    reference_name_map ref_name_map;

    uint8_t query_reference_name_idx = ref_name_map.reference_to_int(query.get_reference_name());

    // Find bin before the first bin start that matches the query
    long start_entry_address = 0;
    bool at_least_one_entry = false;
    struct index_entry entry;
    // enum query_state {
    //     BEFORE_QUERY, IN_QUERY, AFTER_QUERY
    // };
    // enum query_state current_state = query_state::BEFORE_QUERY;
    int bytes_read;
    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif
    while (true) {
        size_t current_entry_address = tellfd(index_fd);
        if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
            debugf("Failed to read a full index_entry from index file");
            close(compressed_fd);
            close(index_fd);
            return;
            // throw std::runtime_error("Failed to read a full index_entry from index file");
        }

        at_least_one_entry = true;
        // This bin entry starts exactly at the query start range
        if (entry.reference_name_idx == query_reference_name_idx
                && entry.position == query.get_start_position()) {
            debugf("Index bin starts exactly at start of query range\n");
            start_entry_address = current_entry_address;
            // current_state = query_state::IN_QUERY; // starting in the query range already
            break;
        }
        // This bin entry is after the query start range, so should start before this
        else if (entry.reference_name_idx >= query_reference_name_idx
                && entry.position >= query.get_start_position()) {
            debugf("Index bin is after the start of the query range\n");
            start_entry_address = current_entry_address;
            // If there is a preceeding bin, start there
            if (current_entry_address != 0) {
                start_entry_address -= struct_index_entry_size;
                lseek(index_fd, start_entry_address, SEEK_SET);
                // re-read previous entry
                debugf("Re-reading previous index entry\n");
                if (read_index_entry_fd(index_fd, &entry, &bytes_read) != 0) {
                    debugf("Failed to read a full index_entry from index file");
                    close(compressed_fd);
                    close(index_fd);
                    return;
                    // throw std::runtime_error("Failed to read a full index_entry from index file");
                }
            }
            break;
        }
    }
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING index_search: %lu\n", duration.count());
    #endif

    debugf("query start_entry_address = %ld\n", start_entry_address);
    #ifdef DEBUG
    int bin_idx = start_entry_address / struct_index_entry_size;
    debugf("bin_idx = %d\n", bin_idx);
    #endif
    // string_t linebuf;
    // string_reserve(&linebuf, 4 * 4096);
    std::string linebuf;
    linebuf.reserve(4 * 4096);

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    if (at_least_one_entry) {
        debugf("entry reference_name_idx = %u, position = %u, byte_offset = %lu\n",
            entry.reference_name_idx, entry.position, entry.byte_offset);

        // Iterate through data file
        lseek(compressed_fd, entry.byte_offset, SEEK_SET);
        // struct compressed_line_length_headers length_headers;

        // Record the number of lines scanned before hitting the desired range
        int before_count = 0;

        while (true) {
            // string_clear(&linebuf);
            linebuf.clear();
            size_t compressed_line_length;
            // status = decompress2_data_line_fd2_string(compressed_fd, schema, linebuf, &compressed_line_length);
            status = decompress2_data_line_FILEwrapper(compressed_fd, schema, linebuf, &compressed_line_length);
            if (status == 0) {
                // EOF
                debugf("End of input file\n");
                break;
            } else if (status < 0) {
                throw VcfValidationError("Error, failed to read next line from compressed input file");
            }
            SplitIterator spi(linebuf, "\t");
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string reference_name = spi.next();
            if (!spi.has_next()) {
                throw VcfValidationError("Line did not match expected schema\n");
            }
            std::string pos_str = spi.next();
            bool conversion_success = false;
            uint64_t pos = str_to_uint64(pos_str, conversion_success);
            if (!conversion_success) {
                throw VcfValidationError(("Failed to parse integer pos from " + pos_str).c_str());
            }
            debugf("Checking reference_name = %s, pos = %lu against query reference_name = %s, start = %lu, end = %lu\n",
                reference_name.c_str(), pos, query.get_reference_name().c_str(), query.get_start_position(), query.get_end_position());
            if (query.matches(reference_name, pos)) {
                debugf("Query matched line, outputting\n");
                // current_state = query_state::IN_QUERY;
                printf(linebuf.c_str());
            // } else if (current_state == query_state::IN_QUERY) {
            } else if (pos > query.get_end_position()) {
                debugf("Query did not match, state is now AFTER_QUERY\n");
                // current_state = query_state::AFTER_QUERY;
                debugf("Reached end of query range\n");
                break;
            } else if (pos < query.get_start_position()) {
                debugf("Query did not match, state is still BEFORE_QUERY\n");
                before_count++;
                continue;
            }
        }
        debugf("lines decompressed before query = %d\n", before_count);

    } else {
        debugf("Index was empty\n");
    }

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("TIMING decompress_iteration: %lu\n", duration.count());
    #endif

    // free(linebuf.buf);
    close(index_fd);
    close(compressed_fd);
}


void query_compressed_file(const std::string& input_filename, VcfCoordinateQuery query) {
    debugf("Querying %s for %s:%lu-%lu\n", input_filename.c_str(),
        query.get_reference_name().c_str(),
        query.get_start_position(),
        query.get_end_position());
    int input_fd = open(input_filename.c_str(), O_RDONLY);
    int status;

    debugf("Parsing metadata lines and header line\n");
    VcfCompressionSchema schema;
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);

    size_t matched_line_count = 0;
    std::string variant_line;
    variant_line.reserve(1024 * 1024); // 1MiB
    // string_t variant_line;
    // string_reserve(&variant_line, 1024 * 1024); // 1 MiB

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
            variant_line.clear();
            // string_clear(&variant_line);
            size_t compressed_line_length = 0;

            int status = decompress2_data_line_FILEwrapper(input_fd, schema, variant_line, &compressed_line_length);
            if (status == 0) {
                // EOF
                throw std::runtime_error("Unexpected EOF");
                // break;
            }
            else if (status < 0) {
                throw std::runtime_error("Failed to decompressed data line\n");
                // break;
            }

            matched_line_count++;
            std::cout << variant_line; // newline is included in decompress2_data_line

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
    // string_t variant_line;
    // string_reserve(&variant_line, 1024 * 1024); // 1 MiB
    std::string variant_line;
    variant_line.reserve(1024 * 1024);

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
        // string_clear(&variant_line);
        variant_line.clear();
        size_t compressed_line_length = 0;
        int status = decompress2_data_line_FILEwrapper(input_fd, schema, variant_line, &compressed_line_length);
        if (status == 0) {
            break;
        } else if (status < 0) {
            throw std::runtime_error("Failed to decompress data line\n");
        }
        // split the line
        std::vector<std::string> line_terms = split_string(variant_line, "\t");
        start_position_fstream << line_terms[1];
        start_position_fstream << " ";
        start_position_fstream << std::to_string(variant_line.size());
        start_position_fstream << " ";
        start_position_fstream << std::to_string(compressed_line_length);
        start_position_fstream << "\n";

    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    start_position_fstream.flush();
    close(input_fd);
}

/**
 * Parses a coordinate query into a VcfCoordinateQuery object. Overwrites existing query input.
 *
 * Expects `s` to be of one of the following forms:
 *
 * - "<ref>"
 *
 * - "<ref>:<start>-<end>"
 *
 * Returns 0 on success, negative on error.
 */
int parse_coordinate_string(const std::string& s, VcfCoordinateQuery& query) {
    // int status;
    size_t colon_idx = s.find_first_of(':');
    if (colon_idx == std::string::npos) {
        // only passed a reference name
        debugf("query reference_name = %s\n", s.c_str());
        query = VcfCoordinateQuery(s);
        return 0;
    } else {
        std::string reference_name = s.substr(0, colon_idx - 0);
        size_t dash_idx = s.find_first_of('-', colon_idx + 1);
        if (dash_idx == std::string::npos) {
            printf("Query must contain a dash character: <ref>:<start>-<end>\n");
            return -1;
        }
        std::string start_str = s.substr(colon_idx + 1, dash_idx - (colon_idx + 1));
        std::string end_str = s.substr(dash_idx + 1);
        bool success = false;
        uint64_t start_pos = str_to_uint64(start_str, success);
        if (!success) {
            printf("Failed to parse int from start position: %s\n", start_str.c_str());
            return -1;
        }
        uint64_t end_pos = str_to_uint64(end_str, success);
        if (!success) {
            printf("Failed to parse int from end position: %s\n", end_str.c_str());
            return -1;
        }
        debugf("query reference_name = %s, start = %lu, end = %lu\n", reference_name.c_str(), start_pos, end_pos);
        query = VcfCoordinateQuery(reference_name, start_pos, end_pos);
    }

    return 0;
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
        if (!file_exists(input_filename.c_str())) {
            printf("Input file does not exist: %s\n", input_filename.c_str());
        }
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
            return 1;
        }
    } else if (action == "query") {
        std::string input_filename(argv[2]);
        if (!file_exists(input_filename.c_str())) {
            printf("Input file does not exist: %s\n", input_filename.c_str());
        }
        std::string query_input(argv[3]);
        VcfCoordinateQuery query;
        status = parse_coordinate_string(query_input, query);
        if (status != 0) {
            printf("Failed to parse query string: %s\n", query_input.c_str());
            return 1;
        }
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
        if (input_filename == output_filename) {
            throw std::runtime_error("input and output file are the same");
        }
        if (!file_exists(input_filename.c_str())) {
            printf("Input file does not exist: %s\n", input_filename.c_str());
        }
        sparsify_file(input_filename, output_filename);
    } else if (action == "sparse-query") {
        std::string input_filename(argv[2]);
        std::string query_input(argv[3]);
        VcfCoordinateQuery query;
        status = parse_coordinate_string(query_input, query);
        if (status != 0) {
            printf("Failed to parse query string: %s\n", query_input.c_str());
            return 1;
        }
        query_sparse_file_fd(input_filename, query);

    } else if (action == "create-binned-index") {
        if (argc != 4) {
            printf("Usage: ./main create-binned-index <bin-size> <compressed-filename>\n");
            return 1;
        }
        std::string bin_size_str(argv[2]);
        std::string input_filename(argv[3]);
        std::string index_filename = input_filename + VCFC_BINNING_INDEX_EXTENSION;
        bool success = false;
        size_t bin_size = str_to_uint64(bin_size_str, success);
        if (!success) {
            printf("bin size must be a positive integer\n");
            return 1;
        }

        VcfPackedBinningIndexConfiguration index_configuration(bin_size);
        // create_binned_index(input_filename, index_filename, index_configuration);
        create_binned_index2(input_filename, index_filename, index_configuration);
    } else if (action == "query-binned-index") {
        if (argc < 4) {
            printf("Usage: ./main query-binned-index <compressed-filename> <region>");
            return 1;
        }
        std::string input_filename(argv[2]);
        std::string index_filename = input_filename + VCFC_BINNING_INDEX_EXTENSION;
        if (!file_exists(input_filename.c_str())) {
            printf("Input file does not exist: %s\n", input_filename.c_str());
            return 1;
        }
        if (!file_exists(index_filename.c_str())) {
            printf("Index file does not exist: %s\n", index_filename.c_str());
            return 1;
        }

        std::string query_input(argv[3]);
        VcfCoordinateQuery query;
        status = parse_coordinate_string(query_input, query);
        if (status != 0) {
            printf("Failed to parse query string: %s\n", query_input.c_str());
            return 1;
        }
        debugf("query reference_name = %s, start = %lu, end = %lu\n",
            query.get_reference_name().c_str(), query.get_start_position(), query.get_end_position());
        // query_binned_index_FILE(input_filename, query);
        query_binned_index_binarysearch(input_filename, query);
    } else if (action == "create-sparse-index") {
        if (argc != 3) {
            printf("Usage: ./main create-sparse-index <compressed-filename>\n");
            return 1;
        }
        std::string input_filename(argv[2]);
        std::string index_filename = input_filename + VCFC_BINNING_INDEX_EXTENSION + "-sparse";

        // Drop multiplication factor from default of 4 to 1 since this is only storing index entries, not whole lines
        SparsificationConfiguration sparse_config;
        sparse_config.multiplication_factor = 1;
        sparse_config.block_size = SPARSE_EXTERNAL_INDEX_BLOCK_SIZE;

        create_sparse_binning_index(input_filename, index_filename, sparse_config);
    } else if (action == "query-sparse-index") {
        if (argc != 4) {
            printf("Usage: ./main query-sparse-index <compressed-filename> <region>\n");
            return 1;
        }
        std::string input_filename(argv[2]);
        std::string index_filename = input_filename + VCFC_BINNING_INDEX_EXTENSION + "-sparse";
        std::string query_input(argv[3]);
        VcfCoordinateQuery query;
        status = parse_coordinate_string(query_input, query);
        if (status != 0) {
            printf("Failed to parse query string: %s\n", query_input.c_str());
            return 1;
        }

        // Drop multiplication factor from default of 4 to 1 since this is only storing index entries, not whole lines
        SparsificationConfiguration sparse_config;
        sparse_config.multiplication_factor = 1;
        sparse_config.block_size = SPARSE_EXTERNAL_INDEX_BLOCK_SIZE;
        query_sparse_binned_index(input_filename, index_filename, query, sparse_config);
    }
    else if (action == "query-sparse-index") {

    }
    else {
        std::cout << "Unknown action name: " << action << std::endl;
    }
}