#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <regex>
#include <stdexcept>

//#define USAGE() {std::cerr << "./main <vcf-file>" << std::endl; return 1;}

int usage() {
    std::cerr << "./main [compress|decompress] <input_file> <output_file>" << std::endl;
    return 1;
}

//#define DEBUG
#ifdef DEBUG
    #define debug(s) fputs(s.c_str(), stderr)
    #define debugf(s, ...) fprintf(stderr, s, __VA_ARGS__)
#else
    #define debug(s) {}
    #define debugf(s, ...) {}
#endif

typedef uint8_t byte_t;

// regex submatch flag to match everything not in the pattern
#define REGEX_SELECT_NOTMATCH -1
#define VCF_REQUIRED_COL_COUNT 8

// std::vector<std::string> split_string_regex(const std::string& s, const std::string& regex_str) {
//     //debugf("splitting string [%s] by delim [%s]\n", s.c_str(), regex_str.c_str());
//     std::vector<std::string> v;
//     std::regex r(regex_str);
//     std::sregex_token_iterator regex_iter(s.begin(), s.end(), r, REGEX_SELECT_NOTMATCH);
//     std::sregex_token_iterator iter_end;
//     for (auto iter = regex_iter; iter != iter_end; iter++) {
//         //std::string m = *iter;
//         //debugf("found match: %s\n", m.c_str());
//         v.push_back(*iter);
//     }
//     return v;
// }

std::vector<std::string> split_string(const std::string& s, const std::string& delim) {
    std::vector<std::string> v;
    size_t idx = 0;
    size_t search_idx = 0;
    // loop through the delimiter instances
    while ((idx = s.find(delim, search_idx)) != std::string::npos) {
        //debugf("idx: %ld, prev_idx: %ld\n", idx, search_idx);
        std::string term = s.substr(search_idx, (idx - search_idx));
        debugf("Found term: %s\n", term.c_str());
        if (term.size() > 0) {
            v.push_back(term);
        }

        // set next search index to be after the current delimiter
        search_idx = idx + delim.size();
    }
    if (search_idx < s.size()) {
        // leftover term after last delimiter
        std::string term = s.substr(search_idx, (s.size() - search_idx));
        debugf("Adding trailing split term: %s\n", term.c_str());
        v.push_back(term);
    }
    debugf("search_idx: %ld, s.size: %ld\n", search_idx, s.size());
    return v;
}

std::string vector_join(std::vector<std::string>& v, std::string delim) {
    std::ostringstream ss;
    bool first = true;
    for (auto iter = v.begin(); iter != v.end(); iter++) {
        if (!first) {
            ss << delim;
        } else {
            first = false;
        }
        ss << *iter;
    }
    return ss.str();
}

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

typedef struct {
    byte_t  *bytes;
    size_t len;
} byte_array;

class VcfValidationError : public std::runtime_error {
public:
    VcfValidationError()
        : runtime_error("VCF Validation Error"){};
    VcfValidationError(const char *message)
        : runtime_error(message){};
};

class VcfCompressionSchema {
public:
    VcfCompressionSchema(){};
    size_t alt_allele_count = 0;
    size_t sample_count = 0;
    std::map<std::string,byte_array> sequence_map;
};

void push_string_to_byte_vector(std::vector<byte_t>& v, const std::string& s) {
    for (size_t i = 0; i < s.size(); i++) {
        v.push_back((byte_t) s.at(i));
    }
}

byte_array byte_vector_to_bytearray(const std::vector<byte_t>& v) {
    byte_array ba;
    ba.bytes = new byte_t[v.size()];
    ba.len = v.size();
    for (size_t i = 0; i < v.size(); i++) {
        ba.bytes[i] = v[i];
    }
    return ba;
}

std::vector<byte_t> compress_data_line(const std::string& line, const VcfCompressionSchema& schema) {
    //byte_array ret;
    std::vector<byte_t> byte_vec;
    std::vector<std::string> terms = split_string(line, "\t");
    const size_t terms_size = terms.size();
    if (terms_size < VCF_REQUIRED_COL_COUNT) {
        throw VcfValidationError("VCF data line did not contain at least 8 terms");
    }
    // get ref to terms (no copy)
    const std::string& ref_name = terms[0];
    const std::string& position = terms[1];
    const std::string& id = terms[2];
    const std::string& ref_bases = terms[3];
    const std::string& alt_bases = terms[4];
    const std::string& quality = terms[5];
    const std::string& filter = terms[6];
    const std::string& info = terms[7];

    // store non-sample columns raw
    push_string_to_byte_vector(byte_vec, ref_name);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, position);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, id);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, ref_bases);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, alt_bases);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, quality);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, filter);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, info);

    debugf("quality: %s, filter: %s\n", quality.c_str(), filter.c_str());

    // handle sample columns
    // first is the FORMAT column
    // TODO parse this and use it to split up sample columns
    //const std::string& format = terms[8];

    std::vector<std::string> samples; // copy of sample data
    const size_t vcf_sample_start_colum = VCF_REQUIRED_COL_COUNT + 1;
    size_t samples_num = terms_size - (vcf_sample_start_colum);
    if (samples_num > 0) {
        byte_vec.push_back('\t');
    }
    debugf("terms: %ld, samples: %ld\n", terms_size, samples_num);
    //debugf("samples in row: %ld\n", samples_num);
    samples.resize(samples_num);
    for (size_t i = terms_size; i > vcf_sample_start_colum; i--) {
        size_t idx = i - vcf_sample_start_colum - 1;
        std::string& val = terms.back();
        //debugf("copying samples idx: %ld, val: %s\n", idx, val.c_str());
        samples[idx] = val;
        terms.pop_back();
    }

    for (size_t i = 0; i < samples.size(); i++) {
        std::string& sample_val = samples.at(i);
        uint8_t max_dedup_00 = 0x7F; // first bit 0 signals this is a 0|0 term (max = 127)
        uint8_t max_dedup_01_10_11 = 0x3F; // first 2 bits reserved, 6 bits left for count (max = 63)
        if (sample_val == "0|0") {
            size_t count = 1;
            for ( ;
                    count <= max_dedup_00
                    && i < samples.size()
                    && samples.at(i) == "0|0";
                    i++) {
                count++;
            }
            debugf("0|0 repeated %ld times\n", count);
            byte_t b = 0x0 | ((uint8_t) count);
            debugf("compressed to: %02X\n", b);
            byte_vec.push_back(b);
        } else if (sample_val == "0|1" || sample_val == "1|0" || sample_val == "1|1") {
            // not 0|0
            size_t count = 1;
            for ( ;
                    count < max_dedup_01_10_11
                    && i < samples.size()
                    && samples.at(i) == sample_val;
                    i++) {
                count++;
            }
            debugf("%s repeated %ld times\n", sample_val.c_str(), count);
            byte_t b = 0x0 | ((uint8_t) count);
            debugf("compressed to: %02X\n", b);
            byte_vec.push_back(b);
        } else {
            // this sample's allele genotype was higher than ALT 1 (>= 2), don't bother compressing
            debugf("sample > 1 (%s), skipping compression\n", sample_val.c_str());
            push_string_to_byte_vector(byte_vec, sample_val);
            byte_vec.push_back('\t');
        }
    }


    return byte_vec;
}

int compress(const std::string& input_filename, const std::string& output_filename) {
    std::ifstream input_fstream(input_filename);
    //size_t readbuf_size = 10 * 1024 * 1024; // 1 MiB
    //char *local_readbuf = new char[readbuf_size];
    //input_fstream.rdbuf()->pubsetbuf(local_readbuf, readbuf_size);
    std::ofstream output_fstream(output_filename);
    //VcfLineStateMachine lineStateMachine;
    std::string linebuf;
    VcfCompressionSchema schema;

    //input_fstream.sync_with_stdio(false);
    //output_fstream.sync_with_stdio(false);

    size_t variant_count = 0;

    while (std::getline(input_fstream, linebuf)) {
        if (linebuf.size() == 0) {
            // empty input line, ignore
            continue;
        } else if (linebuf.substr(0, 2) == "##") {
            //lineStateMachine.to_meta();
            // compress vcf header
            // TODO
            output_fstream << linebuf << "\n";
        } else if (linebuf.substr(0, 1) == "#") {
            //lineStateMachine.to_header();
            // get the number of samples from the header
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
            #ifdef DEBUG
            for (auto iter = line_terms.begin(); iter != line_terms.end(); iter++) {
                std::cout << *iter << " ";
            }
            std::cout << std::endl;
            #endif
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                //delete local_readbuf;
                throw VcfValidationError("VCF Header did not have enough columns");
            }
            schema.sample_count = line_terms.size() - VCF_REQUIRED_COL_COUNT;
            debugf("sample count: %ld\n", schema.sample_count);
            // insert header in raw format
            output_fstream << linebuf << "\n";
        } else {
            // treat line as variant
            variant_count++;
            //lineStateMachine.to_variant();
            std::vector<byte_t> compressed_line = compress_data_line(linebuf, schema);
            for (std::vector<byte_t>::iterator iter = compressed_line.begin(); iter != compressed_line.end(); iter++) {
                //const char c = *iter;
                output_fstream.write((const char*)(&(*iter)), 1);
            }
            // for (size_t ci = 0; ci < compressed_line.size(); ci++) {
            //     output_fstream << compressed_line[ci];
            // }
            output_fstream << "\n";
        }
    }
    debugf("variant count: %ld\n", variant_count);
    //delete local_readbuf;
    return 0;
}

int decompress(const std::string& input_filename, const std::string& output_filename) {
    std::ifstream input_fstream(input_filename);
    std::ofstream output_fstream(output_filename);
    //VcfLineStateMachine lineStateMachine;
    std::string linebuf;
    VcfCompressionSchema schema;
    bool meta = false, header = false;
    int i1, i2;
    char c1, c2;
    size_t meta_count = 0, header_count = 0;
    while (!input_fstream.eof()) {
        c1 = (char) input_fstream.peek();
        //debugf("got char: %c\n", c1);
        if (input_fstream.peek() != '#') {
            if (!meta || !header) {
                throw VcfValidationError("File was missing headers or metadata");
            }
            debugf("%s", "Finished reading metadata and headers\n");
            break;
        }
        i1 = input_fstream.get();
        if (header == true) {
            // got a metadata or header row after already receiving a header row
            throw VcfValidationError("Read a metadata or header row after already reading a header");
        }
        if (input_fstream.eof()) {
            throw VcfValidationError("Invalid format, empty header row");
        }
        // this is a metadata or header line, so read the rest of the line
        std::getline(input_fstream, linebuf);
        //i2 = input_fstream.peek();
        if (linebuf.size() == 0) {
            throw VcfValidationError("Invalid format, header line was has no columns");
        }
        c2 = linebuf.at(0);
        //c2 = (char) i2;
        if (c2 == '#') {
            // is a metadata line
            if (header == true) {
                throw VcfValidationError("Invalid format, metadata line found after header line");
            }
            debugf("%s", "Got metadata line\n");
            meta = true;
            meta_count++;
        } else {
            debugf("%s", "Got header line\n");
            header = true;
            header_count++;
            // determine sample count
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
            #ifdef DEBUG
            for (auto iter = line_terms.begin(); iter != line_terms.end(); iter++) {
                std::cout << *iter << " ";
            }
            std::cout << std::endl;
            #endif
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                std::string msg = "Header line had fewer columns than the minimum " + std::to_string(VCF_REQUIRED_COL_COUNT);
                throw VcfValidationError(msg.c_str());
            }
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT + 1) {
                // FORMAT column was there, but no samples
                throw VcfValidationError("Compressed VCF file had no sample data");
            }
            schema.sample_count = line_terms.size() - (VCF_REQUIRED_COL_COUNT + 1);
        }
        // c1 has been taken off stream
        // c2 is still in the stream
        // emit c1, then rest of line, then line terminator
        output_fstream.write(&c1, 1);

        output_fstream.write(linebuf.c_str(), linebuf.size());
        output_fstream.write("\n", 1);
    }
    debugf("Line counts: metadata = %ld, header = %ld\n", meta_count, header_count);
    debugf("Sample count: %ld\n", schema.sample_count);
    // do variant decompression
    //struct read_buffer {
    //    char *buf;
    //    size_t len;
    //};
    //read_buffer rbuf;
    char vc;
    // keep track of how many tab characters we've seen.
    size_t line_tab_count = 0;
    while (!input_fstream.eof()) {
        input_fstream.read(&vc, 1);


    }


/// OLD

    /*
    while (std::getline(input_fstream, linebuf)) {
        if (linebuf.size() == 0) {
            continue;
        }
        if (linebuf.substr(0, 2) == "##") {
            meta = true;
            output_fstream << linebuf << "\n";
            continue;
        } else if (meta == false) {
            throw VcfValidationError("Compressed VCF file had no metadata lines");
        }
        if (linebuf.substr(0, 1) == "#") {
            header = true;
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                throw VcfValidationError("VCF Header did not have enough columns");
            }
            schema.sample_count = line_terms.size() - VCF_REQUIRED_COL_COUNT;
            debugf("sample count: %ld\n", schema.sample_count);
            // emit header line, was not compressed
            output_fstream << linebuf << "\n";
            continue;
        } else if (header == false) {
            throw VcfValidationError("Compressed VCF file had no header line");
        }
        // Variant line

    }
    */

    return -1;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        return usage();
    }

    std::string action(argv[1]);
    std::string input_filename(argv[2]);
    std::string output_filename(argv[3]);
    if (input_filename == output_filename) {
        throw std::runtime_error("input and output file are the same");
    }
    int status;
    if (action == "compress") {
        status = compress(input_filename, output_filename);
        if (status < 0) {
            std::cerr << "Error in compression of file" << std::endl;
        }
    } else if (action == "decompress") {
        status = decompress(input_filename, output_filename);
        if (status < 0) {
            std::cerr << "Error in decompression of file" << std::endl;
        }
    } else {
        std::cout << "Unknown action name: " << action << std::endl;
    }
}