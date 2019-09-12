#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <regex>
#include <stdexcept>

//#define USAGE() {std::cerr << "./main <vcf-file>" << std::endl; return 1;}

int usage() {
    std::cerr << "./main <vcf-file>" << std::endl;
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

std::vector<std::string> split_string_regex(const std::string& s, const std::string& regex_str) {
    //debugf("splitting string [%s] by delim [%s]\n", s.c_str(), regex_str.c_str());
    std::vector<std::string> v;
    std::regex r(regex_str);
    std::sregex_token_iterator regex_iter(s.begin(), s.end(), r, REGEX_SELECT_NOTMATCH);
    std::sregex_token_iterator iter_end;
    for (auto iter = regex_iter; iter != iter_end; iter++) {
        //std::string m = *iter;
        //debugf("found match: %s\n", m.c_str());
        v.push_back(*iter);
    }
    return v;
}

std::vector<std::string> split_string(const std::string& s, const std::string& delim) {
    std::vector<std::string> v;
    size_t idx = 0;
    size_t search_idx = 0;
    while ((idx = s.find(delim, search_idx)) != std::string::npos) {
        //debugf("idx: %ld, prev_idx: %ld\n", idx, search_idx);
        std::string term = s.substr(search_idx, (idx - search_idx));
        //debugf("Found term: %s\n", term.c_str());
        v.push_back(term);
        search_idx = idx + 1;
    }
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
    size_t sample_count;
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

int main(int argc, char **argv) {
    if (argc < 2) {
        return usage();
    }

    std::string filename(argv[1]);
    std::ifstream input_fstream(filename);
    std::ofstream output_fstream(filename + ".vcfc");
    VcfLineStateMachine lineStateMachine;
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
            lineStateMachine.to_meta();
            // compress vcf header
            // TODO
            output_fstream << linebuf << "\n";
        } else if (linebuf.substr(0, 1) == "#") {
            lineStateMachine.to_header();
            // get the number of samples from the header
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                throw VcfValidationError("VCF Header did not have enough columns");
            }
            schema.sample_count = line_terms.size() - VCF_REQUIRED_COL_COUNT;
            debugf("sample count: %ld\n", schema.sample_count);
            // insert header in raw format
            output_fstream << linebuf << "\n";
        } else {
            // treat line as variant
            variant_count++;
            lineStateMachine.to_variant();
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
}