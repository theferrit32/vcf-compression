#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <regex>
#include <stdexcept>

#include "utils.hpp"

//#define USAGE() {std::cerr << "./main <vcf-file>" << std::endl; return 1;}

int usage() {
    std::cerr << "./main [compress|decompress] <input_file> <output_file>" << std::endl;
    return 1;
}

// regex submatch flag to match everything not in the pattern
#define REGEX_SELECT_NOTMATCH -1
// VCF file format 4.2, 4.3 require 8 tab-separated columns at start of row
// followed by a variable number of columns depending on sample count
#define VCF_REQUIRED_COL_COUNT 8

////////////////////////////////////////////////////////////////
// Byte packing masks and flag values
//
// All uncompressed VCF input bytes are ASCII, all leading bits are 0
// so we can use the value of the first bit as a flag.
// If first bit is zero, we know it is compressed and a 0|0 genotype
#define SAMPLE_MASK_00              0b10000000
#define SAMPLE_MASKED_00            0b00000000
// If first bit is a 1, the first 3 bits are reserved for genotype flag
#define SAMPLE_MASK_01_10_11        0b11100000
#define SAMPLE_MASKED_01            0b10100000
#define SAMPLE_MASKED_10            0b11000000
#define SAMPLE_MASKED_11            0b10000000
// If first bit is a 1 and the first 3 bits are 111, this column is uncompressed
// and this byte is entirely a flag. Everything from the next byte to the next tab
// is the column value.
#define SAMPLE_MASK_UNCOMPRESSED    0b11100000
#define SAMPLE_MASKED_UNCOMPRESSED  0b11100000
// the value of the remaining 5 bits in the 0b111 case are unused
////////////////////////////////////////////////////////////////

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

    // store non-sample columns uncompressed
    // maximum of 9 of these, so fits in a SAMPLE_MASKED_UNCOMPRESSED (max 31)
    // uncompressed flag
    byte_t non_sample_uncompressed_flag = SAMPLE_MASKED_UNCOMPRESSED | 8;
    byte_vec.push_back(non_sample_uncompressed_flag); // must update this later

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
    if (terms.size() > VCF_REQUIRED_COL_COUNT) {
        // update uncompressed column count
        byte_vec[0] = byte_vec[0] + 1; // no overflow check, max=31

        std::string format = terms[8];
        byte_vec.push_back('\t');
        push_string_to_byte_vector(byte_vec, format);
        debugf("pushing format: %s\n", format.c_str());
        //byte_vec.push_back('\t');
    }
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
        std::cout << samples[i] << " ";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < samples.size(); i++) {
        std::string& sample_val = samples.at(i);
        uint8_t max_dedup_00 = 0x7F; // first bit 0 means this is a 0|0 term, 7 bits left for count (max = 127)
        uint8_t max_dedup_01_10_11 = 0x1F; // first 3 bits reserved, 5 bits left for count (max = 31)
        debugf("sample_val: %s\n", sample_val.c_str());
        if (sample_val == "0|0") {
            size_t count = 1;
            i++;
            for ( ;
                    count < max_dedup_00
                    && i < samples.size()
                    && samples.at(i) == "0|0";
                    i++) {
                count++;
            }
            i--;
            debugf("0|0 occurred %ld times\n", count);
            byte_t b = 0x0 | ((uint8_t) count);
            debugf("compressed to: %s\n", char_to_bin_string(b).c_str());
            byte_vec.push_back(b);
        } else if (sample_val == "0|1" || sample_val == "1|0" || sample_val == "1|1") {
            // not 0|0
            size_t count = 1;
            i++;
            for ( ;
                    count < max_dedup_01_10_11
                    && i < samples.size()
                    && samples.at(i) == sample_val;
                    i++ ) {
                count++;
            }
            i--;
            debugf("%s occurred %ld times\n", sample_val.c_str(), count);
            byte_t b;
            if (sample_val == "0|1") {
                b = SAMPLE_MASKED_01 | ((uint8_t) count);
            } else if (sample_val == "1|0") {
                b = SAMPLE_MASKED_10 | ((uint8_t) count);
            } else if (sample_val == "1|1") {
                b = SAMPLE_MASKED_11 | ((uint8_t) count);
            } else {
                throw std::runtime_error("Invalid state! Unknown sample " + sample_val);
            }
            debugf("compressed to: %s\n", char_to_bin_string(b).c_str());
            byte_vec.push_back(b);
        } else {
            // this sample's allele genotype was higher than ALT 1 (>= 2)
            // since this is fairly rare, by the VCF definition, don't bother compressing
            debugf("sample > 1 (%s), skipping compression\n", sample_val.c_str());
            // send a flag to note this column is *not* compressed
            uint8_t uc_val = SAMPLE_MASKED_UNCOMPRESSED | 1;
            debugf("pushing SAMPLE_MASKED_UNCOMPRESSED + count: %s\n", char_to_bin_string(uc_val).c_str());
            byte_vec.push_back(uc_val);
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
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                //delete local_readbuf;
                throw VcfValidationError("VCF Header did not have enough columns");
            }
            schema.sample_count = line_terms.size() - VCF_REQUIRED_COL_COUNT - 1;
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
    //int i1, i2;
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
        input_fstream.get();
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
            //debugf("%s", "Got metadata line\n");
            meta = true;
            meta_count++;
        } else {
            //debugf("%s", "Got header line\n");
            header = true;
            header_count++;
            // determine sample count
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
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


    //char vc;
    char rb;
    //std::string tab("\t");
    const char *tab = "\t";
    const size_t tab_len = 1;
    const std::string GT_00("0|0");
    const std::string GT_01("0|1");
    const std::string GT_10("1|0");
    const std::string GT_11("1|1");
    // keep track of how many columns we've seen
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;
    while (!input_fstream.eof()) {
        input_fstream.read(&rb, 1);
        debugf("Top loop, rb: 0b%s, line_sample_count: %ld, line_tab_count: %ld\n", char_to_bin_string(rb).c_str(), line_sample_count, line_tab_count);

        if (line_tab_count == 0 && rb != '\n') {
            // interpret as first skip flag
            uint8_t count = (rb & (~SAMPLE_MASK_UNCOMPRESSED)) & 0xFF;
            debugf("Skipping first %u columns\n", count);
            while (line_tab_count < count) {
                input_fstream.read(&rb, 1);
                output_fstream.write(&rb, 1);
                printf("%c (%02X)", rb, rb);
                if (rb == '\t') {
                    line_tab_count++;
                }
            }
            printf("\n");
            continue;
        }

        // check if we're within the uncompressed required columns
        else if (line_sample_count >= schema.sample_count && rb == '\n') {
            debugf("Reached end of line\n");
            // handle end of line
            // TODO check col count
            line_tab_count = 0;
            line_sample_count = 0;
            output_fstream.write("\n", 1);
        } else if (line_sample_count >= schema.sample_count && rb != '\n') {
            throw VcfValidationError("Reached expected sample count in line, but did not encounter newline character");
        }

        // else if (line_tab_count < VCF_REQUIRED_COL_COUNT + 1 /*|| (line_tab_count == VCF_REQUIRED_COL_COUNT && rb != '\n')*/) {
        //     // we're before the end of the line, still within uncompressed variant description columns
        //     debugf("Writing VCF required data\n");
        //     output_fstream.write(&rb, 1);
        //     while (true) {
        //         input_fstream.read(&rb, 1);
        //         if (rb == '\t') {
        //             debugf("Finished required VCF column\n");
        //             output_fstream.write(&rb, 1);
        //             line_tab_count++;
        //             break;
        //         }
        //         output_fstream.write(&rb, 1);
        //     }
        // }
        else if ((rb & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (rb & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            debugf("0|0 repeat count: %u\n", count);

            while (count--) {
                output_fstream.write(GT_00.c_str(), GT_00.size());
                line_tab_count++;
                line_sample_count++;
                if (line_sample_count < schema.sample_count) {
                    output_fstream.write(tab, tab_len);
                }
            }
        } else if ((rb & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = rb & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0;
            while (ucounter < uncompressed_count) {
                input_fstream.read(&rb, 1);
                output_fstream.write(&rb, 1);
                if (rb == '\t') {
                    ucounter++;
                    line_tab_count++;
                    line_sample_count++;
                }
            }
        } else {
            // either 0|1, 1|0, or 1|1
            byte_t masked = rb & SAMPLE_MASK_01_10_11;
            const std::string *sample_str = NULL;

            if (masked == SAMPLE_MASKED_01) {
                sample_str = &GT_01;
            } else if (masked == SAMPLE_MASKED_10) {
                sample_str = &GT_10;
            } else if (masked == SAMPLE_MASKED_11) {
                sample_str = &GT_11;
            } else {
                throw std::runtime_error("Error during decompression of compressed file, unrecognized bitmask");
            }

            // write the sample GT and increment counters
            uint8_t count = (rb & (~SAMPLE_MASK_01_10_11)) & 0xFF;
            debugf("Got %s, count: %u\n", sample_str->c_str(), count);
            while (count--) {
                output_fstream.write(sample_str->c_str(), sample_str->size());
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    output_fstream.write(tab, tab_len);
                }
                line_tab_count++;
            }

        }
    }

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