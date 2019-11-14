#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <regex>
#include <stdexcept>
#include <cstdio>
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
    byte_vec.reserve(512);
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
        byte_vec[0] = byte_vec[0] + 1; // no overflow check, max=31, we're only at 9

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
        //std::cout << samples[i] << " ";
        debugf("%s ", samples[i].c_str());
    }
    debugf("\n");
    //std::cout << std::endl;

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

                output_fstream.write((const char*)(&(*iter)), 1);
            }
            output_fstream << "\n";
        }
    }
    debugf("variant count: %ld\n", variant_count);
    //delete local_readbuf;
    return 0;
}

const char *tab = "\t";
const size_t tab_len = 1;
const std::string GT_00("0|0");
const std::string GT_01("0|1");
const std::string GT_10("1|0");
const std::string GT_11("1|1");
const char eof = std::char_traits<char>::eof();

int decompress2_data_line(
        std::ifstream& input_fstream,
        const VcfCompressionSchema& schema,
        std::string& linebuf) {
    // decompress a single variant line
    int ib;
    char b;

    // keep track of how many columns we've seen
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;
    //std::string linebuf;
    linebuf.reserve(1024);


    //line_tab_count = 0;
    //line_sample_count = 0;
    //linebuf.clear();

    // interpret first byte as a skip flag
    ib = input_fstream.get();
    if (ib == eof) {
        debugf("%s, no data in input_fstream\n", __FUNCTION__);
        return -1;
    }
    b = 0xFF & ib;
    std::string msg;
    msg = "first byte: " + char_to_bin_string(b) + ", hex: %02x\n";
    debugf(msg.c_str(), b);
    uint8_t count = (b & (~SAMPLE_MASK_UNCOMPRESSED)) & 0xFF;
    debugf("Skipping first %u columns\n", count);
    while (line_tab_count < count) {
        // TODO check to make sure it is not newline or eof
        input_fstream.read(&b, 1);
        //output_fstream.write(&b, 1);
        linebuf.push_back(b);
        debugf("%c (%02X)", b, b);
        if (b == '\t') {
            line_tab_count++;
        }
    }
    debugf("\n");

    // check to ensure we read in the appropriate number of uncompressed columns
    // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
    if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
        debugf("line_tab_count: %lu\n", line_tab_count);
        throw VcfValidationError("Did not read all uncompressed columns");
        return 0;
    }

    debugf("linebuf: %s\n", linebuf.c_str());


    debugf("Reading sample columns\n");
    // read the sample columns
    while (line_sample_count < schema.sample_count) {
        ib = input_fstream.peek();
        b = 0xFF & ib;
        debugf("b: %02x\n", (unsigned char)b);

        if (ib == eof /*|| (char)ib == '\n' || (char)ib == '\t'*/) {
            std::ostringstream msg;
            msg << "Missing samples, expected " << schema.sample_count
                << ", received " << line_sample_count;
            throw VcfValidationError(msg.str().c_str());
        }
        ib = input_fstream.get();

        if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            debugf("0|0 repeat count: %u\n", count);

            while (count--) {
                //output_fstream.write(GT_00.c_str(), GT_00.size());
                linebuf.append(GT_00.c_str());
                line_tab_count++;
                line_sample_count++;
                if (line_sample_count < schema.sample_count) {
                    //output_fstream.write(tab, tab_len);
                    linebuf.append(tab);
                }
            }
        } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0; // number of uncompressed columns
            while (ucounter < uncompressed_count) {
                input_fstream.read(&b, 1);
                //output_fstream.write(&b, 1);

                if (b == '\t') {
                    // don't push tabs, handled outside if
                    ucounter++;
                    line_tab_count++;
                    line_sample_count++;
                    if (line_sample_count < schema.sample_count) {
                        // if not the last term, include the tab
                        // otherwise, the tab is handled outside the if
                        linebuf.push_back(b);
                    }
                } else {
                    linebuf.push_back(b);
                }
            }
        } else {
            // either 0|1, 1|0, or 1|1
            byte_t masked = b & SAMPLE_MASK_01_10_11;
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
            uint8_t count = (b & (~SAMPLE_MASK_01_10_11)) & 0xFF;
            debugf("Got %s, count: %u\n", sample_str->c_str(), count);
            while (count--) {
                linebuf.append(*sample_str);
                //output_fstream.write(sample_str->c_str(), sample_str->size());
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    //output_fstream.write(tab, tab_len);
                    linebuf.append(tab);
                }
                line_tab_count++;
            }
        } // end flag cases
    } // end sample loop
    debugf("Finished reading samples\n");

    // make sure next byte is a newline
    ib = input_fstream.get();
    if (ib == '\n') {
        b = (char)ib;
        //output_fstream.write(&b, 1);
        linebuf.push_back(b);
    } else {
        throw VcfValidationError("Sample line did not end in a newline\n");
    }

    //debugf("variant_line_count: %lu\n", variant_line_count);
    //output_fstream.flush();
    return 0;
}

int decompress2_metadata_headers(
        std::ifstream& input_fstream,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema) {
    // decompress all metadata and header lines
    //std::string linebuf;
    //VcfCompressionSchema schema;
    bool got_meta = false, got_header = false;
    int i1, i2;
    char c1, c2;
    size_t meta_count = 0, header_count = 0;

    std::string linebuf;
    linebuf.reserve(4096);

    debugf("Parsing metadata lines and header line\n");
    while (true) {
        linebuf.clear();
        if (input_fstream.eof() && (!got_header || !got_meta)) {
            throw VcfValidationError("File ended before a header or metadata line");
        }
        i1 = input_fstream.peek();
        debugf("i1: %02x\n", i1);
        if (i1 != '#') {
            if (!got_meta || !got_header) {
                throw VcfValidationError("File was missing headers or metadata");
            }
            debugf("%s", "Finished reading metadata and headers\n");
            break;
        } else if (got_header == true) {
            // got a metadata or header row after already receiving a header row
            throw VcfValidationError("Read a metadata or header row after already reading a header");
        } else if (input_fstream.eof()) {
            throw VcfValidationError("Invalid format, empty header row");
        }
        i1 = input_fstream.get();
        i2 = input_fstream.get();
        debugf("i2: %02x\n", i2);

        if (i2 == '#') {
            if (got_header) {
                throw VcfValidationError("Got a metadata row after the CSV header");
            }
            debugf("Got a metadata line\n");
            got_meta = true;
            meta_count++;
        } else {
            if (!got_meta) {
                throw VcfValidationError("Got a header line but no metadata lines");
            }
            got_header = true;
            header_count++;
        }

        // this is a metadata or header line, so read the rest of the line
        size_t tab_count = 0;
        linebuf.clear();
        while (true) {
            int i3 = input_fstream.get();
            char i3b = 0xFF & i3;
            debugf("byte: %02x\n", i3b);
            if (i3 == eof) {
                debugf("No data lines were in the file\n");
            }
            if ((char)i3 == '\n') {
                linebuf.push_back('\n');
                break;
            }
            if ((char)i3 == '\t') {
                tab_count++;
                if (tab_count > VCF_REQUIRED_COL_COUNT) {
                    output_schema.sample_count++;
                }
            }
            char ci3 = 0XFF & i3;
            linebuf.push_back(ci3);
        }
        c1 = (char) i1;
        c2 = (char) i2;

        // push the current line data into the output vector
        std::string output_line;

        output_line.push_back(c1);
        output_line.push_back(c2);
        output_line.append(linebuf);
        output_vector.push_back(output_line);
        // output_fstream.write(&c1, 1);
        // output_fstream.write(&c2, 1);
        // output_fstream.write(linebuf.c_str(), linebuf.size());
        debugf("Line: %s\n", output_line.c_str());
    }
    debugf("Line counts: metadata = %ld, header = %ld\n", meta_count, header_count);
    debugf("Sample count: %ld\n", output_schema.sample_count);

    return 0;
}



int decompress2(const std::string& input_filename, const std::string& output_filename) {
    std::ifstream input_fstream(input_filename);
    std::ofstream output_fstream(output_filename);
    //std::string linebuf;
    VcfCompressionSchema schema;
    //bool got_meta = false, got_header = false;
    // int i1, i2;
    // char c1, c2;
    // size_t meta_count = 0, header_count = 0;
    // const char eof = std::char_traits<char>::eof();
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    const char newline = '\n';
    decompress2_metadata_headers(input_fstream, meta_header_lines, schema);
    for (size_t i = 0; i < meta_header_lines.size(); i++) {
        // these lines still have the newline char included
        output_fstream.write(meta_header_lines[i].c_str(), meta_header_lines[i].size());
        //output_fstream.write(&newline, 1);
    }
    //while (true) {
    //     if (input_fstream.eof() && (!got_header || !got_meta)) {
    //         throw VcfValidationError("File ended before a header or metadata line");
    //     }
    //     i1 = input_fstream.peek();
    //     debugf("i1: %02x\n", i1);
    //     if (i1 != '#') {
    //         if (!got_meta || !got_header) {
    //             throw VcfValidationError("File was missing headers or metadata");
    //         }
    //         debugf("%s", "Finished reading metadata and headers\n");
    //         break;
    //     } else if (got_header == true) {
    //         // got a metadata or header row after already receiving a header row
    //         throw VcfValidationError("Read a metadata or header row after already reading a header");
    //     } else if (input_fstream.eof()) {
    //         throw VcfValidationError("Invalid format, empty header row");
    //     }
    //     i1 = input_fstream.get();
    //     i2 = input_fstream.get();
    //     debugf("i2: %02x\n", i2);

    //     if (i2 == '#') {
    //         if (got_header) {
    //             throw VcfValidationError("Got a metadata row after the CSV header");
    //         }
    //         debugf("Got a metadata line\n");
    //         got_meta = true;
    //         meta_count++;
    //     } else {
    //         if (!got_meta) {
    //             throw VcfValidationError("Got a header line but no metadata lines");
    //         }
    //         got_header = true;
    //         header_count++;
    //     }

    //     // this is a metadata or header line, so read the rest of the line
    //     size_t tab_count = 0;
    //     linebuf.clear();
    //     while (true) {
    //         int i3 = input_fstream.get();
    //         char i3b = 0xFF & i3;
    //         debugf("byte: %02x\n", i3b);
    //         if (i3 == eof) {
    //             debugf("No data lines were in the file\n");
    //         }
    //         if ((char)i3 == '\n') {
    //             linebuf.push_back('\n');
    //             break;
    //         }
    //         if ((char)i3 == '\t') {
    //             tab_count++;
    //             if (tab_count > VCF_REQUIRED_COL_COUNT) {
    //                 schema.sample_count++;
    //             }
    //         }
    //         char ci3 = 0XFF & i3;
    //         linebuf.push_back(ci3);
    //     }
    //     c1 = (char) i1;
    //     c2 = (char) i2;
    //     output_fstream.write(&c1, 1);
    //     output_fstream.write(&c2, 1);
    //     output_fstream.write(linebuf.c_str(), linebuf.size());
    //     debugf("Line: %s\n", linebuf.c_str());
    // }
    // debugf("Line counts: metadata = %ld, header = %ld\n", meta_count, header_count);
    // debugf("Sample count: %ld\n", schema.sample_count);


    // TODO decompress variant lines
    // int ib;
    // char b;
    // const char *tab = "\t";
    // const size_t tab_len = 1;
    // const std::string GT_00("0|0");
    // const std::string GT_01("0|1");
    // const std::string GT_10("1|0");
    // const std::string GT_11("1|1");
    // // keep track of how many columns we've seen
    // size_t line_tab_count = 0;
    // size_t line_sample_count = 0;
    size_t variant_line_count = 0;
    //std::string variant_line;

    std::string variant_line;
    variant_line.reserve(1024 * 1024); // 1MiB

    while (true) {
        //std::string variant_line;
        // line_tab_count = 0;
        // line_sample_count = 0;
        // linebuf.clear();

        if (input_fstream.eof()) {
            // done
            debugf("Finished decompressing lines");
            break;
        }
        variant_line_count++;
        variant_line.clear();
        int status = decompress2_data_line(input_fstream, schema, variant_line);
        output_fstream.write(variant_line.c_str(), variant_line.size());
    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    output_fstream.flush();

    return 0;
}



int main(int argc, char **argv) {
    if (argc < 4) {
        return usage();
    }

    // Not using stdio FILE functions, so disable sync
    std::ios_base::sync_with_stdio(false);

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
        status = decompress2(input_filename, output_filename);
        if (status < 0) {
            std::cerr << "Error in decompression of file" << std::endl;
        }
    } else {
        std::cout << "Unknown action name: " << action << std::endl;
    }
}