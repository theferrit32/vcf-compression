#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
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

class VcfCoordinateQuery {
public:
    VcfCoordinateQuery(
            const std::string& reference_name,
            int64_t start_position,
            int64_t end_position) {
        this->reference_name = reference_name;
        this->start_position = start_position;
        this->end_position = end_position;
    }

    VcfCoordinateQuery() {
        this->reference_name = "";
        this->start_position = -1;
        this->end_position = -1;
    }

    bool matches(const std::string& reference_name, int64_t position) {
        if (this->reference_name.size() > 0 && this->reference_name != reference_name) {
            return false;
        }
        if (this->start_position >= 0 && this->start_position > position) {
            return false;
        }
        if (this->end_position >= 0 && this->end_position < position) {
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

    int64_t get_start_position() {
        return this->start_position;
    }

    int64_t get_end_position() {
        return this->end_position;
    }

private:
    std::string reference_name;
    int64_t start_position;
    int64_t end_position;
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

#define LINE_LENGTH_HEADER_MAX_EXTENSION 3
class LineLengthHeader {
public:
    LineLengthHeader() {}

    void set_extension_count(uint8_t count) {
        if (count > LINE_LENGTH_HEADER_MAX_EXTENSION) {
            throw std::runtime_error(string_format(
                "Count exceeded max allowed %d: %d",
                LINE_LENGTH_HEADER_MAX_EXTENSION, count));
        }
        if (count != 3) {
            throw std::runtime_error(string_format(
                "Extension count %d not implemented, must be 3",
                count));
        }
        this->extension_count = count;
    }

    void set_length(uint32_t length) {
        #define LINE_LENGTH_HEADER_MAX_VALUE (uint32_t)0x3FFFFFFF
        if (length > LINE_LENGTH_HEADER_MAX_VALUE) {
            throw std::runtime_error(string_format(
                "Length exceeded max allowed %d: %d",
                LINE_LENGTH_HEADER_MAX_VALUE));
        }
        this->length_bytes[0] = (length >> 24) & 0xFF;
        this->length_bytes[1] = (length >> 16) & 0xFF;
        this->length_bytes[2] = (length >> 8) & 0xFF;
        this->length_bytes[3] = (length >> 0) & 0xFF;
        this->length = length;
        debugf("%s %u extension_count = %u, length = %u, bytes %02X %02X %02X %02X, bin = %s\n",
            __FUNCTION__, length,
            this->extension_count, this->length,
            this->length_bytes[0],
            this->length_bytes[1],
            this->length_bytes[2],
            this->length_bytes[3],
            string_to_bin_string(std::to_string(this->length)).c_str()
        );
    }

    void serialize(uint8_t out[4]) {
        out[0] = this->extension_count | this->length_bytes[0];
        out[1] = this->length_bytes[1];
        out[2] = this->length_bytes[2];
        out[3] = this->length_bytes[3];
        debugf("%s %02X %02X %02X %02X extension_count = %u, length = %u, bytes %02X %02X %02X %02X, bin = %s\n",
            __FUNCTION__, out[0], out[1], out[2], out[3],
            this->extension_count, this->length,
            this->length_bytes[0],
            this->length_bytes[1],
            this->length_bytes[2],
            this->length_bytes[3],
            string_to_bin_string(std::to_string(this->length)).c_str()
        );
    }

    void deserialize(uint8_t in[4]) {
        debugf("%s input bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n", __FUNCTION__, in[0], in[1], in[2], in[3]);
        this->extension_count = (in[0] >> 6) & 0x03;
        if (this->extension_count != 3) {
            throw std::runtime_error(string_format(
                "Extension count %d not implemented, must be 3",
                this->extension_count));
        }
        this->length_bytes[0] = in[0] & 0x3F; // 00111111
        this->length_bytes[1] = in[1];
        this->length_bytes[2] = in[2];
        this->length_bytes[3] = in[3];
        this->length = (this->length_bytes[0] << 24) & (0xFF << 24);
        debugf("length = 0x%08X ", length);
        this->length |= (this->length_bytes[1] << 16) & (0xFF << 16);
        debugf("length = 0x%08X ", length);
        this->length |= (this->length_bytes[2] << 8) & (0xFF << 8);
        debugf("length = 0x%08X ", length);
        this->length |= (this->length_bytes[3] << 0) & (0xFF << 0);
        debugf("length = 0x%08X ", length);
        debugf("\n");

        // this->
        //     ((this->length_bytes[1] << 16) & (0xFF << 16)) |
        //     ((this->length_bytes[2] << 8) & (0xFF << 8)) |
        //     ((this->length_bytes[3] << 24) & (0xFF));
        debugf("%s %02X %02X %02X %02X extension_count = %u, length = %u, bytes %02X %02X %02X %02X, bin = %s\n",
            __FUNCTION__, in[0], in[1], in[2], in[3],
            this->extension_count, this->length,
            this->length_bytes[0],
            this->length_bytes[1],
            this->length_bytes[2],
            this->length_bytes[3],
            string_to_bin_string(std::to_string(this->length)).c_str()
        );
    }

    uint8_t extension_count;
    uint32_t length;
    uint8_t length_bytes[LINE_LENGTH_HEADER_MAX_EXTENSION+1];
};

int compress_data_line(const std::string& line, const VcfCompressionSchema& schema, std::vector<byte_t>& byte_vec, bool add_newline) {
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
    // byte_t non_sample_uncompressed_flag = SAMPLE_MASKED_UNCOMPRESSED | 8;
    // byte_vec.push_back(non_sample_uncompressed_flag); // must update this later

    // TODO estimate the length required for this line so that we can shrink the bytes
    // required for storing length while making it less likely
    // to need to change the number of bytes later, avoiding a vector shift
    LineLengthHeader line_length_header;
    line_length_header.set_extension_count(3);
    uint8_t length_header_bytes[4] = {0xC0, 0, 0, 0};
    line_length_header.deserialize(length_header_bytes);
    byte_vec.push_back(length_header_bytes[0]);
    byte_vec.push_back(length_header_bytes[1]);
    byte_vec.push_back(length_header_bytes[2]);
    byte_vec.push_back(length_header_bytes[3]);

    // include another header for the length of the uncompressed region before the sample columns
    LineLengthHeader required_col_length_header;
    required_col_length_header.set_extension_count(3); // restrict to 3, max = (2 ^ (30) - 1)
    uint8_t required_col_length_header_bytes[4] = {0xC0, 0, 0, 0};
    line_length_header.deserialize(required_col_length_header_bytes);
    byte_vec.push_back(required_col_length_header_bytes[0]);
    byte_vec.push_back(required_col_length_header_bytes[1]);
    byte_vec.push_back(required_col_length_header_bytes[2]);
    byte_vec.push_back(required_col_length_header_bytes[3]);

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

    //debugf("quality: %s, filter: %s\n", quality.c_str(), filter.c_str());

    debugf("reference_name = %s, pos = %s\n", ref_name.c_str(), position.c_str());
    uint64_t required_length = 7 + ref_name.size() + position.size() + id.size() +
        ref_bases.size() + alt_bases.size() + quality.size() + filter.size() + info.size();

    // handle sample columns
    // first is the FORMAT column
    // TODO parse this and use it to split up sample columns
    if (terms.size() > VCF_REQUIRED_COL_COUNT) {
        // update uncompressed column count
        byte_vec[0] = byte_vec[0] + 1; // no overflow check, max=31, we're only at 9

        const std::string& format = terms[8];
        byte_vec.push_back('\t');
        push_string_to_byte_vector(byte_vec, format);
        debugf("pushing format: %s\n", format.c_str());
        //byte_vec.push_back('\t');
        required_length += format.size() + 1;
    }

    const size_t vcf_sample_start_colum = VCF_REQUIRED_COL_COUNT + 1;
    size_t samples_num = terms_size - (vcf_sample_start_colum);
    if (samples_num > 0) {
        byte_vec.push_back('\t');
        required_length += 1;
    }

    debugf("Updating required length to %lu\n", required_length);
    uint32_t required_length32 = (uint32_t) required_length;
    byte_vec[4] = ((required_length32 >> 24) & 0xFF) | 0xC0;
    byte_vec[5] = (required_length32 >> 16) & 0xFF;
    byte_vec[6] = (required_length32 >> 8) & 0xFF;
    byte_vec[7] = (required_length32 >> 0) & 0xFF;
    debugf("Required length header bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n", byte_vec[4], byte_vec[5], byte_vec[6], byte_vec[7]);

    std::vector<std::string> samples; // copy of sample data

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
    #ifdef DEBUG
    debugf("SAMPLES: ");
    for (size_t i = 0; i < samples.size(); i++) {
        debugf("%s ", samples[i].c_str());
    }
    debugf("\n");
    #endif


    for (size_t i = 0; i < samples.size(); i++) {
        const std::string& sample_val = samples.at(i);
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
            // loop goes to first element not matching conditions, so set i back to that one
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
            // loop goes to first element not matching conditions, so set i back to that one
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
            // TODO make this do a lookahead to see if there are multiple uncompressed columns
            // so the count here would not always just be 1
            debugf("sample > 1 (%s), skipping compression\n", sample_val.c_str());
            // send a flag to note this column is *not* compressed
            uint8_t uc_val = SAMPLE_MASKED_UNCOMPRESSED | 1;
            debugf("pushing SAMPLE_MASKED_UNCOMPRESSED + count: %s\n", char_to_bin_string(uc_val).c_str());
            byte_vec.push_back(uc_val);
            push_string_to_byte_vector(byte_vec, sample_val);
            if (i < samples.size() - 1) { // not the last sample
                byte_vec.push_back('\t');
            }
        }
    }

    if (add_newline) {
        byte_vec.push_back('\n');
    }

    // Update the start-of-line line-length header
    // include the required column length header by only subtracting 4 instead of 8
    uint32_t line_length = byte_vec.size() - 4;
    debugf("Updating line length to %u\n", line_length);
    byte_vec[0] = ((line_length >> 24) & 0xFF) | 0xC0;
    byte_vec[1] = (line_length >> 16) & 0xFF;
    byte_vec[2] = (line_length >> 8) & 0xFF;
    byte_vec[3] = (line_length >> 0) & 0xFF;
    debugf("Line length header bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n", byte_vec[0], byte_vec[1], byte_vec[2], byte_vec[3]);

    return 0;
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
    std::vector<byte_t> compressed_line;
    compressed_line.reserve(4096);

    while (std::getline(input_fstream, linebuf)) {
        if (linebuf.size() == 0) {
            // empty input line, ignore
            continue;
        } else if (linebuf[0] == '#' && linebuf[1] == '#' /*linebuf.substr(0, 2) == "##"*/) {
            //lineStateMachine.to_meta();
            // compress vcf header
            // TODO
            output_fstream << linebuf << "\n";
        } else if (linebuf[0] == '#' /*linebuf.substr(0, 1) == "#"*/) {
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
            compressed_line.clear();
            /*int status = */compress_data_line(linebuf, schema, compressed_line, true);
            if (compressed_line.back() != '\n') {
                throw std::runtime_error("No newline at end of compressed line!");
            }
            for (std::vector<byte_t>::iterator iter = compressed_line.begin(); iter != compressed_line.end(); iter++) {
                output_fstream.write((const char*)(&(*iter)), 1);
            }
            //output_fstream.write("\n", 1);
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
const int eof = std::char_traits<char>::eof();

/**
 * Reads from input_fstream to decompress one line from the vcfc file.
 * Schema must match the actual schema of the file.
 *
 * Appends the decompressed line to linebuf.
 */
int decompress2_data_line(
        std::ifstream& input_fstream,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length) {
    debugf("%s decompressing line, expecting %lu samples\n", __FUNCTION__, schema.sample_count);
    // decompress a single variant line
    int ib;
    unsigned char b;

    // keep track of how many columns we've seen
    size_t line_byte_count = 0;
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;

    // interpret first 1 byte as a skip flag, then up to 3 additional bytes
    ib = input_fstream.get();
    debugf("ib = 0x%08X\n", ib);
    if (ib == eof) {
        debugf("%s, no data in input_fstream\n", __FUNCTION__);
        return -1;
    }
    line_byte_count++;
    b = 0xFF & ib;

    // read the line length header
    uint8_t line_length_bytes[4] = {b, 0, 0, 0};
    LineLengthHeader line_length_header;
    line_length_header.deserialize(line_length_bytes);
    uint8_t extension_count = line_length_header.extension_count;
    debugf("Line length extension count: %u\n", extension_count);
    for (uint8_t i = 0; i < extension_count; i++) {
        int ib = input_fstream.get();
        if (ib == eof) {
            debugf("%s, no data in input_fstream\n", __FUNCTION__);
            return -1;
        }
        line_byte_count++;
        line_length_bytes[i+1] = 0xFF & ib;
    }
    line_length_header.deserialize(line_length_bytes);
    //uint32_t expected_line_length = line_length_header.length;


    // read the required column skip length
    ib = input_fstream.get();
    if (ib == eof) {
        debugf("%s, no data in input_fstream\n", __FUNCTION__);
        return -1;
    }
    line_byte_count++;
    b = 0xFF & ib;
    uint8_t required_length_bytes[4] = {b, 0, 0, 0};
    LineLengthHeader required_length_header;
    required_length_header.deserialize(required_length_bytes);
    uint8_t required_extension_count = required_length_header.extension_count;
    for (uint8_t i = 0; i < required_extension_count; i++) {
        int ib = input_fstream.get();
        if (ib == eof) {
            debugf("%s, no data in input_fstream\n", __FUNCTION__);
            return -1;
        }
        line_byte_count++;
        required_length_bytes[i+1] = 0xFF & ib;
    }
    debugf("Deserializing required length header\n");
    required_length_header.deserialize(required_length_bytes);
    uint32_t required_length = required_length_header.length;
    debugf("Skipping %u bytes for required columns section\n", required_length);

    for (size_t i = 0; i < required_length; i++) {
        input_fstream.read((char*)&b, 1);
        char cb = reinterpret_cast<char&>(b);
        linebuf.push_back(cb);
        if (cb == '\t') {
            line_tab_count++;
        }
    }
    debugf("Finished reading required columns: %s\n", linebuf.c_str());
    line_byte_count += required_length;


    // check to ensure we read in the appropriate number of uncompressed columns
    // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
    if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
        if (line_tab_count == VCF_REQUIRED_COL_COUNT && schema.sample_count == 0) {
            // do nothing
        } else {
            debugf("line_tab_count: %lu\n", line_tab_count);
            throw VcfValidationError("Did not read all uncompressed columns");
            return 0;
        }
    }

    //debugf("linebuf: %s\n", linebuf.c_str());

    debugf("Reading sample columns\n");
    // read the sample columns
    while (line_sample_count < schema.sample_count) {
        ib = input_fstream.peek();
        b = 0xFF & ib;
        //debugf("b: %02x\n", (unsigned char)b);

        if (ib == eof /*|| (char)ib == '\n' || (char)ib == '\t'*/) {
            std::ostringstream msg;
            msg << "Missing samples, expected " << schema.sample_count
                << ", received " << line_sample_count;
            throw VcfValidationError(msg.str().c_str());
        }
        input_fstream.get(); // remove peeked char from stream
        line_byte_count++;

        if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            uint8_t counter = count;
            debugf("0|0 repeat count: %u\n", count);

            while (counter--) {
                linebuf.append(GT_00.c_str());
                linebuf.append(tab);
            }

            // update line counts
            line_tab_count += count;
            line_sample_count += count;

            // remove last tab if at end of line
            if (line_sample_count >= schema.sample_count) {
                linebuf.pop_back();
                line_tab_count--;
            }
        } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0; // number of uncompressed columns
            while (ucounter < uncompressed_count) {
                input_fstream.read((char*)(&b), 1);
                line_byte_count++;
                // TODO handle newline ?
                if (b == '\n') {
                    // newline after uncompressed sample value
                    ucounter++;
                    line_sample_count++;
                    if (ucounter != uncompressed_count) {
                        throw VcfValidationError("Reached end of line before reading all decompressed columns");
                    }
                    input_fstream.putback(b);
                }
                else if (b == '\t') {
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
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    linebuf.append(tab);
                }
                line_tab_count++;
            }
        } // end flag cases
    } // end sample loop
    debugf("Finished reading samples\n");

    // make sure next byte is a newline
    ib = input_fstream.get();
    line_byte_count++;
    if (ib == '\n') {
        b = (char)ib;
        linebuf.push_back(b);
    } else {
        throw VcfValidationError("Sample line did not end in a newline\n");
    }

    *compressed_line_length = line_byte_count;

    return 0;
}

/**
 * Reads from input_fstream. Assumes stream position is in the metadata section.
 * Reads all following metadata lines and the header line. If the stream does not conform
 * to the VCF line specification, throws an error.
 *
 * Places all lines read into output_vector, and the schema into output_schema.
 */
int decompress2_metadata_headers(
        std::ifstream& input_fstream,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema) {
    // decompress all metadata and header lines
    bool got_meta = false, got_header = false;
    int i1, i2;
    char c1, c2;
    size_t meta_count = 0, header_count = 0;

    //std::string linebuf;
    //linebuf.reserve(4096);

    debugf("Parsing metadata lines and header line\n");
    while (true) {
        std::string linebuf;
        linebuf.reserve(4096);
        //linebuf.clear();

        if (input_fstream.eof() && (!got_header || !got_meta)) {
            throw VcfValidationError("File ended before a header or metadata line");
        }
        i1 = input_fstream.peek();
        //debugf("i1: %02x\n", i1);
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
        c1 = (char) i1;
        c2 = (char) i2;
        //debugf("i2: %02x\n", i2);

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
        linebuf.push_back(c1);
        linebuf.push_back(c2);

        while (true) {
            int i3 = input_fstream.get();
            char i3b = 0xFF & i3;
            //debugf("byte: %02x\n", i3b);
            if (i3 == eof) {
                debugf("No data lines were in the file\n");
            }
            if (i3b == '\n') {
                linebuf.push_back('\n');
                break;
            }
            if (i3b == '\t') {
                tab_count++;
                if (tab_count > VCF_REQUIRED_COL_COUNT) {
                    output_schema.sample_count++;
                }
            }
            linebuf.push_back(i3b);
        }

        debugf("Line: %s\n", linebuf.c_str());
        output_vector.push_back(std::move(linebuf));
    }
    debugf("Line counts: metadata = %ld, header = %ld\n", meta_count, header_count);
    debugf("Sample count: %ld\n", output_schema.sample_count);

    return 0;
}



int decompress2(const std::string& input_filename, const std::string& output_filename) {
    std::ifstream input_fstream(input_filename);
    std::ofstream output_fstream(output_filename);
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    //const char newline = '\n';
    decompress2_metadata_headers(input_fstream, meta_header_lines, schema);
    for (size_t i = 0; i < meta_header_lines.size(); i++) {
        // these lines still have the newline char included
        std::string& line = meta_header_lines.at(i);
        output_fstream.write(line.c_str(), line.size());
    }
    size_t variant_line_count = 0;

    std::string variant_line;
    variant_line.reserve(16 * 1024); // 16 KiB

    while (true) {
        if (input_fstream.eof() || input_fstream.peek() == eof) {
            // done
            debugf("Finished decompressing lines");
            break;
        }
        variant_line_count++;
        variant_line.clear();
        size_t compressed_line_length = 0;
        int status = decompress2_data_line(input_fstream, schema, variant_line, &compressed_line_length);
        if (status < 0 && !input_fstream.eof()) {
            throw VcfValidationError("Failed to decompress file");
        } else if (status < 0) {
            break;
        } else {
            output_fstream.write(variant_line.c_str(), variant_line.size());
        }

    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    output_fstream.flush();

    return 0;
}

void query_compressed_file(const std::string& input_filename, VcfCoordinateQuery query) {
    debugf("Querying %s for %s:%lu-%lu\n", input_filename.c_str(),
        query.get_reference_name().c_str(),
        query.get_start_position(),
        query.get_end_position());
    std::ifstream input_fstream(input_filename);
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(input_fstream, meta_header_lines, schema);

    size_t matched_line_count = 0;
    std::string variant_line;
    variant_line.reserve(1024 * 1024); // 1MiB

    while (true) {
        if (input_fstream.eof() || input_fstream.peek() == eof) {
            // done
            debugf("Finished querying file");
            break;
        }

        debugf("Start of line, stream positioned so next byte is: 0x%02X, at position %lld (0x%08llx)\n",
                (char) input_fstream.peek(),
                (long long)input_fstream.tellg(),
                (long long)input_fstream.tellg());

        // right now all length headers are 4 bytes
        // TODO update to interpret variable-length length headers
        uint8_t line_length_header_bytes[4] = {0,0,0,0};
        input_fstream.read(reinterpret_cast<char*>(line_length_header_bytes), 4);
        debugf("line_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                line_length_header_bytes[0],
                line_length_header_bytes[1],
                line_length_header_bytes[2],
                line_length_header_bytes[3]);

        uint8_t required_columns_length_header_bytes[4] = {0,0,0,0};
        input_fstream.read(reinterpret_cast<char*>(required_columns_length_header_bytes), 4);
        debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                required_columns_length_header_bytes[0],
                required_columns_length_header_bytes[1],
                required_columns_length_header_bytes[2],
                required_columns_length_header_bytes[3]);

        int64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is: 0x%02X, at position %lld (0x%08llx)\n",
                (char) input_fstream.peek(),
                (long long)input_fstream.tellg(),
                (long long)input_fstream.tellg());

        std::string ref;
        char c;
        while (true) {
            if (input_fstream.eof()) {
                throw VcfValidationError("Failed to read line, unexpected EOF while looking for reference name!");
            }

            input_fstream.read(&c, 1);
            read_bytes++;
            if (c == '\t') {
                break;
            }
            ref.push_back(c);
        }

        std::string pos_str;
        while (true) {
            if (input_fstream.eof() || input_fstream.peek() == eof) {
                throw VcfValidationError("Failed to read line, unexpected EOF while looking for position!");
            }
            input_fstream.read(&c, 1);
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
            int64_t seek_bytes = -1 * read_bytes;
            debugf("Line matches, so seeking %ld bytes\n", seek_bytes);
            input_fstream.seekg(seek_bytes, input_fstream.cur);

            debugf("Now positioned so next byte is: 0x%02X, at position %lld (0x%08llx)\n",
                    (char) input_fstream.peek(),
                    (long long)input_fstream.tellg(),
                    (long long)input_fstream.tellg());
            variant_line.clear();
            size_t compressed_line_length = 0;

            int status = decompress2_data_line(input_fstream, schema, variant_line, &compressed_line_length);
            if (status < 0 && !input_fstream.eof()) {
                throw VcfValidationError("Failed to decompress file");
            } else if (status < 0) {
                break;
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
            debugf("line length = %u, already read = %ld, so moving %u bytes from position 0x%08llx to next line\n",
                line_length, read_bytes - 4, skip_count, (long long) input_fstream.tellg());
            input_fstream.seekg(skip_count, input_fstream.cur);
        }


    } // end line loop
    debugf("matched_line_count: %lu\n", matched_line_count);

}

void gap_analysis(const std::string& input_filename) {
    std::ifstream input_fstream(input_filename);
    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    //const char newline = '\n';
    decompress2_metadata_headers(input_fstream, meta_header_lines, schema);
    // for (size_t i = 0; i < meta_header_lines.size(); i++) {
    //     // these lines still have the newline char included
    //     output_fstream.write(meta_header_lines[i].c_str(), meta_header_lines[i].size());
    // }

    size_t variant_line_count = 0;
    std::string variant_line;
    variant_line.reserve(1024 * 1024); // 1MiB

    std::string start_position_filename("start-positions.txt");
    std::ofstream start_position_fstream(start_position_filename);


    while (true) {
        if (input_fstream.eof() || input_fstream.peek() == eof) {
            // done
            debugf("Finished decompressing lines");
            break;
        }
        variant_line_count++;
        variant_line.clear();
        size_t compressed_line_length = 0;
        int status = decompress2_data_line(input_fstream, schema, variant_line, &compressed_line_length);
        if (status < 0 && !input_fstream.eof()) {
            throw VcfValidationError("Failed to decompress file");
        } else if (status < 0) {
            break;
        }
        // split the line
        std::vector<std::string> line_terms = split_string(variant_line, "\t");
        start_position_fstream << line_terms[1];
        start_position_fstream << " ";
        start_position_fstream << std::to_string(variant_line.size());
        start_position_fstream << " ";
        start_position_fstream << std::to_string(compressed_line_length);
        start_position_fstream << "\n";

        //output_fstream.write(variant_line.c_str(), variant_line.size());
    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    start_position_fstream.flush();
}


int main(int argc, char **argv) {
    // if (argc < 4) {
    //     return usage();
    // }

    // Not using stdio FILE functions, so disable sync
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
            status = decompress2(input_filename, output_filename);
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
    }
    // } else if (action == "decompress") {
    //     status = decompress2(input_filename, output_filename);
    //     if (status < 0) {
    //         std::cerr << "Error in decompression of file" << std::endl;
    //     }
    // } else if (action == "gap-analysis") {
    //     gap_analysis(input_filename);

    // }
    else {
        std::cout << "Unknown action name: " << action << std::endl;
    }
}