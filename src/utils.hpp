#pragma once
#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <sstream>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

typedef uint8_t byte_t;

//#define DEBUG
#ifdef DEBUG
    #define debug(s) fputs(s.c_str(), stderr)
    #define debugf(...) fprintf(stdout, __VA_ARGS__)
#else
    #define debug(s) {}
    #define debugf(...) {}
#endif

#define DEFAULT_FILE_CREATE_FLAGS (O_CREAT | O_TRUNC | O_RDWR)
#define DEFAULT_FILE_CREATE_MODE (S_IRUSR | S_IWUSR)
#define VCFC_BINNING_INDEX_EXTENSION ".vcfci"


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

extern const char *tab;
extern const size_t tab_len;
extern const std::string GT_00;
extern const std::string GT_01;
extern const std::string GT_10;
extern const std::string GT_11;
extern const int eof;

// const char *tab = "\t";
// const int tab_len = 1;
// const std::string GT_00("0|0");
// const std::string GT_01("0|1");
// const std::string GT_10("1|0");
// const std::string GT_11("1|1");
// const int eof = std::char_traits<char>::eof();


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

class reference_name_map {
public:
    reference_name_map();

    uint8_t reference_to_int(const std::string& reference_name);

private:
    const std::vector<std::string> references = {
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "X", "Y", "M"};

    std::map<std::string,uint8_t> n_map;
};

typedef struct {
    byte_t  *bytes;
    size_t len;
} byte_array;

std::string char_to_bin_string(const char c_input);
std::string string_to_bin_string(const std::string& str);
std::string string_format(const char* format, ...);
// Provides istream::peek function for C FILEs
int peek(FILE *stream);
long tellfd(int fd);

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


struct compressed_line_length_headers {
    uint32_t line_length;
    uint32_t required_columns_length;
};
extern size_t compressed_line_length_headers_size;

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
        out[0] = ((this->extension_count << 6) & 0xC0) | this->length_bytes[0];
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
            debugf("Error in deserialize, extension count was %d\n", this->extension_count);
            throw std::runtime_error(string_format(
                "Extension count %d not implemented, must be 3",
                this->extension_count));
        }

        this->length_bytes[0] = in[0] & 0x3F; // 00111111
        this->length_bytes[1] = in[1];
        this->length_bytes[2] = in[2];
        this->length_bytes[3] = in[3];

        // this->length = (this->length_bytes[0] << 24) & (0xFF << 24);
        // debugf("length = 0x%08X ", length);
        // this->length |= (this->length_bytes[1] << 16) & (0xFF << 16);
        // debugf("length = 0x%08X ", length);
        // this->length |= (this->length_bytes[2] << 8) & (0xFF << 8);
        // debugf("length = 0x%08X ", length);
        // this->length |= (this->length_bytes[3] << 0) & (0xFF << 0);
        // debugf("length = 0x%08X ", length);
        // debugf("\n");

        this->length =
            ((this->length_bytes[0] << 24) /*& (0xFF << 24)*/)
            | ((this->length_bytes[1] << 16) /*& (0xFF << 16)*/)
            | ((this->length_bytes[2] << 8) /*& (0xFF << 8)*/)
            | ((this->length_bytes[3] << 0) /*& (0xFF << 0)*/);
        debugf("length = %u, 0x%08X\n", length, length);

        // debugf("%s %02X %02X %02X %02X extension_count = %u, length = %u, bytes %02X %02X %02X %02X, bin = %s\n",
        //     __FUNCTION__, in[0], in[1], in[2], in[3],
        //     this->extension_count, this->length,
        //     this->length_bytes[0],
        //     this->length_bytes[1],
        //     this->length_bytes[2],
        //     this->length_bytes[3],
        //     string_to_bin_string(std::to_string(this->length)).c_str()
        // );
    }

    uint8_t extension_count;
    // Important, this assumes little-endian integer layout
    //union {
        uint32_t length;
        uint8_t length_bytes[LINE_LENGTH_HEADER_MAX_EXTENSION+1]; // 4
    //};
};


std::vector<std::string> split_string(const std::string& s, const std::string& delim, int max_split);
std::vector<std::string> split_string(const std::string& s, const std::string& delim);
std::string vector_join(std::vector<std::string>& v, std::string delim);

void push_string_to_byte_vector(std::vector<byte_t>& v, const std::string& s);
std::string byte_vector_to_string(const std::vector<byte_t>& v);
uint64_t str_to_uint64(const std::string& s, bool& success);

void uint64_to_uint8_array(uint64_t val, uint8_t bytes[8]);
void uint8_array_to_uint64(uint8_t bytes[8], uint64_t *val);
void uint32_to_uint8_array(uint32_t val, uint8_t bytes[4]);

/**
 * This should be avoided as much as possible as it involves a seek back,
 * which depending on underlying kernel and hardware could be expensive if
 * done frequently.
 *
 * Returns the return from read(). Negative on error, zero on EOF, 1 on success.
 */
int peekfd(int fd, unsigned char *c);

bool eof_fd(int fd);
// Provides istream::peek function for C FILEs
int peek_FILE(FILE *stream);
bool file_exists(const char *filename);
long file_size(const char *filename);

#endif