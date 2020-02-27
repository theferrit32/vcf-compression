#include <string>

#include "utils.hpp"

const char *tab = "\t";
const size_t tab_len = 1;
const std::string GT_00("0|0");
const std::string GT_01("0|1");
const std::string GT_10("1|0");
const std::string GT_11("1|1");
const int eof = std::char_traits<char>::eof();

reference_name_map::reference_name_map() {
    uint8_t ref_map_val = 1; // counter val for ref -> int map
    for (size_t ref_idx = 0; ref_idx < references.size(); ref_idx++) {
        n_map[references[ref_idx]] = ref_map_val++;
    }
}

uint8_t reference_name_map::reference_to_int(const std::string& reference_name) {
    return n_map[reference_name];
}


std::string char_to_bin_string(const char c_input) {
    std::string output;
    const uint8_t c = reinterpret_cast<const uint8_t&>(c_input);
    for (size_t i = 0; i < 8; i++) {
        uint8_t mask = 0x1 << (7-i);
        //debugf("mask = %02X\n", mask);
        uint8_t bit = (c & mask);
        if (bit) {
            output += "1";
        } else {
            output += "0";
        }
    }
    //debugf("%s input %02X output %s\n", __FUNCTION__, c, output.c_str());
    return output;
}

std::string string_to_bin_string(const std::string& str) {
    std::string output;
    for (size_t i = 0; i < str.size(); i++) {
        output += char_to_bin_string(str[i]);
    }
    return output;
}

std::string string_format(const char* format, ...) {
    char *buf = NULL;
    va_list va;
    va_start(va, format);

    int buflen = vsnprintf(buf, 0, format, va);
    va_end(va);
    buflen++; // for '\0'
    buf = (char*) calloc(buflen, sizeof(char));

    va_start(va, format);
    vsnprintf(buf, buflen, format, va);
    va_end(va);

    std::string buf_string(buf);
    free(buf);
    return buf_string;
}

int peek(FILE *stream) {
    int c = fgetc(stream);
    ungetc(c, stream);
    return c;
}

long tellfd(int fd) {
    return lseek(fd, 0, SEEK_CUR);
}

std::vector<std::string> split_string(const std::string& s, const std::string& delim, int max_split) {
    std::vector<std::string> v;
    size_t idx = 0;
    size_t search_idx = 0;
    int split_count = 0;
    // loop through the delimiter instances
    while ((idx = s.find(delim, search_idx)) != std::string::npos) {
        if (max_split > 0 && (split_count >= max_split)) {
            break;
        }
        //debugf("idx: %ld, prev_idx: %ld\n", idx, search_idx);
        std::string term = s.substr(search_idx, (idx - search_idx));
        //debugf("Found term: %s\n", term.c_str());
        if (term.size() > 0) {
            v.push_back(term);
            split_count++;
        }
        // set next search index to be after the current delimiter
        search_idx = idx + delim.size();
    }
    if (search_idx < s.size()) {
        if (max_split < 0 || (split_count < max_split)) {
            // leftover term after last delimiter
            std::string term = s.substr(search_idx, (s.size() - search_idx));
            //debugf("Adding trailing split term: %s\n", term.c_str());
            v.push_back(term);
        }
    }
    //debugf("search_idx: %ld, s.size: %ld\n", search_idx, s.size());
    return v;
}

std::vector<std::string> split_string(const std::string& s, const std::string& delim) {
    return split_string(s, delim, -1);
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

void push_string_to_byte_vector(std::vector<byte_t>& v, const std::string& s) {
    for (size_t i = 0; i < s.size(); i++) {
        v.push_back((byte_t) s.at(i));
    }
}

std::string byte_vector_to_string(const std::vector<byte_t>& v) {
    std::ostringstream os;
    bool first = true;
    for (size_t vi = 0; vi < v.size(); vi++) {
        if (!first) {
            os << " ";
        }
        first = false;
        std::string byte_str = string_format("%02X", v[vi]);
        os << byte_str.c_str();
    }
    return os.str();
}

uint64_t str_to_uint64(const std::string& s, bool& success) {
    //char *buf = (char*) calloc(s.size() + 1, sizeof(char));
    //memcpy(buf, s.c_str(), s.size());
    char *endptr = NULL;
    long val = std::strtoul(s.c_str(), &endptr, 10);
    if (endptr != (s.c_str() + s.size())) {
        //free(buf);
        success = false;
        return 0;
    }
    //free(buf);
    success = true;
    return val;
}

// Return value must be deleted by caller
void uint64_to_uint8_array(uint64_t val, uint8_t bytes[8]) {
    //uint8_t *bytes = new uint8_t[8];
    uint64_t FF_low = 0x00000000000000FF;
    bytes[0] = (val >> 56) & (FF_low << 0);
    bytes[1] = (val >> 48) & (FF_low << 0);
    bytes[2] = (val >> 40) & (FF_low << 0);
    bytes[3] = (val >> 32) & (FF_low << 0);
    bytes[4] = (val >> 24) & (FF_low << 0);
    bytes[5] = (val >> 16) & (FF_low << 0);
    bytes[6] = (val >> 8)  & (FF_low << 0);
    bytes[7] = (val >> 0)  & (FF_low << 0);
    //return bytes;
}

void uint8_array_to_uint64(uint8_t bytes[8], uint64_t *val) {
    uint64_t zero64 = 0x0000000000000000;
    uint64_t FF_low = 0x00000000000000FF;
    *val =
        (((bytes[0] | zero64) << 56) & (FF_low << 56)) |
        (((bytes[1] | zero64) << 48) & (FF_low << 48)) |
        (((bytes[2] | zero64) << 40) & (FF_low << 40)) |
        (((bytes[3] | zero64) << 32) & (FF_low << 32)) |
        (((bytes[4] | zero64) << 24) & (FF_low << 24)) |
        (((bytes[5] | zero64) << 16) & (FF_low << 16)) |
        (((bytes[6] | zero64) << 8)  & (FF_low << 8))  |
        (((bytes[7] | zero64) << 0)  & (FF_low << 0));
}

/**
 * This should be avoided as much as possible as it involves a seek back,
 * which depending on underlying kernel and hardware could be expensive if
 * done frequently.
 *
 * Returns the return from read(). Negative on error, zero on EOF, 1 on success.
 */
int peekfd(int fd, unsigned char *c) {
    unsigned char b = 0;
    int status = read(fd, &b, 1);
    *c = b;
    if (status == 1) {
        lseek(fd, -1, SEEK_CUR);
    } else if (status > 1) {
        return -1; // should never happen
    }
    return status;
}

bool eof_fd(int fd) {
    unsigned char c;
    int p = peekfd(fd, &c);
    if (p == 0) {
        return true;
    }
    return false;
}

// Provides istream::peek function for C FILEs
int peek_FILE(FILE *stream) {
    int c = fgetc(stream);
    ungetc(c, stream);
    return c;
}