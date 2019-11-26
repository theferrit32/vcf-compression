#ifndef _UTILS_H
#define _UTILS_H
#include <vector>
#include <string>
#include <sstream>
#include <stdarg.h>
#include <string.h>

typedef uint8_t byte_t;

//#define DEBUG
#ifdef DEBUG
    #define debug(s) fputs(s.c_str(), stderr)
    #define debugf(...) fprintf(stderr, __VA_ARGS__)
#else
    #define debug(s) {}
    #define debugf(...) {}
#endif

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

std::string string_format(const char* format, ...) {
    char *buf = NULL;
    va_list va;
    va_start(va, format);

    int buflen = vsnprintf(buf, 0, format, va);
    buf = (char*) calloc(buflen, sizeof(char));

    va_start(va, format);
    vsnprintf(buf, buflen, format, va);

    std::string buf_string(buf);
    free(buf);
    return buf_string;
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

uint64_t str_to_uint64(const std::string& s, bool& success) {
    char *buf = (char*) calloc(s.size() + 1, sizeof(char));
    memcpy(buf, s.c_str(), s.size());
    char *endptr = NULL;
    long val = std::strtol(buf, &endptr, 10);
    if (endptr != (buf + s.size())) {
        free(buf);
        success = false;
        return 0;
    }
    free(buf);
    success = true;
    return val;
}

#endif