#ifndef _UTILS_H
#define _UTILS_H
#include <vector>
#include <string>
#include <sstream>

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

std::string char_to_bin_string(const char c) {
    std::string output;
    for (size_t i = 0; i < 8; i++) {
        char mask = 0x1 << (7-i);
        char bit = (c & mask);
        if (bit) {
            output += "1";
        } else {
            output += "0";
        }
    }
    return output;
}

std::vector<std::string> split_string(const std::string& s, const std::string& delim) {
    std::vector<std::string> v;
    size_t idx = 0;
    size_t search_idx = 0;
    // loop through the delimiter instances
    while ((idx = s.find(delim, search_idx)) != std::string::npos) {
        //debugf("idx: %ld, prev_idx: %ld\n", idx, search_idx);
        std::string term = s.substr(search_idx, (idx - search_idx));
        //debugf("Found term: %s\n", term.c_str());
        if (term.size() > 0) {
            v.push_back(term);
        }
        // set next search index to be after the current delimiter
        search_idx = idx + delim.size();
    }
    if (search_idx < s.size()) {
        // leftover term after last delimiter
        std::string term = s.substr(search_idx, (s.size() - search_idx));
        //debugf("Adding trailing split term: %s\n", term.c_str());
        v.push_back(term);
    }
    //debugf("search_idx: %ld, s.size: %ld\n", search_idx, s.size());
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

void push_string_to_byte_vector(std::vector<byte_t>& v, const std::string& s) {
    for (size_t i = 0; i < s.size(); i++) {
        v.push_back((byte_t) s.at(i));
    }
}



#endif