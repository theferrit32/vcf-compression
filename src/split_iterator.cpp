#include <stdexcept>

#include "split_iterator.hpp"

SplitIterator::SplitIterator(const std::string& str, const std::string& delim):
    str(str), delim(delim) {}

SplitIterator::SplitIterator(const char *s, const char *delim) {
    this->str = std::string(s);
    this->delim = std::string(delim);
}

bool SplitIterator::has_next() {
    return cur_idx <= str.length();
}

std::string SplitIterator::next() {
    if (has_next()) {
        size_t next_idx = str.find(delim, cur_idx);
        if (next_idx != std::string::npos) {
            std::string r = str.substr(cur_idx, (next_idx - cur_idx));
            cur_idx = next_idx + delim.length();
            return r;
        } else {
            // no remaining delims, rest of str is last term
            std::string r = str.substr(cur_idx, str.length() - cur_idx);
            cur_idx = str.length() + 1;
            return r;
        }
    }
    throw SplitIteratorNoSuchElementError("No next element");
}
