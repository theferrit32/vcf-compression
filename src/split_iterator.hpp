#ifndef _SPLIT_ITERATOR_H
#define _SPLIT_ITERATOR_H
#include <stdexcept>

class SplitIterator {
public:
    class SplitIteratorNoSuchElementError : std::logic_error {
    public:
        explicit SplitIteratorNoSuchElementError(const std::string& message):
        std::logic_error(message.c_str())
        {}
    };

    SplitIterator(const std::string& str, const std::string& delim);

    SplitIterator(const char *s, const char *delim);

    bool has_next();

    std::string next();

private:
    std::string str;
    std::string delim;
    size_t cur_idx = 0;
};

#endif