#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <regex>


using namespace std;

#define DEBUG
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

vector<string> split_string_regex(const string& s, const string& regex_str) {
    //debugf("splitting string [%s] by delim [%s]\n", s.c_str(), regex_str.c_str());
    vector<string> v;
    regex r(regex_str);
    sregex_token_iterator regex_iter(s.begin(), s.end(), r, REGEX_SELECT_NOTMATCH);
    sregex_token_iterator iter_end;
    for (auto iter = regex_iter; iter != iter_end; iter++) {
        string m = *iter;
        //debugf("found match: %s\n", m.c_str());
        v.push_back(m);
    }
    return v;
}

vector<string> split_string(const string& s, const string& delim) {
    vector<string> v;
    size_t idx = 0;
    size_t prev_idx = 0;
    while ((idx = s.find(delim, prev_idx)) != std::string::npos) {
        string term = s.substr(prev_idx, (idx - prev_idx));

        v.push_back(term);
        prev_idx = idx;

    }
    return v;
}

string vector_join(vector<string>& v, string delim) {
    ostringstream ss;
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

class VcfValidationError : public runtime_error {
public:
    VcfValidationError(char const* const message) throw();
    virtual char const* what() const throw();
private:
    char const* const message;
};


class VcfMetadataLine {
public:
    VcfMetadataLine(string line) {
        this->line = line;
    }
private:
    string line;
};

class VcfHeaderLine {
public:
    VcfHeaderLine(string line) {
        this->line = line;
    }
private:
    string line;
};

class VcfDataLine {
public:
    VcfDataLine(string line) {
        vector<string> terms = split_string(line, "\t");
        if (terms.size() < 8) {

        }
    }

    string reference_name;
    unsigned long start_position;
    unsigned long end_position;
    string reference_bases;
    string alternate_base;
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
        if (current_state == HEADER || current_state == VARIANT) {
            throw std::runtime_error("Cannot move to line state META");
        }
        current_state = META;
    }

    void to_header() {
        // can't go backwards, and can't repeat the header line
        if (current_state == VARIANT || current_state == HEADER) {
            throw std::runtime_error("Cannot move to line state VARIANT");
        }
        current_state = HEADER;
    }

    bool to_variant() {
        if (current_state == UNINITIALIZED || current_state == META) {
            throw std::runtime_error("Cannot move to line state VARIANT");
        }
        current_state = VARIANT;
    }
private:
    State current_state;
};


class VcfCompressor {
public:
    VcfCompressor() {
        this->vcfLineStateMachine = VcfLineStateMachine();
    }

    byte_array process_line(string line) {

    }

private:
    VcfLineStateMachine vcfLineStateMachine;
};

int main(int argc, char **argv) {
    string filename = "/home/me/dev/stanford/1000genomes/"
        "ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf";

    ifstream input_stream(filename);
    string line;
    while (getline(input_stream, line)) {
        vector<string> line_terms = split_string_regex(line, "\\s+");
        string joined = vector_join(line_terms, ", ");
        cout << "line: [" << joined << "]" << endl;
    }
}