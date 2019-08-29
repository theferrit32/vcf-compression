#include <iostream>
#include <string>
#include <fstream>

//#define USAGE() {std::cerr << "./main <vcf-file>" << std::endl; return 1;}

int usage() {
    std::cerr << "./main <vcf-file>" << std::endl;
    return 1;

}

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


int main(int argc, char **argv) {
    if (argc < 2) {
        return usage();
    }

    std::string filename(argv[1]);
    std::ifstream input_fstream(filename);
    std::ofstream output_fstream(filename + ".vcfz");
    VcfLineStateMachine lineStateMachine;
    std::string linebuf;
    while (std::getline(input_fstream, linebuf)) {
        if (linebuf.size() == 0) {
            // empty input line, ignore
            continue;
        } else if (linebuf.substr(0, 2) == "##") {
            lineStateMachine.to_meta();
            // compress vcf header
        } else if (linebuf.substr(0, 1) == "#") {
            lineStateMachine.to_header();
            // insert header in raw format
            output_fstream << linebuf << std::endl;
        } else {
            // treat line as variant
        }
    }
}