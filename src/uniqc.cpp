#include <iostream>
#include <map>
#include <string>

#include "utils.hpp"
#include "split_iterator.hpp"

int do_counts() {
    std::map<std::string,long> counts;
    std::string line;
    while (std::cin >> line) {
        if (counts.count(line) == 0) {
            counts[line] = 1;
        } else {
            counts[line] += 1;
        }
    }

    for (auto iter = counts.begin(); iter != counts.end(); iter++) {
        std::cout << iter->first << " " << iter->second << std::endl;
    }
}

int do_counts_by_line() {
    // std::vector<std::pair<std::string,long>> counts;
    std::string line;

    while (std::getline(std::cin, line)) {
        std::map<std::string,long> line_counts;

        std::istringstream line_stream(line);
        std::string term;
        while (line_stream >> term) {
            if (line_counts.count(term) == 0) {
                line_counts[term] = 1;
            } else {
                line_counts[term] += 1;
            }
        }

        for (auto iter = line_counts.begin(); iter != line_counts.end(); iter++) {
            std::cout << iter->first << " " << iter->second << std::endl;
        }

        // for (auto iter = line_counts.begin(); iter != line_counts.end(); iter++) {
        //     counts.push_back(std::pair<std::string,long>(iter->first, iter->second));
        // }
    }

    // for (auto iter = counts.begin(); iter != counts.end(); iter++) {
    //     std::cout << iter->first << " " << iter->second << std::endl;
    // }

    return 0;
}

int do_runs_by_line() {
    std::string line;

    while (std::getline(std::cin, line)) {
        std::vector<std::pair<std::string,long>> runs;
        std::map<std::string,long> line_counts;

        std::istringstream line_stream(line);
        std::string term;

        std::string run_term;
        long run_length = 0;
        bool first = true;
        while (line_stream >> term) {
            if (first) {
                run_term = term;
                run_length = 1;
            } else if (term == run_term) {
                run_length += 1;
            } else {
                runs.push_back(std::pair<std::string,long>(
                    run_term, run_length));
                run_length = 1;
                run_term = term;
            }
            first = false;
        }

        if (run_length > 0) {
            // trailing run
            runs.push_back(std::pair<std::string,long>(run_term, run_length));
        }

        for (auto iter = runs.begin(); iter != runs.end(); iter++) {
            std::cout << iter->first << " " << iter->second << std::endl;
        }
    }

    return 0;
}

int main(int argc, char **argv) {
    // read from stdin
    std::string line;
    std::map<std::string,long> counts;

    std::string cmd = argv[1];

    if (cmd == "counts") {
        return do_counts();
    } else if (cmd == "counts-by-line") {
        return do_counts_by_line();
    } else if (cmd == "runs-by-line") {
        return do_runs_by_line();
    } else {
        throw std::runtime_error("Unknown command: " + cmd);
    }




    // while (std::cin >> line) {
    //     if (counts.count(line) == 0) {
    //         counts[line] = 1;
    //     } else {
    //         counts[line] += 1;
    //     }
    // }

    // std::string run_line;
    // long run_length = 0;
    // bool first = true;
    // while (std::cin >> line) {
    //     if (first) {
    //         run_line = line;
    //         run_length = 1;
    //     } else if (line == run_line) {
    //         run_length += 1;
    //     } else {
    //         if (counts.count(run_line) == 0) {
    //             counts[line] = 0;
    //         }
    //         counts[run_line] += run_length;
    //         run_length = 1;
    //         run_line = line;
    //     }
    // }

    // if (run_length > 0) {
    //     // trailing run
    //     counts[run_line] += run_length;
    // }

    for (auto iter = counts.begin(); iter != counts.end(); iter++) {
        const std::string& k = iter->first;
        const long& v = iter->second;
        std::cout << k << " " << std::to_string(v) << std::endl;
    }
}