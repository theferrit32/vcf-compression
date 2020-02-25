#pragma once
#ifndef _SPARSE_H
#define _SPARSE_H

#include <string>

#include <fcntl.h>

#include "compress.hpp"

class SparsificationConfiguration {
public:
    SparsificationConfiguration();

    size_t compute_sparse_offset(
            const std::string& reference_name,
            size_t pos);

    uint8_t reference_to_int(std::string reference_name);

    // sparsification constant values
    const int multiplication_factor = 4;  // F: offset block multiplier, dependent on VCF file, number of samples
    const int block_size = 4096;          // B: 4k
    //int min_position = 1;                 // min vcf pos. VCFv4.3 defines this as 1
    const int max_position = 300000000;   // L: 300 million, should be size of largest reference
    const std::vector<std::string> references = {
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "X", "Y", "M"};

private:
    std::map<std::string,uint8_t> n_map;
};

void sparsify_file_fd(const std::string& compressed_input_filename, const std::string& sparse_filename);

#endif