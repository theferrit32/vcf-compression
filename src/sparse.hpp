#pragma once
#ifndef _SPARSE_H
#define _SPARSE_H

#include <string>

#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>

#include <limits.h>

#include "compress.hpp"

#define VCFC_SPARSE_MULTIPLE_REF_PER_FILE false


class SparsificationConfiguration {
public:
    SparsificationConfiguration();

    size_t compute_sparse_offset(
            const std::string& reference_name,
            size_t pos);

    uint8_t reference_to_int(const std::string& reference_name);

    // sparsification constant values
    int multiplication_factor = 4;  // F: offset block multiplier, dependent on VCF file, number of samples
    int block_size = 4096;          // B: 4k
    //int min_position = 1;                 // min vcf pos. VCFv4.3 defines this as 1
    int max_position = 300000000;   // L: 300 million, should be size of largest reference

    reference_name_map name_map;

private:
    //std::map<std::string,uint8_t> n_map;
};

// void sparsify_file_fd(const std::string& compressed_input_filename, const std::string& sparse_filename);
void sparsify_file(const std::string& compressed_input_filename, const std::string& sparse_filename);

#endif
