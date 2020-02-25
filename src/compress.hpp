#ifndef _COMPRESS_H
#define _COMPRESS_H

#include <string>
#include <fstream>
#include <chrono>

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "utils.hpp"
#include "string_t.h"

int compress(const std::string& input_filename, const std::string& output_filename);
int compress_data_line(const std::string& line, const VcfCompressionSchema& schema, std::vector<byte_t>& byte_vec, bool add_newline);

int decompress2_fd(const std::string& input_filename, const std::string& output_filename);
int decompress2_metadata_headers_fd(
        int input_fd,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema);
int decompress2_data_line_fd2(
        int input_fd,
        const VcfCompressionSchema& schema,
        string_t *linebuf,
        size_t *compressed_line_length);
int decompress2_data_line(
        std::ifstream& input_fstream,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length);

#endif