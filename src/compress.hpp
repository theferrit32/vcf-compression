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

/** Compression **/
int compress(
        const std::string& input_filename,
        const std::string& output_filename);
int compress_data_line(
        const std::string& line,
        const VcfCompressionSchema& schema,
        std::vector<byte_t>& byte_vec, bool add_newline);

/** Decompression **/
int decompress2_fd(
        const std::string& input_filename,
        const std::string& output_filename);
int decompress2_metadata_headers_fd(
        int input_fd,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema);
int decompress2_metadata_headers(
        FILE *input_file,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema);

int decompress2_data_line(
        FILE *input_file,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length);
int decompress2_data_line_FILEwrapper(
        int input_fd,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length);
/**
 * Helper function to read just the line length headers at start of line.
 */
int read_compressed_line_length_headers(
        FILE *input_file,
        struct compressed_line_length_headers *length_headers);
int read_compressed_line_length_headers_fd(
        int input_fd,
        struct compressed_line_length_headers *length_headers);
#endif