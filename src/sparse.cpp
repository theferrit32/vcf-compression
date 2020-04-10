#include "sparse.hpp"

SparsificationConfiguration::SparsificationConfiguration() {}

size_t SparsificationConfiguration::compute_sparse_offset(
        const std::string& reference_name,
        size_t pos) {
    size_t block_offset =
        this->block_size
        * this->name_map.reference_to_int(reference_name)
        *  this->max_position;
    size_t in_block_offset =  pos * (this->multiplication_factor * this->block_size);
    size_t offset = block_offset + in_block_offset;
    return offset;
}

uint8_t SparsificationConfiguration::reference_to_int(const std::string& reference_name) {
    return name_map.reference_to_int(reference_name);
}


void sparsify_file(const std::string& compressed_input_filename, const std::string& sparse_filename) {
    debugf("Creating sparse indexed file %s from %s\n", sparse_filename.c_str(), compressed_input_filename.c_str());
    // int input_fd = open(compressed_input_filename.c_str(), O_RDONLY);
    FILE *input_file = fopen(compressed_input_filename.c_str(), "r");
    if (input_file == NULL) {
        perror("fopen");
        throw std::runtime_error("Failed to open file " + compressed_input_filename);
    }
    int output_fd = open(sparse_filename.c_str(), DEFAULT_FILE_CREATE_FLAGS, DEFAULT_FILE_CREATE_MODE);
    int status;

    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers(input_file, meta_header_lines, schema);

    for (auto iter = meta_header_lines.begin(); iter != meta_header_lines.end(); iter++) {
        write(output_fd, iter->c_str(), iter->size());
    }

    std::string variant_line;
    variant_line.reserve(1024 * 1024);
    std::vector<uint8_t> line_bytes;

    // sparsification configuration
    SparsificationConfiguration sparse_config;

    // placeholder for first skip count from data_start_offset to first line in data
    for (size_t initial_count_i = 0; initial_count_i < 8; initial_count_i++) {
        const char zero = 0;
        write(output_fd, &zero, 1);
    }

    long data_start_offset = tellfd(output_fd);
    debugf("data_start_offset = %lu\n", data_start_offset);
    bool is_first_line = true;
    long previous_offset = data_start_offset;

    while (true) {
        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                ftell(input_file),
                ftell(input_file));

        // right now all length headers are 4 bytes
        // TODO update to interpret variable-length length headers
        uint8_t line_length_header_bytes[4] = {0,0,0,0};

        // status = read(input_fd, &line_length_header_bytes, 4);
        status = fread(&line_length_header_bytes, sizeof(uint8_t), 4, input_file);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("line_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                line_length_header_bytes[0],
                line_length_header_bytes[1],
                line_length_header_bytes[2],
                line_length_header_bytes[3]);

        uint8_t required_columns_length_header_bytes[4] = {0,0,0,0};

        // status = read(input_fd, &required_columns_length_header_bytes, 4);
        status = fread(&required_columns_length_header_bytes, sizeof(uint8_t), 4, input_file);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                required_columns_length_header_bytes[0],
                required_columns_length_header_bytes[1],
                required_columns_length_header_bytes[2],
                required_columns_length_header_bytes[3]);

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                ftell(input_file),
                ftell(input_file));

        LineLengthHeader line_length_header;
        line_length_header.set_extension_count(3); // TODO interpret
        line_length_header.deserialize(line_length_header_bytes);

        debugf("Line length: %d\n", line_length_header.length);

        // collect bytes to end of line
        size_t i = 0;
        line_bytes.clear();
        if (line_bytes.capacity() < line_length_header.length + read_bytes) {
            line_bytes.reserve(line_length_header.length + read_bytes);
        }

        // two uint64 placeholders for diffs to previous, next line
        for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
            line_bytes.push_back(0);
        }

        line_bytes.push_back(line_length_header_bytes[0]);
        line_bytes.push_back(line_length_header_bytes[1]);
        line_bytes.push_back(line_length_header_bytes[2]);
        line_bytes.push_back(line_length_header_bytes[3]);

        line_bytes.push_back(required_columns_length_header_bytes[0]);
        line_bytes.push_back(required_columns_length_header_bytes[1]);
        line_bytes.push_back(required_columns_length_header_bytes[2]);
        line_bytes.push_back(required_columns_length_header_bytes[3]);

        debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        bool got_reference_name = false;
        std::string reference_name;
        reference_name.reserve(32);

        bool got_pos = false;
        std::string pos_str;
        pos_str.reserve(32);

        size_t pos = 0;

        // iterate over the rest of the line
        // subtract 4 due to length of required_columns_length_header_bytes, which are already read
        while (i++ < line_length_header.length - 4) {
            unsigned char b;
            // if (read(input_fd, &b, 1) <= 0) {
            if (fread(&b, 1, 1, input_file) < 1) {
                std::string msg = string_format(
                    "Unexpectedly reached end of compressed file, line header said %d, but only read %d bytes from line",
                    line_length_header.length, i);
                throw VcfValidationError(msg.c_str());
            }
            line_bytes.push_back(b);

            if (!got_reference_name) {
                if (b != '\t') {
                    reference_name.push_back(b);
                } else {
                    if (reference_name.size() == 0) {
                        throw std::runtime_error("Line did not contain a reference name");
                    } else {
                        debugf("Got reference name: %s\n", reference_name.c_str());
                        got_reference_name = true;
                    }
                }
            } else if (!got_pos) {
                if (b != '\t') {
                    pos_str.push_back(b);
                } else {
                    if (pos_str.size() == 0) {
                        throw std::runtime_error("Line did not contain a position value");
                    } else {
                        debugf("Got position: %s\n", pos_str.c_str());
                        got_pos = true;
                        char *endptr = NULL;
                        pos = strtoul(pos_str.c_str(), &endptr, 10);
                        if (endptr != pos_str.c_str() + pos_str.size()) {
                            throw std::runtime_error("Failed to parse full position value to long: " + pos_str);
                        }
                    }
                }
            }
        }

        // Compute sparse file offset for this line
        size_t variant_offset = sparse_config.compute_sparse_offset(reference_name, pos);
        size_t file_offset = variant_offset + data_start_offset;
        debugf("variant_offset = %lu, file_offset = %lu\n", variant_offset, file_offset);

        // Store uint64 number bytes back to previous line start
        uint64_t count_to_prev = file_offset - previous_offset;
        debugf("Updating line_bytes distance_to_previous = %lu\n", count_to_prev);
        uint8_t offset_bytes[8];
        uint64_to_uint8_array(count_to_prev, offset_bytes);
        for (int offi = 0; offi < 8; offi++) {
            debugf("offset_bytes[%d] = %u\n", offi, offset_bytes[offi]);
        }
        for (int offi = 0; offi < 8; offi++) {
            line_bytes[offi] = offset_bytes[offi];
        }
        debugf("line_bytes (size=%lu): %s\n", line_bytes.size(), byte_vector_to_string(line_bytes).c_str());


        // If this is not the first line, update previous next byte offset
        // with int32 number of bytes forward to this line start
        if (is_first_line) {
            // Write out the first count representing the number of bytes from the start of variant data
            // to the start of the first line. Reader should skip ahead this many bytes.
            lseek(output_fd, data_start_offset - 8, SEEK_SET);

            debugf("Writing first skip uint64_t length: %lu to file address: %lu\n",
                variant_offset, data_start_offset - 8);
            write(output_fd, &variant_offset, 8);

            is_first_line = false;
        } else {
            debugf("Updating previous line next byte diff\n");
            // long current_offset = tellfd(output_fd);
            long prev_distance_to_next_address = previous_offset + 8;
            lseek(output_fd, prev_distance_to_next_address, SEEK_SET);
            uint64_t prev_distance_to_next = file_offset - previous_offset;
            debugf("Updating previous distance_to_next at address %lu to %lu\n",
                    prev_distance_to_next_address,
                    prev_distance_to_next);
            write(output_fd, &prev_distance_to_next, 8);
            lseek(output_fd, prev_distance_to_next_address, SEEK_SET);
        }


        // Seek to this offset in the output file
        debugf("max_position = %d, block_size = %d, multiplication_factor = %d, ref_int = %d\n",
            sparse_config.max_position, sparse_config.block_size,
            sparse_config.multiplication_factor, sparse_config.reference_to_int(reference_name));
        debugf("Seeking to output file_offset: %lu\n", file_offset);
        lseek(output_fd, file_offset, SEEK_SET);

        previous_offset = file_offset; // update prev address pointer

        // Write the compressed line bytes to the sparse file
        for (auto iter = line_bytes.begin(); iter != line_bytes.end(); iter++) {
            debugf("%02X ", *iter);
            write(output_fd, &(*iter), 1);
        }
        debugf("\n");
    }
    fclose(input_file);
    close(output_fd);
}

void sparsify_file_fd(const std::string& compressed_input_filename, const std::string& sparse_filename) {
    debugf("Creating sparse indexed file %s from %s\n", sparse_filename.c_str(), compressed_input_filename.c_str());
    int input_fd = open(compressed_input_filename.c_str(), O_RDONLY);
    int output_fd = open(sparse_filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY, DEFAULT_FILE_CREATE_MODE);
    int status;

    VcfCompressionSchema schema;
    debugf("Parsing metadata lines and header line\n");
    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);

    for (auto iter = meta_header_lines.begin(); iter != meta_header_lines.end(); iter++) {
        write(output_fd, iter->c_str(), iter->size());
    }

    std::string variant_line;
    variant_line.reserve(1024 * 1024);
    std::vector<uint8_t> line_bytes;

    // sparsification configuration
    SparsificationConfiguration sparse_config;

    // placeholder for first skip count from data_start_offset to first line in data
    for (size_t initial_count_i = 0; initial_count_i < 8; initial_count_i++) {
        const char zero = 0;
        write(output_fd, &zero, 1);
    }

    long data_start_offset = tellfd(output_fd);
    debugf("data_start_offset = %lu\n", data_start_offset);
    bool is_first_line = true;
    long previous_offset = data_start_offset;

    while (true) {
        debugf("Start of line, stream positioned so next byte is at position %ld (0x%08lx)\n",
                tellfd(input_fd),
                tellfd(input_fd));

        // right now all length headers are 4 bytes
        // TODO update to interpret variable-length length headers
        uint8_t line_length_header_bytes[4] = {0,0,0,0};

        status = read(input_fd, &line_length_header_bytes, 4);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("line_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                line_length_header_bytes[0],
                line_length_header_bytes[1],
                line_length_header_bytes[2],
                line_length_header_bytes[3]);

        uint8_t required_columns_length_header_bytes[4] = {0,0,0,0};

        status = read(input_fd, &required_columns_length_header_bytes, 4);
        if (status < 0) {
            throw std::runtime_error("Error reading from file");
        } else if (status == 0) {
            debugf("Finished sparsifying file\n");
            break;
        } else if (status < 4) {
            throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
        }

        debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
                required_columns_length_header_bytes[0],
                required_columns_length_header_bytes[1],
                required_columns_length_header_bytes[2],
                required_columns_length_header_bytes[3]);

        uint64_t read_bytes = 4 + 4; // length headers

        debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
                tellfd(input_fd),
                tellfd(input_fd));

        LineLengthHeader line_length_header;
        line_length_header.set_extension_count(3); // TODO interpret
        line_length_header.deserialize(line_length_header_bytes);

        debugf("Line length: %d\n", line_length_header.length);

        // collect bytes to end of line
        size_t i = 0;
        line_bytes.clear();
        if (line_bytes.capacity() < line_length_header.length + read_bytes) {
            line_bytes.reserve(line_length_header.length + read_bytes);
        }

        // two uint64 placeholders for diffs to previous, next line
        for (size_t placeholder_i = 0; placeholder_i < 16; placeholder_i++) {
            line_bytes.push_back(0);
        }

        line_bytes.push_back(line_length_header_bytes[0]);
        line_bytes.push_back(line_length_header_bytes[1]);
        line_bytes.push_back(line_length_header_bytes[2]);
        line_bytes.push_back(line_length_header_bytes[3]);

        line_bytes.push_back(required_columns_length_header_bytes[0]);
        line_bytes.push_back(required_columns_length_header_bytes[1]);
        line_bytes.push_back(required_columns_length_header_bytes[2]);
        line_bytes.push_back(required_columns_length_header_bytes[3]);

        debugf("line_bytes with headers only: %s\n", byte_vector_to_string(line_bytes).c_str());

        // keep track of some column values as we go across, for offset calculation
        bool got_reference_name = false;
        // std::string reference_name;
        string_t reference_name;
        string_reserve(&reference_name, 32);

        bool got_pos = false;
        // std::string pos_str;
        string_t pos_str;
        string_reserve(&pos_str, 32);

        size_t pos = 0;

        // iterate over the rest of the line
        // subtract 4 due to length of required_columns_length_header_bytes, which are already read
        while (i++ < line_length_header.length - 4) {
            unsigned char b;
            if (read(input_fd, &b, 1) <= 0) {
                std::string msg = string_format(
                    "Unexpectedly reached end of compressed file, line header said %d, but only read %d bytes from line",
                    line_length_header.length, i);
                throw VcfValidationError(msg.c_str());
            }
            line_bytes.push_back(b);

            if (!got_reference_name) {
                if (b != '\t') {
                    // reference_name += b;
                    string_appendc(&reference_name, b);
                } else {
                    // if (reference_name.size() == 0) {
                    if (reference_name.size == 0) {
                        throw std::runtime_error("Line did not contain a reference name");
                    } else {
                        // debugf("Got reference name: %s\n", reference_name.c_str());
                        debugf("Got reference name: %s\n", reference_name.buf);
                        got_reference_name = true;
                    }
                }
            } else if (!got_pos) {
                if (b != '\t') {
                    // pos_str += b;
                    string_appendc(&pos_str, b);
                } else {
                    // if (pos_str.size() == 0) {
                    if (pos_str.size == 0) {
                        throw std::runtime_error("Line did not contain a position value");
                    } else {
                        debugf("Got position: %s\n", pos_str.buf);
                        got_pos = true;
                        char *endptr = NULL;
                        // pos = strtoul(pos_str.c_str(), &endptr, 10);
                        // if (endptr != pos_str.c_str() + pos_str.size()) {
                        //     throw std::runtime_error("Failed to parse full position value to long: " + pos_str);
                        // }
                        pos = strtoul(pos_str.buf, &endptr, 10);
                        if (endptr != pos_str.buf + pos_str.size) {
                            throw std::runtime_error("Failed to parse full position value to long: "
                                + std::string(pos_str.buf));
                        }
                    }
                }
            }
        }

        // Compute sparse file offset for this line
        size_t variant_offset = sparse_config.compute_sparse_offset(reference_name.buf, pos);
        size_t file_offset = variant_offset + data_start_offset;
        debugf("variant_offset = %lu, file_offset = %lu\n", variant_offset, file_offset);

        // Store uint64 number bytes back to previous line start
        uint64_t count_to_prev = file_offset - previous_offset;
        debugf("Updating line_bytes distance_to_previous = %lu\n", count_to_prev);
        uint8_t offset_bytes[8];
        uint64_to_uint8_array(count_to_prev, offset_bytes);
        for (int offi = 0; offi < 8; offi++) {
            debugf("offset_bytes[%d] = %u\n", offi, offset_bytes[offi]);
        }
        for (int offi = 0; offi < 8; offi++) {
            line_bytes[offi] = offset_bytes[offi];
        }
        debugf("line_bytes (size=%lu): %s\n", line_bytes.size(), byte_vector_to_string(line_bytes).c_str());


        // If this is not the first line, update previous next byte offset
        // with int32 number of bytes forward to this line start
        if (is_first_line) {
            // Write out the first count representing the number of bytes from the start of variant data
            // to the start of the first line. Reader should skip ahead this many bytes.
            lseek(output_fd, data_start_offset - 8, SEEK_SET);

            debugf("Writing first skip uint64_t length: %lu to file address: %lu\n",
                variant_offset, data_start_offset - 8);
            write(output_fd, &variant_offset, 8);

            is_first_line = false;
        } else {
            debugf("Updating previous line next byte diff\n");
            // long current_offset = tellfd(output_fd);
            long prev_distance_to_next_address = previous_offset + 8;
            lseek(output_fd, prev_distance_to_next_address, SEEK_SET);
            uint64_t prev_distance_to_next = file_offset - previous_offset;
            debugf("Updating previous distance_to_next at address %lu to %lu\n",
                    prev_distance_to_next_address,
                    prev_distance_to_next);
            write(output_fd, &prev_distance_to_next, 8);
            lseek(output_fd, prev_distance_to_next_address, SEEK_SET);
        }


        // Seek to this offset in the output file
        debugf("max_position = %d, block_size = %d, multiplication_factor = %d, ref_int = %d\n",
            sparse_config.max_position, sparse_config.block_size,
            sparse_config.multiplication_factor, sparse_config.reference_to_int(reference_name.buf));
        debugf("Seeking to output file_offset: %lu\n", file_offset);
        lseek(output_fd, file_offset, SEEK_SET);

        previous_offset = file_offset; // update prev address pointer

        // Write the compressed line bytes to the sparse file
        for (auto iter = line_bytes.begin(); iter != line_bytes.end(); iter++) {
            debugf("%02X ", *iter);
            write(output_fd, &(*iter), 1);
        }
        debugf("\n");
    }
    close(input_fd);
    close(output_fd);
}
