#include <unordered_map>
#include <array>
#include "compress.hpp"

int compress_data_line(const std::string& line, const VcfCompressionSchema& schema, std::vector<byte_t>& byte_vec, bool add_newline) {
    std::vector<std::string> terms = split_string(line, "\t");
    const size_t terms_size = terms.size();
    if (terms_size < VCF_REQUIRED_COL_COUNT) {
        throw VcfValidationError("VCF data line did not contain at least 8 terms");
    }
    // get ref to terms (no copy)
    const std::string& ref_name = terms[0];
    const std::string& position = terms[1];
    const std::string& id = terms[2];
    const std::string& ref_bases = terms[3];
    const std::string& alt_bases = terms[4];
    const std::string& quality = terms[5];
    const std::string& filter = terms[6];
    const std::string& info = terms[7];


    // store non-sample columns uncompressed
    // maximum of 9 of these, so fits in a SAMPLE_MASKED_UNCOMPRESSED (max 31)
    // uncompressed flag
    // byte_t non_sample_uncompressed_flag = SAMPLE_MASKED_UNCOMPRESSED | 8;
    // byte_vec.push_back(non_sample_uncompressed_flag); // must update this later

    // TODO estimate the length required for this line so that we can shrink the bytes
    // required for storing length while making it less likely
    // to need to change the number of bytes later, avoiding a vector shift
    LineLengthHeader line_length_header;
    line_length_header.set_extension_count(3);
    uint8_t length_header_bytes[4] = {0xC0, 0, 0, 0};
    line_length_header.deserialize(length_header_bytes);
    byte_vec.push_back(length_header_bytes[0]);
    byte_vec.push_back(length_header_bytes[1]);
    byte_vec.push_back(length_header_bytes[2]);
    byte_vec.push_back(length_header_bytes[3]);

    // include another header for the length of the uncompressed region before the sample columns
    LineLengthHeader required_col_length_header;
    required_col_length_header.set_extension_count(3); // restrict to 3, max = (2 ^ (30) - 1)
    uint8_t required_col_length_header_bytes[4] = {0xC0, 0, 0, 0};
    line_length_header.deserialize(required_col_length_header_bytes);
    byte_vec.push_back(required_col_length_header_bytes[0]);
    byte_vec.push_back(required_col_length_header_bytes[1]);
    byte_vec.push_back(required_col_length_header_bytes[2]);
    byte_vec.push_back(required_col_length_header_bytes[3]);

    push_string_to_byte_vector(byte_vec, ref_name);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, position);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, id);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, ref_bases);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, alt_bases);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, quality);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, filter);
    byte_vec.push_back('\t');
    push_string_to_byte_vector(byte_vec, info);

    //debugf("quality: %s, filter: %s\n", quality.c_str(), filter.c_str());

    debugf("reference_name = %s, pos = %s\n", ref_name.c_str(), position.c_str());
    uint64_t required_length = 7 + ref_name.size() + position.size() + id.size() +
        ref_bases.size() + alt_bases.size() + quality.size() + filter.size() + info.size();

    // handle sample columns
    // first is the FORMAT column
    // TODO parse this and use it to split up sample columns
    if (terms.size() > VCF_REQUIRED_COL_COUNT) {
        // update uncompressed column count
        byte_vec[0] = byte_vec[0] + 1; // no overflow check, max=31, we're only at 9

        const std::string& format = terms[8];
        byte_vec.push_back('\t');
        push_string_to_byte_vector(byte_vec, format);
        debugf("pushing format: %s\n", format.c_str());
        //byte_vec.push_back('\t');
        required_length += format.size() + 1;
    }

    const size_t vcf_sample_start_colum = VCF_REQUIRED_COL_COUNT + 1;
    size_t samples_num = terms_size - (vcf_sample_start_colum);
    if (samples_num > 0) {
        byte_vec.push_back('\t');
        required_length += 1;
    }

    debugf("Updating required length to %lu\n", required_length);
    uint32_t required_length32 = (uint32_t) required_length;
    byte_vec[4] = ((required_length32 >> 24) & 0xFF) | 0xC0;
    byte_vec[5] = (required_length32 >> 16) & 0xFF;
    byte_vec[6] = (required_length32 >> 8) & 0xFF;
    byte_vec[7] = (required_length32 >> 0) & 0xFF;
    debugf("Required length header bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n", byte_vec[4], byte_vec[5], byte_vec[6], byte_vec[7]);

    std::vector<std::string> samples; // copy of sample data

    debugf("terms: %ld, samples: %ld\n", terms_size, samples_num);
    //debugf("samples in row: %ld\n", samples_num);
    samples.resize(samples_num);
    for (size_t i = terms_size; i > vcf_sample_start_colum; i--) {
        size_t idx = i - vcf_sample_start_colum - 1;
        std::string& val = terms.back();
        //debugf("copying samples idx: %ld, val: %s\n", idx, val.c_str());
        samples[idx] = val;
        terms.pop_back();
    }
    #ifdef DEBUG
    debugf("SAMPLES: ");
    for (size_t i = 0; i < samples.size(); i++) {
        debugf("%s ", samples[i].c_str());
    }
    debugf("\n");
    #endif


    for (size_t i = 0; i < samples.size(); i++) {
        const std::string& sample_val = samples.at(i);
        uint8_t max_dedup_00 = 0x7F; // first bit 0 means this is a 0|0 term, 7 bits left for count (max = 127)
        uint8_t max_dedup_01_10_11 = 0x1F; // first 3 bits reserved, 5 bits left for count (max = 31)
        debugf("sample_val: %s\n", sample_val.c_str());
        if (sample_val == "0|0") {
            size_t count = 1;
            i++;
            for ( ;
                    count < max_dedup_00
                    && i < samples.size()
                    && samples.at(i) == "0|0";
                    i++) {
                count++;
            }
            // loop goes to first element not matching conditions, so set i back to that one
            i--;
            debugf("0|0 occurred %ld times\n", count);
            byte_t b = 0x0 | ((uint8_t) count);
            debugf("compressed to: %s\n", char_to_bin_string(b).c_str());
            byte_vec.push_back(b);
        } else if (sample_val == "0|1" || sample_val == "1|0" || sample_val == "1|1") {
            // not 0|0
            size_t count = 1;
            i++;
            for ( ;
                    count < max_dedup_01_10_11
                    && i < samples.size()
                    && samples.at(i) == sample_val;
                    i++ ) {
                count++;
            }
            // loop goes to first element not matching conditions, so set i back to that one
            i--;
            debugf("%s occurred %ld times\n", sample_val.c_str(), count);
            byte_t b;
            if (sample_val == "0|1") {
                b = SAMPLE_MASKED_01 | ((uint8_t) count);
            } else if (sample_val == "1|0") {
                b = SAMPLE_MASKED_10 | ((uint8_t) count);
            } else if (sample_val == "1|1") {
                b = SAMPLE_MASKED_11 | ((uint8_t) count);
            } else {
                throw std::runtime_error("Invalid state! Unknown sample " + sample_val);
            }
            debugf("compressed to: %s\n", char_to_bin_string(b).c_str());
            byte_vec.push_back(b);
        } else {
            // this sample's allele genotype was higher than ALT 1 (>= 2)
            // since this is fairly rare, by the VCF definition, don't bother compressing
            // TODO make this do a lookahead to see if there are multiple uncompressed columns
            // so the count here would not always just be 1
            debugf("sample > 1 (%s), skipping compression\n", sample_val.c_str());
            // send a flag to note this column is *not* compressed
            uint8_t uc_val = SAMPLE_MASKED_UNCOMPRESSED | 1;
            debugf("pushing SAMPLE_MASKED_UNCOMPRESSED + count: %s\n", char_to_bin_string(uc_val).c_str());
            byte_vec.push_back(uc_val);
            push_string_to_byte_vector(byte_vec, sample_val);
            if (i < samples.size() - 1) { // not the last sample
                byte_vec.push_back('\t');
            }
        }
    }

    if (add_newline) {
        byte_vec.push_back('\n');
    }

    // Update the start-of-line line-length header
    // include the required column length header by only subtracting 4 instead of 8
    uint32_t line_length = byte_vec.size() - 4;
    debugf("Updating line length to %u\n", line_length);
    byte_vec[0] = ((line_length >> 24) & 0xFF) | 0xC0;
    byte_vec[1] = (line_length >> 16) & 0xFF;
    byte_vec[2] = (line_length >> 8) & 0xFF;
    byte_vec[3] = (line_length >> 0) & 0xFF;
    debugf("Line length header bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n", byte_vec[0], byte_vec[1], byte_vec[2], byte_vec[3]);

    return 0;
}

int compress(const std::string& input_filename, const std::string& output_filename) {
    std::ifstream input_fstream(input_filename);
    //size_t readbuf_size = 10 * 1024 * 1024; // 1 MiB
    //char *local_readbuf = new char[readbuf_size];
    //input_fstream.rdbuf()->pubsetbuf(local_readbuf, readbuf_size);
    std::ofstream output_fstream(output_filename);
    //VcfLineStateMachine lineStateMachine;
    std::string linebuf;
    VcfCompressionSchema schema;
    size_t variant_count = 0;
    std::vector<byte_t> compressed_line;
    compressed_line.reserve(4096);

    while (std::getline(input_fstream, linebuf)) {
        if (linebuf.size() == 0) {
            // empty input line, ignore
            continue;
        } else if (linebuf[0] == '#' && linebuf[1] == '#' /*linebuf.substr(0, 2) == "##"*/) {
            //lineStateMachine.to_meta();
            // compress vcf header
            // TODO
            output_fstream << linebuf << "\n";
        } else if (linebuf[0] == '#' /*linebuf.substr(0, 1) == "#"*/) {
            //lineStateMachine.to_header();
            // get the number of samples from the header
            std::vector<std::string> line_terms = split_string(linebuf, "\t");
            if (line_terms.size() < VCF_REQUIRED_COL_COUNT) {
                //delete local_readbuf;
                throw VcfValidationError("VCF Header did not have enough columns");
            }
            schema.sample_count = line_terms.size() - VCF_REQUIRED_COL_COUNT - 1;
            debugf("sample count: %ld\n", schema.sample_count);
            // insert header in raw format
            output_fstream << linebuf << "\n";
        } else {
            // treat line as variant
            variant_count++;
            //lineStateMachine.to_variant();
            compressed_line.clear();
            /*int status = */compress_data_line(linebuf, schema, compressed_line, true);
            if (compressed_line.back() != '\n') {
                throw std::runtime_error("No newline at end of compressed line!");
            }
            for (std::vector<byte_t>::iterator iter = compressed_line.begin(); iter != compressed_line.end(); iter++) {
                output_fstream.write((const char*)(&(*iter)), 1);
            }
            //output_fstream.write("\n", 1);
        }
    }
    debugf("variant count: %ld\n", variant_count);
    //delete local_readbuf;
    return 0;
}


/**
 * Read the line length headers from the start of a compressed data line, store in `length_headers`.
 *
 * On success returns a positive integer indicating the number of bytes read.
 *
 * If EOF, returns 0.
 *
 * If negative, is the error return status from read(2).
 */
int read_compressed_line_length_headers(int input_fd, struct compressed_line_length_headers *length_headers) {
    int status;
    // int total_bytes = 0;

    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    uint8_t line_length_header_bytes[8] = {0,0,0,0,0,0,0,0};
    status = read(input_fd, &line_length_header_bytes, 8);
    if (status != 8) {
        if (status < 0) {
            return status;
        } else if (status == 0) {
            debugf("Finished querying file\n");
            return 0;
        } else {
            return status;
        }
    }

    uint8_t *required_columns_length_header_bytes = line_length_header_bytes + 4;

    // uint8_t line_length_header_bytes[4] = {0,0,0,0};
    // status = read(input_fd, &line_length_header_bytes, 4);
    // if (status != 4) {
    //     if (status < 0) {
    //         return status;
    //         // throw std::runtime_error("Error reading from file");
    //     } else if (status == 0) {
    //         debugf("Finished querying file\n");
    //         return 0;
    //     } else if (status < 4) {
    //         return status;
    //         // throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
    //     }
    // }

    debugf("line_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
            line_length_header_bytes[0], line_length_header_bytes[1],
            line_length_header_bytes[2], line_length_header_bytes[3]);

    // uint8_t required_columns_length_header_bytes[4] = {0,0,0,0};
    // status = read(input_fd, &required_columns_length_header_bytes, 4);
    // if (status != 4) {
    //     if (status < 0) {
    //         return status;
    //     } else if (status < 4) {
    //         // throw std::runtime_error(string_format("Only read %d bytes, expected 4", status));
    //         return 4 + status;
    //     }
    // }

    debugf("required_columns_length_header_bytes: 0x%02X 0x%02X 0x%02X 0x%02X\n",
            required_columns_length_header_bytes[0], required_columns_length_header_bytes[1],
            required_columns_length_header_bytes[2], required_columns_length_header_bytes[3]);

    int read_bytes = 4 + 4; // length headers

    debugf("After length headers, stream positioned so next byte is at position %ld (0x%08lx)\n",
            tellfd(input_fd),
            tellfd(input_fd));

    LineLengthHeader line_length_header;
    // line_length_header.set_extension_count(3); // TODO interpret

    line_length_header.deserialize(line_length_header_bytes);
    length_headers->line_length = line_length_header.length;

    line_length_header.deserialize(required_columns_length_header_bytes);
    length_headers->required_columns_length = line_length_header.length;

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("read_compressed_line_length_headers time: %lu\n", duration.count());
    #endif

    return read_bytes;
}


int decompress2_data_line(
        std::ifstream& input_fstream,
        const VcfCompressionSchema& schema,
        // std::string& linebuf,
        string_t *linebuf,
        size_t *compressed_line_length) {
    debugf("%s decompressing line, expecting %lu samples\n", __FUNCTION__, schema.sample_count);
    // decompress a single variant line
    int ib;
    unsigned char b;

    // keep track of how many columns we've seen
    size_t line_byte_count = 0;
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;

    // interpret first 1 byte as a skip flag, then up to 3 additional bytes
    ib = input_fstream.get();
    debugf("ib = 0x%08X\n", ib);
    if (ib == eof) {
        debugf("%s, no data in input_fstream\n", __FUNCTION__);
        return -1;
    }
    line_byte_count++;
    b = 0xFF & ib;

    // read the line length header
    uint8_t line_length_bytes[4] = {b, 0, 0, 0};
    LineLengthHeader line_length_header;
    line_length_header.deserialize(line_length_bytes);
    uint8_t extension_count = line_length_header.extension_count;
    debugf("Line length extension count: %u\n", extension_count);
    for (uint8_t i = 1; i <= extension_count; i++) {
        int ib = input_fstream.get();
        if (ib == eof) {
            debugf("%s, no data in input_fstream\n", __FUNCTION__);
            return -1;
        }
        // line_byte_count++;
        line_length_bytes[i] = 0xFF & ib;
    }
    line_byte_count += extension_count;
    line_length_header.deserialize(line_length_bytes);


    // read the required column skip length
    ib = input_fstream.get();
    if (ib == eof) {
        debugf("%s, no data in input_fstream\n", __FUNCTION__);
        return -1;
    }
    line_byte_count++;
    b = 0xFF & ib;
    uint8_t required_length_bytes[4] = {b, 0, 0, 0};
    LineLengthHeader required_length_header;
    required_length_header.deserialize(required_length_bytes);
    uint8_t required_extension_count = required_length_header.extension_count;
    for (uint8_t i = 0; i < required_extension_count; i++) {
        int ib = input_fstream.get();
        if (ib == eof) {
            debugf("%s, no data in input_fstream\n", __FUNCTION__);
            return -1;
        }
        //line_byte_count++;
        required_length_bytes[i+1] = 0xFF & ib;
    }
    line_byte_count += required_extension_count;
    debugf("Deserializing required length header\n");
    required_length_header.deserialize(required_length_bytes);
    uint32_t required_length = required_length_header.length;
    debugf("Skipping %u bytes for required columns section\n", required_length);

    for (size_t i = 0; i < required_length; i++) {
        input_fstream.read((char*)&b, 1);
        char cb = reinterpret_cast<char&>(b);
        // linebuf.push_back(cb);
        string_appendc(linebuf, cb);
        if (cb == '\t') {
            line_tab_count++;
        }
    }
    debugf("Finished reading required columns: %s\n", linebuf->buf);
    line_byte_count += required_length;


    // check to ensure we read in the appropriate number of uncompressed columns
    // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
    if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
        if (line_tab_count == VCF_REQUIRED_COL_COUNT && schema.sample_count == 0) {
            // do nothing
        } else {
            debugf("line_tab_count: %lu\n", line_tab_count);
            throw VcfValidationError("Did not read all uncompressed columns");
            return 0;
        }
    }

    debugf("Reading sample columns\n");
    // read the sample columns
    while (line_sample_count < schema.sample_count) {
        ib = input_fstream.peek();
        b = 0xFF & ib;
        //debugf("b: %02x\n", (unsigned char)b);

        if (ib == eof /*|| (char)ib == '\n' || (char)ib == '\t'*/) {
            std::ostringstream msg;
            msg << "Missing samples, expected " << schema.sample_count
                << ", received " << line_sample_count;
            throw VcfValidationError(msg.str().c_str());
        }
        input_fstream.get(); // remove peeked char from stream
        line_byte_count++;

        if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            uint8_t counter = count;
            debugf("0|0 repeat count: %u\n", count);

            while (counter--) {
                // linebuf.append(GT_00.c_str());
                // linebuf.append(tab);
                string_appends(linebuf, GT_00.c_str());
                string_appendc(linebuf, tab[0]);
            }

            // update line counts
            line_tab_count += count;
            line_sample_count += count;

            // remove last tab if at end of line
            if (line_sample_count >= schema.sample_count) {
                // linebuf.pop_back();
                string_pop(linebuf);
                line_tab_count--;
            }
        } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0; // number of uncompressed columns
            while (ucounter < uncompressed_count) {
                input_fstream.read((char*)(&b), 1);
                line_byte_count++;
                // TODO handle newline ?
                if (b == '\n') {
                    // newline after uncompressed sample value
                    ucounter++;
                    line_sample_count++;
                    if (ucounter != uncompressed_count) {
                        throw VcfValidationError("Reached end of line before reading all decompressed columns");
                    }
                    debugf("Putting 0x%02X back into input stream\n", b);
                    input_fstream.putback(b);
                }
                else if (b == '\t') {
                    // don't push tabs, handled outside if
                    ucounter++;
                    line_tab_count++;
                    line_sample_count++;
                    if (line_sample_count < schema.sample_count) {
                        // if not the last term, include the tab
                        // otherwise, the tab is handled outside the if
                        // linebuf.push_back(b);
                        string_appendc(linebuf, b);
                    }
                } else {
                    // linebuf.push_back(b);
                    string_appendc(linebuf, b);
                }
            }
        } else {
            // either 0|1, 1|0, or 1|1
            byte_t masked = b & SAMPLE_MASK_01_10_11;
            const std::string *sample_str = NULL;

            if (masked == SAMPLE_MASKED_01) {
                sample_str = &GT_01;
            } else if (masked == SAMPLE_MASKED_10) {
                sample_str = &GT_10;
            } else if (masked == SAMPLE_MASKED_11) {
                sample_str = &GT_11;
            } else {
                throw std::runtime_error("Error during decompression of compressed file, unrecognized bitmask");
            }

            // write the sample GT and increment counters
            uint8_t count = (b & (~SAMPLE_MASK_01_10_11)) & 0xFF;
            debugf("Got %s, count: %u\n", sample_str->c_str(), count);
            while (count--) {
                // linebuf.append(*sample_str);
                string_appends(linebuf, sample_str->c_str());
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    // linebuf.append(tab);
                    string_appendc(linebuf, tab[0]);
                }
                line_tab_count++;
            }
        } // end flag cases
    } // end sample loop
    debugf("Finished reading samples\n");

    // make sure next byte is a newline
    ib = input_fstream.get();
    line_byte_count++;
    if (ib == '\n') {
        b = (char)ib;
        // linebuf.push_back(b);
        string_appendc(linebuf, b);
    } else {
        throw VcfValidationError("Sample line did not end in a newline\n");
    }

    *compressed_line_length = line_byte_count;

    return 0;
}


/**
 * Reads from input_fstream to decompress one line from the vcfc file.
 * Schema must match the actual schema of the file.
 *
 * NOTE: Appends the decompressed line to linebuf. Usually caller
 * should send an empty linebuf.
 *
 * If linebuf needs to be reallocated to a larger buffer,
 * the value of linebuf and linebuf_capacity will be updated.
 *
 * Returns 0 on EOF, 1 on success, -1 on error.
 */
// int decompress2_data_line_fd2(
//         int input_fd,
//         const VcfCompressionSchema& schema,
//         string_t *linebuf,
//         size_t *compressed_line_length) {
//     debugf("%s decompressing line, expecting %lu samples\n", __FUNCTION__, schema.sample_count);
//     // decompress a single variant line
//     unsigned char b;
//     #ifdef TIMING
//     std::chrono::time_point<std::chrono::steady_clock> start;
//     std::chrono::time_point<std::chrono::steady_clock> end;
//     std::chrono::nanoseconds duration;
//     #endif

//     // keep track of how many columns we've seen
//     size_t line_byte_count = 0;
//     size_t line_tab_count = 0;
//     size_t line_sample_count = 0;

//     if (eof_fd(input_fd)) {
//         debugf("%s, no data in input_fd\n", __FUNCTION__);
//         return 0;
//     }
//     struct compressed_line_length_headers line_length_headers;
//     memset(&line_length_headers, 0, sizeof(struct compressed_line_length_headers));

//     read_compressed_line_length_headers(input_fd, &line_length_headers);
//     line_byte_count += compressed_line_length_headers_size;

//     uint32_t required_length = line_length_headers.required_columns_length;
//     debugf("Skipping %u bytes for required columns section\n", required_length);

//     for (size_t i = 0; i < required_length; i++) {
//         char cb;
//         if (read(input_fd, &cb, 1) <= 0) {
//             throw std::runtime_error("fread error");
//         }

//         string_appendc(linebuf, cb);

//         if (cb == '\t') {
//             line_tab_count++;
//         }
//     }
//     debugf("Finished reading required columns\n");
//     line_byte_count += required_length;


//     // check to ensure we read in the appropriate number of uncompressed columns
//     // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
//     if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
//         if (line_tab_count == VCF_REQUIRED_COL_COUNT && schema.sample_count == 0) {
//             // do nothing
//         } else {
//             debugf("line_tab_count: %lu\n", line_tab_count);
//             throw VcfValidationError("Did not read all uncompressed columns");
//             // return -1;
//         }
//     }

//     #ifdef TIMING
//     start = std::chrono::steady_clock::now();
//     #endif

//     debugf("Reading sample columns\n");
//     // read the sample columns
//     while (line_sample_count < schema.sample_count) {
//         //debugf("Trying to read a sample column\n");
//         if (read(input_fd, &b, 1) <= 0) {
//             std::ostringstream msg;
//             msg << "Missing samples, expected " << schema.sample_count
//                 << ", received " << line_sample_count;
//             throw VcfValidationError(msg.str().c_str());
//         }
//         line_byte_count++;

//         if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
//             // is a 0|0 column
//             uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
//             uint8_t counter = count;
//             debugf("0|0 repeat count: %u\n", count);

//             while (counter--) {
//                 // linebuf.append(GT_00.c_str());
//                 // linebuf.append(tab);
//                 string_appends(linebuf, GT_00.c_str());
//                 string_appends(linebuf, tab);
//             }

//             // update line counts
//             line_tab_count += count;
//             line_sample_count += count;

//             // remove last tab if at end of line
//             if (line_sample_count >= schema.sample_count) {
//                 // linebuf.pop_back();
//                 string_pop(linebuf);
//                 line_tab_count--;
//             }
//         } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
//             uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
//             debugf("%u uncompressed columns follow\n", uncompressed_count);
//             // uncompressed columns follow
//             uint8_t ucounter = 0; // number of uncompressed columns
//             while (ucounter < uncompressed_count) {
//                 if (read(input_fd, &b, 1) <= 0) {
//                     throw std::runtime_error("Couldn't read from input_fd");
//                 }
//                 // fread(&b, sizeof(char), 1, input_stream);
//                 line_byte_count++;
//                 // TODO handle newline ?
//                 if (b == '\n') {
//                     // newline after uncompressed sample value
//                     ucounter++;
//                     line_sample_count++;
//                     if (ucounter != uncompressed_count) {
//                         throw VcfValidationError("Reached end of line before reading all decompressed columns");
//                     }
//                     // ending newline handled outside loop
//                     debugf("got ending newline\n");
//                     lseek(input_fd, -1, SEEK_CUR);
//                 }
//                 else if (b == '\t') {
//                     // don't push tabs, handled outside if
//                     ucounter++;
//                     line_tab_count++;
//                     line_sample_count++;
//                     if (line_sample_count < schema.sample_count) {
//                         // if not the last term, include the tab
//                         // otherwise, the tab is handled outside the if
//                         // linebuf.push_back(b);
//                         string_appendc(linebuf, b);
//                     }
//                 } else {
//                     // linebuf.push_back(b);
//                     string_appendc(linebuf, b);
//                 }
//             }
//         } else {
//             // either 0|1, 1|0, or 1|1
//             byte_t masked = b & SAMPLE_MASK_01_10_11;
//             const std::string *sample_str = NULL;

//             if (masked == SAMPLE_MASKED_01) {
//                 sample_str = &GT_01;
//             } else if (masked == SAMPLE_MASKED_10) {
//                 sample_str = &GT_10;
//             } else if (masked == SAMPLE_MASKED_11) {
//                 sample_str = &GT_11;
//             } else {
//                 throw std::runtime_error("Error during decompression of compressed file, unrecognized bitmask");
//             }

//             // write the sample GT and increment counters
//             uint8_t count = (b & (~SAMPLE_MASK_01_10_11)) & 0xFF;
//             debugf("Got %s, count: %u\n", sample_str->c_str(), count);
//             while (count--) {
//                 // linebuf.append(*sample_str);
//                 string_appends(linebuf, sample_str->c_str());
//                 line_sample_count++;
//                 if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
//                     // linebuf.append(tab);
//                     string_appendc(linebuf, tab[0]);
//                 }
//                 line_tab_count++;
//             }
//         } // end flag cases
//     } // end sample loop
//     debugf("Finished reading samples\n");


//     #ifdef TIMING
//     end = std::chrono::steady_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
//     printf("decompress-sample-loop time: %lu\n", duration.count());
//     #endif

//     if (read(input_fd, &b, 1) <= 0) {
//         throw std::runtime_error("Failed to read line ending");
//     }
//     if (b == '\n') {
//         // linebuf.push_back('\n');
//         string_appendc(linebuf, '\n');
//     } else {
//         throw VcfValidationError("Sample line did not end in a newline\n");
//     }

//     *compressed_line_length = line_byte_count;

//     return 1;
// }


/**
 * Memoize decompressed strings
 */
// Memoize instances of runs already decompressed
// map sample identifier to a lookup table of run length, corresponding to decompressed strings
// std::unordered_map<std::string,std::unordered_map<uint8_t,std::string>> decompressed_cache;

std::map<std::string,std::array<std::string,UINT8_MAX>> decompressed_cache;

const std::string generate_cache_line(const std::string& sample_value, const uint8_t run_length) {
    uint8_t counter = run_length;
    debugf("Caching decompressed sample = %s, run_length = %d\n", sample_value.c_str(), run_length);
    std::string cache_line;
    cache_line.reserve(
        (sample_value.size() + 1) * run_length
    );
    while (counter--) {
        cache_line.append(sample_value);
        cache_line.push_back(tab[0]);
    }
    // decompressed_cache[sample_value][run_length] = cache_line;
    // debugf("stored in cached as: %s\n", decompressed_cache.at(sample_value).at(run_length).c_str());
    // return decompressed_cache[sample_value][run_length];
    return cache_line;
};

const std::string get_or_set_decompressed_cache(
        const std::string& sample_value,
        const uint8_t run_length) {
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif
    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    // size_t map_sample_count = decompressed_cache.count(sample_value);
    bool map_sample_exists = decompressed_cache.find(sample_value) != decompressed_cache.end();

    if (!map_sample_exists) {
        std::array<std::string,UINT8_MAX> run_vector;
        const std::string cache_line = generate_cache_line(sample_value, run_length);
        run_vector[run_length] = cache_line;
        decompressed_cache[sample_value] = run_vector;
        debugf("stored in cached as: %s\n", decompressed_cache.at(sample_value).at(run_length).c_str());
        #ifdef TIMING
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        printf("get_or_set_decompressed_cache time: %lu\n", duration.count());
        #endif
        return cache_line;
    }
    else {
        std::array<std::string,UINT8_MAX>& run_vector = decompressed_cache.at(sample_value);

        if (run_vector[run_length].size() == 0) {
            const std::string& cache_line = generate_cache_line(sample_value, run_length);
            decompressed_cache[sample_value][run_length] = cache_line;
            debugf("stored in cached as: %s\n", decompressed_cache.at(sample_value).at(run_length).c_str());
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("get_or_set_decompressed_cache time: %lu\n", duration.count());
            #endif
            return cache_line;
        } else {
            #ifdef TIMING
            end = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            printf("get_or_set_decompressed_cache time: %lu\n", duration.count());
            #endif
            return run_vector[run_length];
        }
    }
    // std::string& val = decompressed_cache[sample_value][run_length];
    // return val;
}

/**
 * Reads from input_fstream to decompress one line from the vcfc file.
 * Schema must match the actual schema of the file.
 *
 * NOTE: Appends the decompressed line to linebuf. Usually caller
 * should send an empty linebuf.
 *
 * If linebuf needs to be reallocated to a larger buffer,
 * the value of linebuf and linebuf_capacity will be updated.
 *
 * Returns 0 on EOF, 1 on success, -1 on error.
 */
int decompress2_data_line_fd2_string(
        int input_fd,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length) {
    debugf("%s decompressing line, expecting %lu samples\n", __FUNCTION__, schema.sample_count);
    // decompress a single variant line
    unsigned char b;
    int status;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;

    std::chrono::time_point<std::chrono::steady_clock> start2;
    std::chrono::time_point<std::chrono::steady_clock> end2;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    // keep track of how many columns we've seen
    size_t line_byte_count = 0;
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;

    if (eof_fd(input_fd)) {
        debugf("%s, no data in input_fd\n", __FUNCTION__);
        return 0;
    }
    struct compressed_line_length_headers line_length_headers;
    memset(&line_length_headers, 0, sizeof(struct compressed_line_length_headers));

    read_compressed_line_length_headers(input_fd, &line_length_headers);
    line_byte_count += compressed_line_length_headers_size;

    uint32_t required_length = line_length_headers.required_columns_length;
    debugf("Skipping %u bytes for required columns section\n", required_length);

    #ifdef TIMING
    start2 = std::chrono::steady_clock::now();
    #endif
    // 19 seconds
    // for (size_t i = 0; i < required_length; i++) {
    //     char cb;
    //     if (read(input_fd, &cb, 1) <= 0) {
    //         throw std::runtime_error("fread error");
    //     }

    //     // string_appendc(linebuf, cb);
    //     linebuf.push_back(cb);

    //     if (cb == '\t') {
    //         line_tab_count++;
    //     }
    // }

    // 175 milliseconds
    // read the skipped bytes
    char *buf = (char*) calloc(required_length + 1, sizeof(char));
    debugf("Reading required columns\n");
    ssize_t read_n = read(input_fd, buf, (size_t)required_length);
    if (read_n <= 0) {
        throw std::runtime_error("error reading required columns");
    }
    if (read_n < required_length) {
        debugf("While reading required columns expected %u bytes but got %ld\n", required_length, read_n);
        throw std::runtime_error("");
    }

    // count the tabs in the read string
    debugf("Counting tabs\n");
    for (uint32_t i = 0; i < required_length; i++) {
        if (buf[i] == '\t') {
            line_tab_count++;
        }
    }
    linebuf.append(buf);

    #ifdef TIMING
    end2 = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2);
    printf("skipping-required-columns time: %lu\n", duration.count());
    #endif
    debugf("Finished reading required columns\n");
    line_byte_count += required_length;


    // check to ensure we read in the appropriate number of uncompressed columns
    // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
    if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
        if (line_tab_count == VCF_REQUIRED_COL_COUNT && schema.sample_count == 0) {
            // do nothing
        } else {
            debugf("line_tab_count: %lu\n", line_tab_count);
            throw VcfValidationError("Did not read all uncompressed columns");
            // return -1;
        }
    }

    debugf("Reading sample columns\n");
    // read the sample columns
    while (line_sample_count < schema.sample_count) {
        //debugf("Trying to read a sample column\n");
        if (read(input_fd, &b, 1) <= 0) {
            std::ostringstream msg;
            msg << "Missing samples, expected " << schema.sample_count
                << ", received " << line_sample_count;
            throw VcfValidationError(msg.str().c_str());
        }
        line_byte_count++;

        if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            // uint8_t counter = count;
            debugf("0|0 repeat count: %u\n", count);

            const std::string& cache_line = get_or_set_decompressed_cache(GT_00, count);
            linebuf.append(cache_line);

            // update line counts
            line_tab_count += count;
            line_sample_count += count;

            // remove last tab if at end of line
            if (line_sample_count >= schema.sample_count) {
                linebuf.pop_back();
                // string_pop(linebuf);
                line_tab_count--;
            }
        } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0; // number of uncompressed columns
            while (ucounter < uncompressed_count) {
                if (read(input_fd, &b, 1) <= 0) {
                    throw std::runtime_error("Couldn't read from input_fd");
                }
                // fread(&b, sizeof(char), 1, input_stream);
                line_byte_count++;
                // TODO handle newline ?
                if (b == '\n') {
                    // newline after uncompressed sample value
                    ucounter++;
                    line_sample_count++;
                    if (ucounter != uncompressed_count) {
                        throw VcfValidationError("Reached end of line before reading all decompressed columns");
                    }
                    // ending newline handled outside loop
                    debugf("got ending newline\n");
                    lseek(input_fd, -1, SEEK_CUR);
                }
                else if (b == '\t') {
                    // don't push tabs, handled outside if
                    ucounter++;
                    line_tab_count++;
                    line_sample_count++;
                    if (line_sample_count < schema.sample_count) {
                        // if not the last term, include the tab
                        // otherwise, the tab is handled outside the if
                        linebuf.push_back(b);
                        // string_appendc(linebuf, b);
                    }
                } else {
                    linebuf.push_back(b);
                    // string_appendc(linebuf, b);
                }
            }
        } else {
            // either 0|1, 1|0, or 1|1
            byte_t masked = b & SAMPLE_MASK_01_10_11;
            const std::string *sample_str = NULL;

            if (masked == SAMPLE_MASKED_01) {
                sample_str = &GT_01;
            } else if (masked == SAMPLE_MASKED_10) {
                sample_str = &GT_10;
            } else if (masked == SAMPLE_MASKED_11) {
                sample_str = &GT_11;
            } else {
                throw std::runtime_error("Error during decompression of compressed file, unrecognized bitmask");
            }

            // write the sample GT and increment counters
            uint8_t count = (b & (~SAMPLE_MASK_01_10_11)) & 0xFF;
            debugf("Got %s, count: %u\n", sample_str->c_str(), count);
            while (count--) {
                linebuf.append(*sample_str);
                // string_appends(linebuf, sample_str->c_str());
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    linebuf.push_back(tab[0]);
                    // string_appendc(linebuf, tab[0]);
                }
                line_tab_count++;
            }
        } // end flag cases
    } // end sample loop
    debugf("Finished reading samples\n");




    if (read(input_fd, &b, 1) <= 0) {
        throw std::runtime_error("Failed to read line ending");
    }
    if (b == '\n') {
        linebuf.push_back(b);
        // string_appendc(linebuf, '\n');
    } else {
        throw VcfValidationError("Sample line did not end in a newline\n");
    }

    *compressed_line_length = line_byte_count;

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("decompress2_data_line_fd2_string time: %lu\n", duration.count());
    #endif
    return 1;
}


int decompress2_data_line_FILE(
        int input_fd,
        const VcfCompressionSchema& schema,
        std::string& linebuf,
        size_t *compressed_line_length) {
    debugf("%s decompressing line, expecting %lu samples\n", __FUNCTION__, schema.sample_count);
    // decompress a single variant line
    unsigned char b;
    int status;
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;

    std::chrono::time_point<std::chrono::steady_clock> start2;
    std::chrono::time_point<std::chrono::steady_clock> end2;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    // Create FILE* wrapper of input_fd to improve small read performance.
    // input_fd_dup is automatically closed when input_file is closed.
    int input_fd_dup = dup(input_fd);
    FILE *input_file = fdopen(input_fd_dup, "r");


    // keep track of how many columns we've seen
    size_t line_byte_count = 0;
    size_t line_tab_count = 0;
    size_t line_sample_count = 0;

    if (eof_fd(input_fd)) {
        debugf("%s, no data in input_fd\n", __FUNCTION__);
        return 0;
    }
    struct compressed_line_length_headers line_length_headers;
    memset(&line_length_headers, 0, sizeof(struct compressed_line_length_headers));

    read_compressed_line_length_headers(input_fd, &line_length_headers);
    line_byte_count += compressed_line_length_headers_size;

    uint32_t required_length = line_length_headers.required_columns_length;
    debugf("Skipping %u bytes for required columns section\n", required_length);

    #ifdef TIMING
    start2 = std::chrono::steady_clock::now();
    #endif
    // 19 seconds
    // for (size_t i = 0; i < required_length; i++) {
    //     char cb;
    //     if (read(input_fd, &cb, 1) <= 0) {
    //         throw std::runtime_error("fread error");
    //     }

    //     // string_appendc(linebuf, cb);
    //     linebuf.push_back(cb);

    //     if (cb == '\t') {
    //         line_tab_count++;
    //     }
    // }

    // 175 milliseconds
    // read the skipped bytes
    char *buf = (char*) calloc(required_length + 1, sizeof(char));
    debugf("Reading required columns\n");
    // ssize_t read_n = read(input_fd, buf, (size_t)required_length);
    size_t read_n = fread(buf, sizeof(char), required_length, input_file);
    if (read_n <= 0) {
        throw std::runtime_error("error reading required columns");
    }
    if (read_n < required_length) {
        debugf("While reading required columns expected %u bytes but got %ld\n", required_length, read_n);
        throw std::runtime_error("");
    }

    // count the tabs in the read string
    debugf("Counting tabs\n");
    for (uint32_t i = 0; i < required_length; i++) {
        if (buf[i] == '\t') {
            line_tab_count++;
        }
    }
    linebuf.append(buf);

    #ifdef TIMING
    end2 = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2);
    printf("skipping-required-columns time: %lu\n", duration.count());
    #endif
    debugf("Finished reading required columns\n");
    line_byte_count += required_length;


    // check to ensure we read in the appropriate number of uncompressed columns
    // here it expects VCF_REQUIRED_COL_COUNT + 1 because it skips the format column as well
    if (line_tab_count != VCF_REQUIRED_COL_COUNT + 1) {
        if (line_tab_count == VCF_REQUIRED_COL_COUNT && schema.sample_count == 0) {
            // do nothing
        } else {
            debugf("line_tab_count: %lu\n", line_tab_count);
            throw VcfValidationError("Did not read all uncompressed columns");
            // return -1;
        }
    }

    debugf("Reading sample columns\n");
    // read the sample columns
    while (line_sample_count < schema.sample_count) {
        //debugf("Trying to read a sample column\n");
        // if (read(input_fd, &b, 1) <= 0) {
        if (fread(&b, 1, 1, input_file) < 1) {
            std::ostringstream msg;
            msg << "Missing samples, expected " << schema.sample_count
                << ", received " << line_sample_count;
            throw VcfValidationError(msg.str().c_str());
        }
        line_byte_count++;

        if ((b & SAMPLE_MASK_00) == SAMPLE_MASKED_00) {
            // is a 0|0 column
            uint8_t count = (b & (~SAMPLE_MASK_00)) & 0xFF; // and with inverse of flag mask
            debugf("0|0 repeat count: %u\n", count);

            // uint8_t counter = count;
            // while (counter--) {
            //     linebuf.append(GT_00);
            //     linebuf.append(tab);
            //     // string_appends(linebuf, GT_00.c_str());
            //     // string_appendc(linebuf, tab[0]);
            // }


            const std::string& cache_line = get_or_set_decompressed_cache(GT_00, count);
            linebuf.append(cache_line);

            // update line counts
            line_tab_count += count;
            line_sample_count += count;

            // remove last tab if at end of line
            if (line_sample_count >= schema.sample_count) {
                linebuf.pop_back();
                // string_pop(linebuf);
                line_tab_count--;
            }
        } else if ((b & SAMPLE_MASK_UNCOMPRESSED) == SAMPLE_MASKED_UNCOMPRESSED) {
            uint8_t uncompressed_count = b & (~SAMPLE_MASK_UNCOMPRESSED);
            debugf("%u uncompressed columns follow\n", uncompressed_count);
            // uncompressed columns follow
            uint8_t ucounter = 0; // number of uncompressed columns
            while (ucounter < uncompressed_count) {
                // if (read(input_fd, &b, 1) <= 0) {
                if (fread(&b, 1, 1, input_file) < 1) {
                    throw std::runtime_error("Couldn't read from input_fd");
                }
                // fread(&b, sizeof(char), 1, input_stream);
                line_byte_count++;
                // TODO handle newline ?
                if (b == '\n') {
                    // newline after uncompressed sample value
                    ucounter++;
                    line_sample_count++;
                    if (ucounter != uncompressed_count) {
                        throw VcfValidationError("Reached end of line before reading all decompressed columns");
                    }
                    // ending newline handled outside loop
                    debugf("got ending newline\n");
                    // lseek(input_fd, -1, SEEK_CUR);
                    fseek(input_file, -1, SEEK_CUR);
                }
                else if (b == '\t') {
                    // don't push tabs, handled outside if
                    ucounter++;
                    line_tab_count++;
                    line_sample_count++;
                    if (line_sample_count < schema.sample_count) {
                        // if not the last term, include the tab
                        // otherwise, the tab is handled outside the if
                        linebuf.push_back(b);
                        // string_appendc(linebuf, b);
                    }
                } else {
                    linebuf.push_back(b);
                    // string_appendc(linebuf, b);
                }
            }
        } else {
            // either 0|1, 1|0, or 1|1
            byte_t masked = b & SAMPLE_MASK_01_10_11;
            const std::string *sample_str = NULL;

            if (masked == SAMPLE_MASKED_01) {
                sample_str = &GT_01;
            } else if (masked == SAMPLE_MASKED_10) {
                sample_str = &GT_10;
            } else if (masked == SAMPLE_MASKED_11) {
                sample_str = &GT_11;
            } else {
                throw std::runtime_error("Error during decompression of compressed file, unrecognized bitmask");
            }

            // write the sample GT and increment counters
            uint8_t count = (b & (~SAMPLE_MASK_01_10_11)) & 0xFF;
            debugf("Got %s, count: %u\n", sample_str->c_str(), count);



            // const std::string& cache_line = get_or_set_decompressed_cache(*sample_str, count);
            // linebuf.append(cache_line);

            // // update line counts
            // line_tab_count += count;
            // line_sample_count += count;

            // // remove last tab if at end of line
            // if (line_sample_count >= schema.sample_count) {
            //     linebuf.pop_back();
            //     line_tab_count--;
            // }



            while (count--) {
                linebuf.append(*sample_str);
                // string_appends(linebuf, sample_str->c_str());
                line_sample_count++;
                if (line_sample_count < schema.sample_count /*input_fstream.peek() != '\n'*/) {
                    linebuf.push_back(tab[0]);
                    // string_appendc(linebuf, tab[0]);
                }
                line_tab_count++;
            }



        } // end flag cases
    } // end sample loop
    debugf("Finished reading samples\n");

    // if (read(input_fd, &b, 1) <= 0) {
    if (fread(&b, 1, 1, input_file) < 1) {
        throw std::runtime_error("Failed to read line ending");
    }
    if (b == '\n') {
        linebuf.push_back(b);
        // string_appendc(linebuf, '\n');
    } else {
        throw VcfValidationError("Sample line did not end in a newline\n");
    }

    *compressed_line_length = line_byte_count;

    debugf("input_fd offset: %ld\n", tellfd(input_fd));
    debugf("input_fd_dup offset: %ld\n", tellfd(input_fd_dup));
    debugf("input_file offset: %ld\n", ftell(input_file));

    // Reset input fd stream to offset of input_file to account for readahead
    lseek(input_fd, ftell(input_file), SEEK_SET);
    debugf("input_fd offset reset to: %ld\n", tellfd(input_fd));


    fclose(input_file);

    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("decompress2_data_line_fd2_string time: %lu\n", duration.count());
    #endif
    return 1;
}


/**
 * Reads from input_fstream. Assumes stream position is in the metadata section.
 * Reads all following metadata lines and the header line. If the stream does not conform
 * to the VCF line specification, throws an error.
 *
 * Places all lines read into output_vector, and the schema into output_schema.
 */
int decompress2_metadata_headers_fd(
        int input_fd,
        std::vector<std::string>& output_vector,
        VcfCompressionSchema& output_schema) {
    #ifdef TIMING
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::nanoseconds duration;
    #endif

    #ifdef TIMING
    start = std::chrono::steady_clock::now();
    #endif

    // decompress all metadata and header lines
    bool got_meta = false, got_header = false;
    //int i1, i2;
    char c1, c2;
    size_t meta_count = 0, header_count = 0;

    std::string linebuf;
    linebuf.reserve(4 * 4096);
    // string_t linebuf;
    // string_reserve(&linebuf, 4096);

    while (true) {
        debugf("Reading next line\n");
        linebuf.clear();
        // string_clear(&linebuf);

        if (read(input_fd, &c1, 1) <= 0 && (!got_header || !got_meta)) {
            throw VcfValidationError("File ended before a header or metadata line");
        }
        debugf("c1: %02x\n", c1);
        if (c1 != '#') {
            if (!got_meta || !got_header) {
                throw VcfValidationError("File was missing headers or metadata");
            }
            debugf("%s", "Finished reading metadata and headers\n");
            // need to undo read
            debugf("Moving file offset back one\n");
            lseek(input_fd, -1, SEEK_CUR);
            break;
        } else if (got_header == true) {
            // got a metadata or header row after already receiving a header row
            throw VcfValidationError("Read a metadata or header row after already reading a header");
        }

        if (read(input_fd, &c2, 1) <= 0) {
            throw VcfValidationError("Invalid format, empty header row");
        }

        if (c2 == '#') {
            if (got_header) {
                throw VcfValidationError("Got a metadata row after the CSV header");
            }
            debugf("Got a metadata line\n");
            got_meta = true;
            meta_count++;
        } else {
            if (!got_meta) {
                throw VcfValidationError("Got a header line but no metadata lines");
            }
            got_header = true;
            header_count++;
        }

        // this is a metadata or header line, so read the rest of the line
        size_t tab_count = 0;
        linebuf.push_back(c1);
        linebuf.push_back(c2);
        // string_appendc(&linebuf, c1);
        // string_appendc(&linebuf, c2);

        while (true) {
            char c3;
            if (read(input_fd, &c3, 1) <= 0) {
                throw VcfValidationError("Failed to read the rest of the metadata or header row!");
            } else if (c3 == '\n') {
                linebuf.push_back(c3);
                // string_appendc(&linebuf, '\n');
                break;
            } else if (got_header && c3 == '\t') { // got_header is only true once we are parsing a TSV header line
                tab_count++;
                if (tab_count > VCF_REQUIRED_COL_COUNT) {
                    output_schema.sample_count++;
                }
            }
            linebuf.push_back(c3);
            // string_appendc(&linebuf, c3);
        }

        debugf("Line: %s\n", linebuf.c_str());
        output_vector.push_back(linebuf);
    }
    debugf("Line counts: metadata = %ld, header = %ld\n", meta_count, header_count);
    debugf("Sample count: %ld\n", output_schema.sample_count);
    #ifdef TIMING
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("decompress2_metadata_headers_fd time: %lu\n", duration.count());
    #endif
    return 0;
}


int decompress2_fd(const std::string& input_filename, const std::string& output_filename) {
    debugf("Decompressing %s to %s\n", input_filename.c_str(), output_filename.c_str());
    int input_fd = open(input_filename.c_str(), O_RDONLY);
    int output_fd = open(output_filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY, DEFAULT_FILE_CREATE_MODE);
    VcfCompressionSchema schema;

    std::vector<std::string> meta_header_lines;
    meta_header_lines.reserve(256);
    decompress2_metadata_headers_fd(input_fd, meta_header_lines, schema);
    for (size_t i = 0; i < meta_header_lines.size(); i++) {
        // these lines still have the newline char included
        std::string& line = meta_header_lines.at(i);
        write(output_fd, line.c_str(), line.size());
    }
    size_t variant_line_count = 0;

    std::string variant_line;
    variant_line.reserve(16 * 1024); // 16 KiB
    // string_t variant_line;
    // string_init(&variant_line);

    while (true) {
        variant_line_count++;
        // string_clear(&variant_line);
        variant_line.clear();

        size_t compressed_line_length = 0;
        int status = decompress2_data_line_fd2_string(input_fd, schema, variant_line, &compressed_line_length);
        if (status == 0) {
            debugf("Finished reading file\n");
            break;
        } else if (status < 0) {
            throw VcfValidationError("Failed to decompress file");
        } else {
            write(output_fd, variant_line.c_str(), variant_line.size());
        }

    } // end line loop
    debugf("variant_line_count: %lu\n", variant_line_count);
    close(input_fd);
    close(output_fd);

    return 0;
}
