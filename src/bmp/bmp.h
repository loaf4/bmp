#ifndef BMP_H
#define BMP_H

#include <iostream>
#include <cstdint>
#include <vector>
#include <string>

class BMP {
#pragma pack(push)
#pragma pack(1)
    struct BMPFileHeader {
        uint16_t bf_type; // must be 0x4D42
        uint32_t bf_size;
        uint16_t bf_reserved1;
        uint16_t bf_reserved2;
        uint32_t bf_off_bits;
    } _file_header;
    struct BMPInfoHeader {
        uint32_t bi_size;
        int32_t bi_width;
        int32_t bi_height;
        uint16_t bi_planes;
        uint16_t bi_bit_count;
        uint32_t bi_compression;
        uint32_t bi_size_image;
        int32_t bi_x_pixels_per_meter;
        int32_t bi_y_pixels_per_meter;
        uint32_t bi_colors_used; 
        uint32_t bi_colors_important;
    } _info_header;
#pragma pack(pop)
    std::vector<uint8_t> _data;


public:
    BMP() = default;
    BMP(const std::string& fname);

    void save_file(const std::string& fname);
    void save_file(const std::string& fname, const std::vector<uint8_t>& data);
    void save_file_by_component(const std::string& fname, const char mod);

    void print_header() {
        std::cout << 
            _file_header.bf_type << " " << 
            _file_header.bf_size << " " << 
            _file_header.bf_reserved1 << " " <<
            _file_header.bf_reserved2 << " " <<
            _file_header.bf_off_bits << " " << std::endl;

        std::cout <<
            _info_header.bi_size << " " << 
            _info_header.bi_width << " " <<
            _info_header.bi_height << " " <<
            _info_header.bi_planes << " " <<
            _info_header.bi_bit_count << " " <<
            _info_header.bi_compression << " " <<
            _info_header.bi_size_image << " " <<
            _info_header.bi_x_pixels_per_meter << " " <<
            _info_header.bi_y_pixels_per_meter << " " <<
            _info_header.bi_colors_used << " " << 
            _info_header.bi_colors_important << " " << std::endl;
    }
};

#endif

