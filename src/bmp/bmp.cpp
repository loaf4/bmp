#include "bmp.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

namespace fs = std::filesystem;

BMP::BMP(const std::string& fname) {
    fs::path fpath(fs::current_path().parent_path().parent_path() / "img" / fname);
    std::ifstream fin(fpath, std::ios::binary);
    if (fin) {
        // read header of BMP file
        fin.read(reinterpret_cast<char*>(&_file_header), sizeof(_file_header));
        if (_file_header.bf_type != 0x4D42) {
            throw std::runtime_error("error: unknown file format");
        }
        fin.read(reinterpret_cast<char*>(&_info_header), sizeof(_info_header));
        // get type of BMP file
        _bmp_type = static_cast<BMPType>(_info_header.bi_bit_count);
        // read the palette
        if (_info_header.bi_colors_used != 0) {
            _palette.resize(_info_header.bi_colors_used);
            fin.read(reinterpret_cast<char*>(_palette.data()), _info_header.bi_colors_used * sizeof(RGBQuad));
        }

        // fseek to start of the pixels
        fin.seekg(0);
        fin.seekg(_file_header.bf_off_bits);

        // read pixels
        _data.resize(_info_header.bi_width * _info_header.bi_height * _info_header.bi_bit_count / 8);
        int multiplier {};
        switch (_info_header.bi_bit_count) {
            case 8:
                multiplier = 3;
                break;
            case 16:
                multiplier = 2;
                break;
            case 24:
                multiplier = 1;
                break;
        }
        int step {_info_header.bi_width * _info_header.bi_bit_count / 8};
        if (step % 4 == 0) {
            fin.read(reinterpret_cast<char*>(_data.data()), _data.size());
        } else {
            size_t padding {static_cast<size_t>(multiplier * _info_header.bi_width % 4)};
            std::vector<uint8_t> garbage(padding);

            for (int y {}; y < _info_header.bi_height; ++y) {
                fin.read(reinterpret_cast<char*>(_data.data()) + y * step, step);
                fin.read(reinterpret_cast<char*>(garbage.data()), garbage.size());
            }
        }
    }
}

void BMP::save_file(const std::string& fname) const {
    save_file(fname, _data, _palette);
}

void BMP::save_file(const std::string& fname, const std::vector<uint8_t>& data, const std::vector<RGBQuad>& palette) const {
    fs::path fpath(fs::current_path().parent_path().parent_path() / "img");
    if (!fs::is_directory(fpath)) {
        fs::create_directory(fpath);
    }
    std::ofstream fout(fpath / fname, std::ios::binary);
    if (fout) {
        // write header
        fout.write(reinterpret_cast<char*>(const_cast<struct BMPFileHeader*>(&_file_header)), sizeof(_file_header));
        fout.write(reinterpret_cast<char*>(const_cast<struct BMPInfoHeader*>(&_info_header)), sizeof(_info_header));

        // write palette and offset garbage
        if (_file_header.bf_off_bits > 54) {
            int offset {static_cast<int>(_file_header.bf_off_bits - sizeof(_file_header) - sizeof(_info_header))};
            if (_info_header.bi_colors_used != 0 && !_palette.empty()) {
                offset -= _palette.size();
                fout.write(reinterpret_cast<char*>(const_cast<RGBQuad*>(palette.data())), _palette.size());
            }
            std::vector<uint8_t> garbage(offset);
            fout.write(reinterpret_cast<char*>(garbage.data()), garbage.size());
        }
        
        // write data
        int multiplier {};
        switch (_info_header.bi_bit_count) {
            case 8:
                multiplier = 3;
                break;
            case 16:
                multiplier = 2;
                break;
            case 24:
                multiplier = 1;
                break;
        }
        int step {_info_header.bi_width * _info_header.bi_bit_count / 8};
        if (step % 4 == 0) {
            fout.write(reinterpret_cast<char*>(const_cast<uint8_t*>(data.data())), data.size());
        } else {
            size_t padding {static_cast<size_t>(multiplier * _info_header.bi_width % 4)};
            std::vector<uint8_t> garbage(padding);

            for (int y {}; y < _info_header.bi_height; ++y) {
                fout.write(reinterpret_cast<char*>(const_cast<uint8_t*>(data.data())) + y * step, step);
                fout.write(reinterpret_cast<char*>(garbage.data()), garbage.size());
            }
        }
    }
}

void BMP::save_file_by_component(const std::string& fname, const char mod) const {
    std::vector<uint8_t> tmp_data {_data};
    switch (mod) {
        case 'r':
            for (size_t i {}; i < tmp_data.size(); i += 3) {
                tmp_data[i    ] &= 0x00;
                tmp_data[i + 1] &= 0x00;
            }
            break;
        case 'g':
            for (size_t i {}; i < tmp_data.size(); i += 3) {
                tmp_data[i    ] &= 0x00;
                tmp_data[i + 2] &= 0x00;
            }
            break;
        case 'b':
            for (size_t i {}; i < tmp_data.size(); i += 3) {
                tmp_data[i + 1] &= 0x00;
                tmp_data[i + 2] &= 0x00;
            }
            break;
        default:
            break;
    }
    save_file(fname + "_by_" + mod + "_component.bmp", tmp_data);
}
