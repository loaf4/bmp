#include "bmp.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

namespace fs = std::filesystem;

BMP::BMP(const std::string& fname) {
    fs::path fpath(fs::current_path().parent_path().parent_path() / "img" / (fname + ".bmp"));
    std::ifstream fin(fpath, std::ios::binary);
    if (fin) {
        fin.read(reinterpret_cast<char*>(&_file_header), sizeof(_file_header));
        if (_file_header.bf_type != 0x4D42) {
            throw std::runtime_error("error: unknown file format");
        }
        fin.read(reinterpret_cast<char*>(&_info_header), sizeof(_info_header));

        _data.resize(_info_header.bi_width * _info_header.bi_height * _info_header.bi_bit_count / 8);
        fin.read(reinterpret_cast<char*>(_data.data()), _data.size());
    }
}

void BMP::save_file(const std::string& fname) {
    save_file(fname, _data);
}

void BMP::save_file(const std::string& fname, const std::vector<uint8_t>& data) {
    fs::path fpath(fs::current_path().parent_path().parent_path() / "img");
    if (!fs::is_directory(fpath)) {
        fs::create_directory(fpath);
    }
    std::ofstream fout(fpath / (fname + ".bmp"), std::ios::binary);
    if (fout) {
        fout.write(reinterpret_cast<char*>(&_file_header), sizeof(_file_header));
        fout.write(reinterpret_cast<char*>(&_info_header), sizeof(_info_header));
        fout.write(reinterpret_cast<char*>(const_cast<std::vector<uint8_t>&>(data).data()), data.size());
    }
}

void BMP::save_file_by_component(const std::string& fname, const char mod) {
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
    save_file(fname + "_by_" + mod + "_component", tmp_data);
}
