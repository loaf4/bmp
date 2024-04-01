#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bmp.h"

#include <map>
#include <vector>
#include <cstdint>
#include <cmath>

void rgb_component_frequency(const std::vector<uint8_t> &data);
void ycbcr_component_frequency(const std::vector<uint8_t> &data);

// correlation
double correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b);
void rgb_correlation(const BMP& file);
void ycbcr_correlation(const std::vector<uint8_t> &data);

// auto correlation
std::map<int, double> auto_correlation(const std::vector<uint8_t> &data, int32_t w, int32_t h, int y);
void auto_correlation_by_channel_along_y(const std::vector<uint8_t> &data, int32_t w, int32_t h, int start, int end, int step, const std::string &filename);
void rgb_auto_correlation(const BMP& file, const std::string &filename);
void ycbcr_auto_correlation(const BMP& file, const std::string &filename, const std::vector<uint8_t> &data);

// ycbcr
std::vector<uint8_t> rgb_to_ycbcr(const BMP &file);
std::vector<uint8_t> ycbcr_to_rgb(const std::vector<uint8_t> &data);
void save_rgb_to_ycbcr(const std::string &fname, const BMP &file, const std::vector<uint8_t> &data);
uint8_t saturation(double x, int x_min, int x_max);
double psnr(const std::vector<uint8_t> &a, const std::vector<uint8_t> &b);
void rgb_ycbcr_psnr(const std::vector<uint8_t> &rgb_data, const std::vector<uint8_t> &ycbcr_data, int w, int h);
void decimation_even(const BMP &file, int factor);
void decimation_square(const BMP &file, int factor);

#endif

