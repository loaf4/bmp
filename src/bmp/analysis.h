#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bmp.h"

#include <map>
#include <vector>
#include <cstdint>
#include <cmath>

// correlation
double correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b);
void rgb_correlation(const BMP& file);
void ycbcr_correlation(const BMP& file);
double math_expectation(const std::vector<uint8_t>& data, const uint32_t width, const uint32_t height, const char component);
double standard_deviation(const std::vector<uint8_t>& data, double m_e, uint32_t width, uint32_t height, const char component);
// auto correlation
std::map<int, double> auto_correlation(std::vector<uint8_t> data, int32_t w, int32_t h, int y);
void auto_correlation_by_channel_along_y(const std::vector<uint8_t> &data, int32_t w, int32_t h, int start, int end, int step, const std::string &filename);
void rgb_auto_correlation(const BMP& file, const std::string &filename);
void ycbcr_auto_correlation(const BMP& file, const std::string &filename);
void save_rgb_to_ycbcr(const std::string &fname, const BMP &file);
uint8_t saturation(double x, int x_min, int x_max);
 

#endif

