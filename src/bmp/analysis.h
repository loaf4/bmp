#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bmp.h"

#include <map>
#include <vector>
#include <cstdint>
#include <cmath>

// correlation
void calc_correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b);
void rgb_correlation(const BMP& file);
double math_expectation(const std::vector<uint8_t>& data, const uint32_t width, const uint32_t height, const char component);
double standard_deviation(const std::vector<uint8_t>& data, double m_e, uint32_t width, uint32_t height, const char component);
// auto correlation
std::map<int, double> auto_correlation(std::vector<uint8_t> data, int32_t w, int32_t h, int y);
void auto_correlation_by_channel_along_y(std::vector<uint8_t> &data, int32_t w, int32_t h, int start, int end, int step, const std::string &filename);
void rgb_auto_correlation(const BMP& file, const std::string &filename);

#endif
