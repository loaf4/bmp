#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bmp.h"

#include <vector>
#include <cstdint>
#include <cmath>

// correlation
void calc_correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b);
void rgb_correlation(const BMP& file);
double math_expectation(const std::vector<uint8_t>& data, const uint32_t width, const uint32_t height, const char component);
double standard_deviation(const std::vector<uint8_t>& data, double m_e, uint32_t width, uint32_t height, const char component);
void auto_correlation(const BMP& file, int y);

#endif

