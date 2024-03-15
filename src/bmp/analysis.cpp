#include "analysis.h"

#include <cmath>
#include <cstdint>
#include <vector>
#include <chrono>

void calc_correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    double m_a {}, d_a {}, m_b {}, d_b {};
    int c_size {static_cast<int>(a.size())};

    // math expectation
    for (int i {}; i < c_size; ++i) {
        m_a += a[i];
        m_b += b[i];
    }
    m_a /= c_size;
    m_b /= c_size;

    // standard deviation
    for (int i {}; i < c_size; ++i) {
        d_a += std::pow(a[i] - m_a, 2);
        d_b += std::pow(b[i] - m_b, 2);
    }
    d_a = std::sqrt(d_a / (c_size - 1));
    d_b = std::sqrt(d_b / (c_size - 1));

    // correlation
    double r {};
    for (int i {}; i < c_size; ++i) {
        r += (a[i] - m_a) * (b[i] - m_b);
    }
    r /= c_size;
    r /= d_a * d_b;

    std::cout << "Correlation of components = " << r << std::endl;
}

void rgb_correlation(const BMP& file) {

    std::vector<uint8_t> r(file.get_size_image()), g(file.get_size_image()), b(file.get_size_image());
    const std::vector<uint8_t> ddata {file.get_data()};
    for (size_t i {}, j {}; i <= ddata.size(); i += 3, ++j) {
        r[j] = ddata[i + 2];
        g[j] = ddata[i + 1];
        b[j] = ddata[i    ];
    }
    std::cout << "Correlation of image components: [R,G]" << std::endl;
    calc_correlation(r, g);
    std::cout << "Correlation of image components: [G,B]" << std::endl;
    calc_correlation(g, b);
    std::cout << "Correlation of image components: [B,R]" << std::endl;
    calc_correlation(b, r);

    // previous version
    /*
    double m_r {}, d_r {}, m_g {}, d_g {}, m_b {}, d_b {};
    const std::vector<uint8_t> data {file.get_data()}; 
    int img_size {static_cast<int>(file.get_size_image())};

    // math expectation
    for (size_t i {}; i < img_size; i += 3) {
        m_r += data[i + 2];
        m_g += data[i + 1];
        m_b += data[i    ];
    }
    m_r /= (img_size);
    m_g /= (img_size);
    m_b /= (img_size);

    // standard deviation
    for (size_t i {}; i < img_size; i += 3) {
        d_r += std::pow(data[i + 2] - m_r, 2);
        d_g += std::pow(data[i + 1] - m_g, 2);
        d_b += std::pow(data[i    ] - m_b, 2);
    }
    d_r = std::sqrt(d_r / (img_size - 1));
    d_g = std::sqrt(d_g / (img_size - 1));
    d_b = std::sqrt(d_b / (img_size - 1));

    // correlation
    double r_rg {}, r_gb {}, r_br {};
    for (size_t i {}; i < img_size; i += 3) {
        r_rg += (data[i + 2] - m_r) * (data[i + 1] - m_g);
        r_gb += (data[i + 1] - m_g) * (data[i    ] - m_b);
        r_br += (data[i    ] - m_b) * (data[i + 2] - m_r);
    }
    r_rg /= img_size;
    r_gb /= img_size;
    r_br /= img_size;
    r_rg /= d_r * d_g;
    r_gb /= d_g * d_b;
    r_br /= d_b * d_r;

    std::cout << "Correlation of imgae components: [R,G] = " << r_rg << ", [G,B] = " << r_gb << ", [B,R] = " << r_br << std::endl;
    */
}

double math_expectation(const std::vector<uint8_t>& data, const uint32_t width, const uint32_t height, const char component) {
    double m {};
    switch (component) {
        case 'r':
            for (size_t i {}; i < data.size(); i += 3) {
                m += data[i];
            }
            break;
        case 'g':
            for (size_t i {1}; i < data.size(); i += 3) {
                m += data[i];
            }
            break;
        case 'b':
            for (size_t i {2}; i < data.size(); i += 3) {
                m += data[i];
            }
            break;
        default:
            break;
    }
    return m / (width * height);
}

double standard_deviation(const std::vector<uint8_t>& data, double m_e, uint32_t width, uint32_t height, const char component) {
    double d {};
    switch (component) {
        case 'r':
            for (size_t i {}; i < data.size(); i += 3) {
                d += std::pow(m_e - data[i], 2);
            }
            break;
        case 'g':
            for (size_t i {1}; i < data.size(); i += 3) {
                d += std::pow(m_e - data[i], 2);
            }
            break;
        case 'b':
            for (size_t i {2}; i < data.size(); i += 3) {
                d += std::pow(m_e - data[i], 2);
            }
            break;
        default:
            break;
    }
    return std::sqrt(d / (width * height - 1));
}

void auto_correlation(const BMP &file, int y) {
    for (int x {-(file.get_width() / 4 - 1)}; x < file.get_width() / 4; ++x) {

    }
}
