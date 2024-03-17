#include "analysis.h"

#include <algorithm>
#include <future>
#include <thread>
#include <fstream>
#include <utility>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include <chrono>
#include <map>
#include <filesystem>

namespace fs = std::filesystem; 

double correlation(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    double m_a {}, d_a {}, m_b {}, d_b {};
    int c_size {static_cast<int>(std::min(a.size(), b.size()))};

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

    return r;
}

void rgb_correlation(const BMP& file) {

    std::vector<uint8_t> r(file.get_size_image()), g(file.get_size_image()), b(file.get_size_image());
    const std::vector<uint8_t> ddata {file.get_data()};
    for (size_t i {}, j {}; i <= ddata.size(); i += 3, ++j) {
        r[j] = ddata[i + 2];
        g[j] = ddata[i + 1];
        b[j] = ddata[i    ];
    }
    std::cout << "Correlation of image components: [R,G] = " << correlation(r, g) << std::endl;
    std::cout << "Correlation of image components: [G,B] = " << correlation(g, b) << std::endl;
    std::cout << "Correlation of image components: [B,R] = " << correlation(b, r) << std::endl;

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

std::map<int, double> auto_correlation(std::vector<uint8_t> data, int32_t w, int32_t h, int y) {
    std::map<int, double> res;

    /*
    for (int x {-(w / 4 - 1)}; x < w / 4; ++x) {
        std::vector<uint8_t> a;
        std::vector<uint8_t> b;
        for (int i {1}; i < h - y; ++i) {
            for (int j {1}; j < w - x; ++j) {
                a.push_back(data[i * w + j]);
            }
        }
        for (int m {y + 1}; m < h; ++m) {
            for (int n {x + 1}; n < w; ++n) {
                b.push_back(data[m * w + n]);
            }
        }
        res.insert(std::make_pair(x, correlation(a, b)));
    }
    */
    int upper_bound {h}, lower_bound {},
        upper_bound_slide {h}, lower_bound_slide {};
    if (y < 0) {
        lower_bound_slide -= y;
        upper_bound += y; 
    } else {
        upper_bound_slide -= y;
        lower_bound += y;
    }
    for (int x {-(w / 4 - 1)}; x < 0; ++x) {
        std::vector<uint8_t> a;
        std::vector<uint8_t> b;
        // static image
        for (int i {lower_bound}; i < upper_bound; ++i) {
            for (int j {0}; j < w + x; ++j) {
                a.push_back(data[i * w + j]);
            }
        }
        // sliding image
        for (int m {lower_bound_slide}; m < upper_bound_slide; ++m) {
            for (int n {+x}; n < w; ++n) {
                b.push_back(data[m * w + n]);
            }
        }
        res.insert(std::make_pair(x, correlation(a, b)));
    }
    for (int x {0}; x < w / 4; ++x) {
        std::vector<uint8_t> a;
        std::vector<uint8_t> b;
        // static image
        for (int i {lower_bound}; i < upper_bound; ++i) {
            for (int j {x}; j < w; ++j) {
                a.push_back(data[i * w + j]);
            }
        }
        // sliding image
        for (int m {lower_bound_slide}; m < upper_bound_slide; ++m) {
            for (int n {0}; n < w - x; ++n) {
                b.push_back(data[m * w + n]);
            }
        }
        res.insert(std::make_pair(x, correlation(a, b)));
    }
    return res;
}

void auto_correlation_by_channel_along_y(std::vector<uint8_t> &data, int32_t w, int32_t h, int start, int end, int step, const std::string &filename) {
    fs::path fpath(fs::current_path().parent_path().parent_path() / "tmp");
    if (!fs::is_directory(fpath)) {
        fs::create_directory(fpath);
    }
    std::ofstream fout(fpath / filename);

    int thread_count {(end - start) / step + 1};
    std::vector<std::map<int, double>> au_cr_results(thread_count);
    std::vector<std::future<std::map<int, double>>> au_cr_futures(thread_count);

    int i {};
    for (int y {start}; y <= end; y += step) {
        au_cr_futures[i++] = std::async(std::launch::async, auto_correlation, data, w, h, y);
    }
    for (i = 0; i < thread_count; ++i) {
        au_cr_results[i] = au_cr_futures[i].get();
    }
    for (int k {-(w / 4 - 1)}; k < w / 4; ++k) {
        fout << k; 
        for (int j {}; j < thread_count; ++j) {
            fout << " " << au_cr_results[j].at(k);
        }
        fout << "\n";
    }
}

void rgb_auto_correlation(const BMP& file, const std::string &filename) {
    int32_t im_size {file.get_height() * file.get_width()},
            w {file.get_width()},
            h {file.get_height()};
    std::vector<uint8_t> channel(im_size);
    const std::vector<uint8_t> &data = file.get_data();

    // r channel auto correlation
    for (int i {}, j {}; i < data.size(); i += 3, ++j) {
        channel[j] = data[i + 2];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_r.txt"));
    // g channel auto correlation
    for (int i {}, j {}; i < data.size(); i += 3, ++j) {
        channel[j] = data[i + 1];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_g.txt"));
    // b channel auto correlation
    for (int i {}, j {}; i < data.size(); i += 3, ++j) {
        channel[j] = data[i];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_b.txt"));
}
