#include "analysis.h"
#include "bmp.h"

#include <algorithm>
#include <future>
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
    const std::vector<uint8_t> &data {file.get_data()};
    for (size_t i {}, j {}; i <= data.size(); i += 3, ++j) {
        r[j] = data[i + 2];
        g[j] = data[i + 1];
        b[j] = data[i    ];
    }
    std::cout << "Correlation of image components: [R,G] = " << correlation(r, g) << std::endl;
    std::cout << "Correlation of image components: [G,B] = " << correlation(g, b) << std::endl;
    std::cout << "Correlation of image components: [B,R] = " << correlation(b, r) << std::endl;
}

void ycbcr_correlation(const std::vector<uint8_t> &data) {
    std::vector<uint8_t> y_data(data.size()), cb_data(data.size()), cr_data(data.size());

    for (size_t i {}, j {}; j != data.size() && i != data.size(); i += 3, ++j) {
        y_data[j] = data[i];
        cb_data[j] = data[i + 1];
        cr_data[j] = data[i + 2];
    }

    std::cout << "Correlation of image components: [Y,Cb] = " << correlation(y_data, cb_data) << std::endl;
    std::cout << "Correlation of image components: [Cb,Cr] = " << correlation(cb_data, cr_data) << std::endl;
    std::cout << "Correlation of image components: [Cr,Y] = " << correlation(cr_data, y_data) << std::endl;
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

std::map<int, double> auto_correlation(const std::vector<uint8_t> &data, int32_t w, int32_t h, int y) {
    std::map<int, double> res;

    int upper_bound {h}, lower_bound {},
        upper_bound_slide {h}, lower_bound_slide {};
    if (y < 0) {
        lower_bound_slide -= y;
        upper_bound += y; 
    } else {
        upper_bound_slide -= y;
        lower_bound += y;
    }

    // -(w / 4 - 1) < x < 0
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
            for (int n {std::abs(x)}; n < w; ++n) {
                b.push_back(data[m * w + n]);
            }
        }
        res.insert(std::make_pair(x, correlation(a, b)));
    }
    // 0 < x < (w / 4)
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

void auto_correlation_by_channel_along_y(const std::vector<uint8_t> &data, int32_t w, int32_t h, int start, int end, int step, const std::string &filename) {
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
    int32_t size {file.get_height() * file.get_width()},
            w {file.get_width()},
            h {file.get_height()};
    std::vector<uint8_t> channel(size);
    const std::vector<uint8_t> &data = file.get_data();

    // r channel auto correlation
    for (size_t i {}, j {}; i != data.size(); i += 3, ++j) {
        channel[j] = data[i + 2];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_r.txt"));
    // g channel auto correlation
    for (size_t i {}, j {}; i != data.size(); i += 3, ++j) {
        channel[j] = data[i + 1];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_g.txt"));
    // b channel auto correlation
    for (size_t i {}, j {}; i != data.size(); i += 3, ++j) {
        channel[j] = data[i];
    }
    auto_correlation_by_channel_along_y(channel, w, h, -10, 10, 5, (filename + "_b.txt"));
}

void ycbcr_auto_correlation(const BMP& file, const std::string &filename, const std::vector<uint8_t> &data) {
    int32_t size {file.get_height() * file.get_width()},
            w {file.get_width()},
            h {file.get_height()},
            j {};
    std::vector<uint8_t> y_data(size), cb_data(size), cr_data(size);

    for (size_t i {}; j < size && i != data.size(); i += 3, ++j) {
        y_data[j] = data[i];
        cb_data[j] = data[i + 1];
        cr_data[j] = data[i + 2];
    }

    // y channel auto correlation
    auto_correlation_by_channel_along_y(y_data, w, h, -10, 10, 5, (filename + "_y.txt"));
    // cb channel auto correlation
    auto_correlation_by_channel_along_y(cb_data, w, h, -10, 10, 5, (filename + "_cb.txt"));
    // cr channel auto correlation
    auto_correlation_by_channel_along_y(cr_data, w, h, -10, 10, 5, (filename + "_cr.txt"));
}

std::vector<uint8_t> rgb_to_ycbcr(const BMP &file) {
    const std::vector<uint8_t> &data {file.get_data()};
    std::vector<uint8_t> ycbcr_data(data.size());

    for (size_t i {}; i != data.size(); i += 3) {
        uint8_t y = saturation(0.299 * data[i + 2] + 0.587 * data[i + 1] + 0.114 * data[i], 0, 255);
        uint8_t cb = saturation(0.5643 * (data[i] - y) + 128, 0, 255);
        uint8_t cr = saturation(0.7132 * (data[i + 2] - y) + 128, 0, 255);

        ycbcr_data[i    ] = y;
        ycbcr_data[i + 1] = cb;
        ycbcr_data[i + 2] = cr;
    }

    return ycbcr_data;
}

std::vector<uint8_t> ycbcr_to_rgb(const std::vector<uint8_t> &data) {
    std::vector<uint8_t> rgb_data(data.size());

    for (size_t i {}; i != data.size(); i += 3) {
        uint8_t g = saturation(data[i] - 0.714 * (data[i + 2] - 128) - 0.334 * (data[i + 1] - 128), 0, 255);
        uint8_t r = saturation(data[i] + 1.402 * (data[i + 2] - 128), 0, 255);
        uint8_t b = saturation(data[i] + 1.772 * (data[i + 1] - 128), 0, 255);

        rgb_data[i    ] = b;
        rgb_data[i + 1] = g;
        rgb_data[i + 2] = r;
    }

    return rgb_data;
}

void save_rgb_to_ycbcr(const std::string &fname, const BMP &file, const std::vector<uint8_t> &data) {
    std::vector<uint8_t> y_data (file.get_size_image());
    std::vector<uint8_t> cb_data (file.get_size_image());
    std::vector<uint8_t> cr_data (file.get_size_image());

    for (size_t i {}; i != data.size(); i += 3) {

        y_data[i] = y_data[i + 1] = y_data[i + 2] = data[i];
        cb_data[i] = cb_data[i + 1] = cb_data[i + 2] = data[i + 1];
        cr_data[i] = cr_data[i + 1] = cr_data[i + 2] = data[i + 2];
    }

    file.save_file(fname + "_by_y_component.bmp", y_data);
    file.save_file(fname + "_by_cb_component.bmp", cb_data);
    file.save_file(fname + "_by_cr_component.bmp", cr_data);
}

uint8_t saturation(double x, int x_min, int x_max) {
    if (x < x_min) {
        return x_min;
    } else if (x > x_max) {
        return x_max;
    }
    return static_cast<uint8_t>(x);
} 

double psnr(const std::vector<uint8_t> &a, const std::vector<uint8_t> &b) {
    double res {};
    int size {static_cast<int>(std::min(a.size(), b.size()))};

    for (int i {}; i < size; ++i) {
        res += std::pow(a[i] - b[i], 2);
    }
    res = 10 * std::log10(size * std::pow((std::pow(2, 8) - 1), 2) / res);

    return res;
}

void rgb_ycbcr_psnr(const std::vector<uint8_t> &rgb_data, const std::vector<uint8_t> &ycbcr_data, int w, int h) {
    std::vector<uint8_t> recovered_rgb_data {ycbcr_to_rgb(ycbcr_data)};

    std::vector<uint8_t> a(w * h);
    std::vector<uint8_t> b(w * h);

    // PSNR of r
    for (size_t i {}; i < rgb_data.size(); i += 3) {
        a[i] = rgb_data[i + 2];
        b[i] = recovered_rgb_data[i + 2];
    }
    std::cout << "RGB -> YCbCr -> RGB; R PSNR = " << psnr(a, b) << std::endl;
    // PSNR of g
    for (size_t i {}; i < rgb_data.size(); i += 3) {
        a[i] = rgb_data[i + 1];
        b[i] = recovered_rgb_data[i + 1];
    }
    std::cout << "RGB -> YCbCr -> RGB; G PSNR = " << psnr(a, b) << std::endl;
    // PSNR of b
    for (size_t i {}; i < rgb_data.size(); i += 3) {
        a[i] = rgb_data[i];
        b[i] = recovered_rgb_data[i];
    }
    std::cout << "RGB -> YCbCr -> RGB; B PSNR = " << psnr(a, b) << std::endl;
}
