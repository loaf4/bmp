#include "bmp.h"
#include "analysis.h"

#include <iostream>
#include <vector>
#include <cstdint>

int main(int argc, char* argv[]) {
    BMP f("goldhill.bmp");
    const std::vector<uint8_t> &rgb_data {f.get_data()};
    std::vector<uint8_t> ycbcr_data {rgb_to_ycbcr(f)};

    /*
    // output file header
    f.print_header();

    // get .bmp files by r, g, b components ( 3 )
    f.save_file_by_component("goldhill", 'r');
    f.save_file_by_component("goldhill", 'g');
    f.save_file_by_component("goldhill", 'b');

    // output correlation of r, g, b components ( 4 )
    rgb_correlation(f);

    // create files with auto correlation with y=[-10, -5, 0, 5, 10] ( 4 )
    rgb_auto_correlation(f, "goldhill_au_cr");

    // get .bmp files by y, cb, cr components ( 6 )
    save_rgb_to_ycbcr("goldhill", f, ycbcr_data);

    // output correlation of y, cb, cr components ( 5 )
    ycbcr_correlation(ycbcr_data);

    // create files with auto correlation with y=[-10, -5, 0, 5, 10] ( 5 )
    ycbcr_auto_correlation(f, "goldhill_au_cr", ycbcr_data);
    
    // save recovered YCbCr to RGB ( 7 )
    f.save_file("goldhill_after_ycbcr", rgb_data);

    // calculate PSNR of original and recovered RGB data
    rgb_ycbcr_psnr(f.get_data(), ycbcr_data, f.get_width(), f.get_height());

    // decimation ( 8 - 11 )
    decimation_even(f, 2);
    decimation_even(f, 4);
    decimation_square(f, 2);
    decimation_square(f, 4);

    // component frequency ( 12 )
    rgb_component_frequency(rgb_data);
    ycbcr_component_frequency(ycbcr_data);

    // entropy ( 13 )
    rgb_entropy(rgb_data);
    ycbcr_entropy(ycbcr_data);

    // DPCM ( 14 - 16 )
    rgb_DPCM(rgb_data, f.get_width(), f.get_height());
    ycbcr_DPCM(ycbcr_data, f.get_width(), f.get_height());
    */

    return 0;
}

