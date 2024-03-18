#include "bmp.h"
#include "analysis.h"
#include "bmp/analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    BMP f("goldhill.bmp");

    // output file header
    f.print_header();

    // get .bmp files by r, g, b components
    f.save_file_by_component("goldhill", 'r');
    f.save_file_by_component("goldhill", 'g');
    f.save_file_by_component("goldhill", 'b');

    /*
    // output correlation of r, g, b components
    rgb_correlation(f);

    // create files with auto correlation with y=[-10, -5, 0, 5, 10]
    rgb_auto_correlation(f, "goldhill_au_cr");
    */

    // get .bmp files by y, cb, cr components
    save_rgb_to_ycbcr("goldhill", f);

    return 0;
}

