#include "bmp.h"
#include "analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    BMP f("goldhill.bmp");

    f.print_header();
    rgb_correlation(f);

    return 0;
}

