#include "bmp.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    BMP f("parrot.bmp");
    f.save_file("erere.bmp");

    f.print_header();

    return 0;
}

