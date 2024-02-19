#include "bmp.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    BMP f("goldhill");
    f.save_file_by_component("111", 'r');

    return 0;
}

