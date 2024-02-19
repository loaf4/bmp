#include "bmp.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    BMP f("goldhill");
    f.save_file_by_component("goldhill", 'r');
    f.save_file_by_component("goldhill", 'g');
    f.save_file_by_component("goldhill", 'b');

    return 0;
}

