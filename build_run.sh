#!/bin/bash
set -e
./build.sh
cd build/src
./BMP
cd ../../
