#!/usr/bin/env bash
set -e
cd cpp/release
rm -rf CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../..

cp cpp/release/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so
python pykingas/test_inside.py