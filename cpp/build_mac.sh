#!/usr/bin/env bash
set -e
cd cpp/release
cmake -O3 -Wall -Wextra -Release ..
make
cd ../..

cp cpp/release/KineticGas.cpython-39-darwin.so python/KineticGas.cpython-39-darwin.so