#!/usr/bin/env bash
set -e
cd cpp/release
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../..

cp cpp/release/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so