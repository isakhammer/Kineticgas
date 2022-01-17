#!/usr/bin/env bash
set -e
cd cpp/release
cmake -O3 -Wall -Wextra -Release ..
make
cd ../..