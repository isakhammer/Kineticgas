#!/usr/bin/env bash
set -e
cd cpp/release
if [ $1 -e "-cleancache" ]; then
    echo $1
    echo "Cleaning!"
    rm -rf CMakeCache.txt
fi
if [ $1 -e "-fullclean" ]; then
    echo $1
    echo "Wiping!"
    cd ..
    rm -rf release
    mkdir release
    cd release
fi
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../..

cp cpp/release/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so
python -m pykingas -test