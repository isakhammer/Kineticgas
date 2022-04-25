#!/usr/bin/env bash
set -e
cd cpp/release
if [ $1 -e "-cleancache" ] || [ $2 -e "-cleancache"]; then
    echo $1
    echo "Cleaning!"
    rm -rf CMakeCache.txt
fi
if [ $1 -e "-fullclean" ] || [ $2 -e "-fullclean"]; then
    echo $1
    echo "Wiping!"
    cd ..
    rm -rf release
    mkdir release
    cd release
fi
if [ $1 -e "-Debug" ]; then
    cmake -DCMAKE_BUILD_TYPE=Debug
else
    cmake -DCMAKE_BUILD_TYPE=Release  ..
fi

make
cd ../..

cp cpp/release/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so
python -m pykingas -test