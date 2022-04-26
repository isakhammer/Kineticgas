#!/usr/bin/env bash
set -e
cd cpp/release

while test $# -gt 0
do
    case "$1" in
        --fullclean)
            echo "Wiping release directory"
            cd ..
            rm -rf release
            mkdir release
            cd release
            ;;
        --cleancache)
            echo "Cleaning cache"
            rm -rf CMakeCache.txt
            ;;
        --Debug)
            echo "Building Debug"
            cd ../debug
            cmake -DCMAKE_BUILD_TYPE=Debug ..
            make
            cd ../..
            cp cpp/debug/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so
            cp cpp/debug/KineticGas.cpython-39-darwin.so cpp/debug/KineticGas.so
            python -m pykingas -test -debug
            exit 0
            ;;
    esac
    shift
done

echo "Building Release"
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../..
cp cpp/release/KineticGas.cpython-39-darwin.so pykingas/KineticGas.so
cp cpp/release/KineticGas.cpython-39-darwin.so cpp/release/KineticGas.so
python -m pykingas -test