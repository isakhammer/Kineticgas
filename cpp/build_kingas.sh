#!/usr/bin/env bash
set -e

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
            bash cpp/build_integration.sh --Debug
            cd cpp/debug
            cmake -DCMAKE_BUILD_TYPE=Debug ..
            make
            cd ../..
            DIR=/tmp; [ -d "$DIR" ] && rm -rf tmp
            mkdir tmp
            mv pykingas/KineticGas_d.so tmp/KineticGas_d.so
            cp cpp/debug/KineticGas_d.cpython-39-darwin.so pykingas/KineticGas_d.so
            python -m pykingas -test -debug
            status=$?
            if [[ $status -eq 0 ]]; then
                rm -rf tmp
            else
                mv tmp/KineticGas_d.so pykingas/KineticGas_d.so
                rm -rf tmp
            fi
            exit $status
            ;;
        -*)
            echo "Bad option $1"
            exit -1
    esac
    shift
done

echo "Building Release"
bash cpp/build_integration.sh
cd cpp/release
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../..
DIR=/tmp; [ -d "$DIR" ] && rm -rf tmp
mkdir tmp
mv pykingas/KineticGas_r.so tmp/KineticGas_r.so
cp cpp/release/KineticGas_r.cpython-39-darwin.so pykingas/KineticGas_r.so

python -m pykingas -test -release
status=$?
echo "Test Exited with status $status"
if [[ $status -eq 0 ]]; then
    echo "Overwriting .so file"
    rm -rf tmp
else
    mv tmp/KineticGas_r.so pykingas/KineticGas_r.so
    rm -rf tmp
fi
exit $status