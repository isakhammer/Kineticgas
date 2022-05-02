#!/usr/bin/env bash
set -e
cd cpp/Integration/release

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
            cd ../../..
            DIR=/tmp; [ -d "$DIR" ] && rm -rf tmp
            mkdir tmp
            mv integration/Integration_d.so tmp/Integration_d.so
            cp cpp/Integration/debug/Integration_d.cpython-39-darwin.so integration/Integration_d.so
            python integration/integrator_unittests.py -d
            status=$?
            if [[ $status -eq 0 ]]; then
                rm -rf tmp
            else
                mv tmp/Integration_d.so integration/Integration_d.so
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
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
cd ../../..

DIR=/tmp; [ -d "$DIR" ] && rm -rf tmp
mkdir tmp
mv integration/Integration_r.so tmp/Integration_r.so
cp cpp/Integration/release/Integration_r.cpython-39-darwin.so integration/Integration_r.so
python integration/integrator_unittests.py
status=$?
if [[ $status -eq 0 ]]; then
    rm -rf tmp
else
    mv tmp/Integration_r.so integration/Integration_r.so
    rm -rf tmp
fi
exit $status