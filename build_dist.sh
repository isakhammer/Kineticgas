#!/usr/bin/env bash

# Script to build and publish pypi-distribution
# TODO: Add code that bumps version number

echo "This script builds and distributes to pypi! Comment out the next line when you are ready to distribute."
exit 0

bash cpp/build_kingas.sh --fullclean
status=$?
if [[ $status -ne 0 ]]; then
    echo "pykingas tests failed with exit code $status"
    exit $status
else
    bash cpp/build_kingas.sh --Debug
    status=$?
    if [[ $status -ne 0 ]]; then
        echo "pykingas -debug tests failed with exit code $status"
        exit $status
    else
        echo "pykingas -debug tests were successful"
    fi
fi

python setup.py sdist bdist_wheel

while test $# -gt 0
do
    case "$1" in
        -dist)
            echo "Publishing to PyPI!"
            twine upload -r pypi dist/*
            status=$?
            exit $status
            ;;
        -*)
            echo "Invalid option $1"
            exit -1
    esac
    shift
done
echo "Publishing to TestPyPI"
twine upload --repository testpypi dist/*