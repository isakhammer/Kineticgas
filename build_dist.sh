#!/usr/bin/env bash

# Script to build and publish pypi-distribution
# TODO: Add code that bumps version number

bash cpp/build_mac.sh
python -m pykingas -test
status=$?
if [[ $status -ne 0 ]]; then
    echo "pykingas-tests failed with exit code $status"
    exit $status
else
    echo "pykingas-tests were successful"
fi
python setup.py sdist bdist_wheel
twine upload -r pypi dist/*