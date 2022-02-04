# Script to build and publish pypi-distribution
# TODO: Add code that bumps version number

python setup.py sdist bdist_wheel
twine upload -r pypi dist/*