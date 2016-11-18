if [ "${TRAVIS_OS_NAME}" = "linux" ]; then
    uname -a
else
    uname -a
    brew update
    brew install python
    brew tap homebrew/science
    brew install libxml2
    brew link --overwrite python
    /usr/local/bin/pip install numpy
fi