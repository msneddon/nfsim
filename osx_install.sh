if [ "${TRAVIS_OS_NAME}" = "linux" ]; then
    uname -a
else
    uname -a
    brew update
    brew install python
    brew tap homebrew/science
    brew install libxml2
    brew link --overwrite python
    brew unlink python && brew link python
    pip install numpy
fi
