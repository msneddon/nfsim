# `NFsim` – the network free stochastic simulator

This version is forked from https://github.com/msneddon/nfsim.
We also took the `CMakeLists.txt` file from  https://github.com/RuleWorld/nfsim  for compiling using CMake (and not having to use Eclipse for auto-generating `Make` files).

To install `NFsim`, clone this repository to your =~/git= repo:

```bash
git clone https://github.com/rasilab/nfsim.git
```

Then compile by (see https://github.com/RuleWorld/nfsim/blob/master/.travis.yml#L41):

```bash
mkdir build
cd build
cmake ..
make
```

You need cmake, make and gcc for above. gcc-9.3.0 in Ubuntu 20.04 worked without any issues for me.

Once it is compiled, add the binary `NFsim` in the `build` directory to your `PATH`.
For example, I add the following to my `~/.bashrc` script:

```bash
export PATH=$PATH":/home/rasi/git/nfsim/build"
```
