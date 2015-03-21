

Nauty v2.4 is an open-source package used for canonical graph labeling in the populations feature of NFsim.
see: http://cs.anu.edu.au/~bdm/nauty/

This code has been patched by renaming the typedef 'set' to 'nset', as set conflicts
with the c++ std::set and generates compiler errors in the most recent releases of
x code on OSX. Earlier versions of x code and current versions of Linux gcc have no
such proble.  There is likely a better fix, but as a quick patch to get things
running in OSX, this works perfectly fine.

