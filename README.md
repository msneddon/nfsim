# NFsim - the network free stochastic simulator, v1.12.1

[![Build Status](https://travis-ci.org/RuleWorld/nfsim.svg?branch=master)](https://travis-ci.org/RuleWorld/nfsim)
[![Build status](https://ci.appveyor.com/api/projects/status/423rdsbq0dak9gtt?svg=true)](https://ci.appveyor.com/project/jczech/nfsim)
<a href="https://scan.coverity.com/projects/nfsim">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/15734/badge.svg"/>
</a>

- michael w. sneddon
- justin s. hogg
- jose-juan tapia
- james r. faeder
- thierry emonet

Yale University  
University of Pittsburgh  
funded by the National Science Foundation  

## Overview

NFsim is a free, open-source, biochemical reaction simulator designed to handle
systems that have a large or even infinite number of possible molecular
interactions or states. NFsim also has advanced and flexible options for
simulating coarse-grained representations of complex nonlinear reaction
mechanisms.

NFsim is ideal for modeling polymerization, aggregation, and cooperative
reactions that cannot be handled with traditional stochastic or ODE simulators.
Models are specified in the BioNetGen Langauge, providing a powerful model
building environment.

If you just want to download and use NFsim, you should simply download a
preconfigured packaged release from http://emonet.biology.yale.edu/nfsim. If
you want to hack on the code or make contributions, please create a fork and
submit pull requests to the dev branch.

If you use NFsim for your research or work, please cite NFsim as: Sneddon MW,
Faeder JR & Emonet T. Efficient modeling, simulation and coarse-graining of
biological complexity with NFsim. Nature Methods,(2011) 8(2):177-83.

## Repository Contents

NFsim is released under the MIT License. See LICENSE.txt for more details about
redistribution restrictions.
  
For help with running NFsim, see the user manual, NFsim\_manual\_[version].pdf,
and open the example model "simple\_system.bngl".

Source code is in the "src" directory. Example models are in the "models"
directory, with README files. BioNetGen and the ODE and SSA solvers used by
BioNetGen are in the BNG and Network2 directories.

Enjoy your new network-free world!

## Building NFsim

### Linux and OSX

Make sure you have a recent version of CMake installed and run the following at
a terminal:

    mkdir build
    cd build
    cmake ..
    make

### Windows

Make sure you have a recent version of CMake, Cygwin (or MinGW), and Ninja
installed. Then run the following at a terminal:

    mkdir build
    cd build
    cmake -G "Ninja"
    ninja

If you've built it with Cygwin and you want to move or package up the
executable, you'll need to copy the the following DLLs along with it:

* cygwin1.dll
* cygstdc++-6.dll
* cygz.dll
* cyggcc_s-seh-1.dll.

## Download Latest Test Builds

These builds are the from the head of master and are not guaranteed to be
stable. Use at your own risk.

* [Linux](https://bintray.com/jczech/nfsim/download_file?file_path=NFsim-linux)
* [OSX](https://bintray.com/jczech/nfsim/download_file?file_path=NFsim-osx)
* [Windows x64](https://bintray.com/jczech/nfsim/download_file?file_path=NFsim-Win64.zip) 
* [Windows x86](https://bintray.com/jczech/nfsim/download_file?file_path=NFsim-Winx86.zip)


## Release Notes

### v1.12.1 Aug, 2016

(a) Bugfix release. Addresses an error dealing with local function and species
labels. The error dealt with the way mappingSets where created and passed to
the localFunction evaluation. A test case (v19.bngl) was added that addresses
this case.

### v1.12   Dec, 2015

(a) Changes to how molecule instances are mapped to BasicRxn's
(BasicRxnClass::tryToAdd()). It was possible for certain kinds of rules that
the mappingSets were not updated correctly because the head molecule matched a
reactant pattern before and after a reaction event BUT the mapping was
different after the firing. In such cases, the mappingSet was not updated
properly. The logic was changed to fix this problem. (b) Further changes to
how molecule instances are mapped to a ReactionClass object (BasicRxn and
DORReactions). In particular it is often the case that graph symmetry leads to
a complex being able to map to a ReactionClass multiple times. Symmetry
considerations were being made on the reaction center but not on the context
components which led to an undercounting of the number of times a pattern agent
could match a rule instance in some edge cases (see v17.bngl in the validation
suite). This led to incorrect results or even NFSim crashes. NOTE: the changes
made in points (a) and (b) may cause some models to execute less efficiently.
(c) Fixed index bound checking in MoleculeType::getComponentStateName(). (d)
Updated the validation models to reflect current BNGL formatting standards.
Added a few new validation models that address the bugfixes included in this
version. This update also includes a Python version of the validation script.

### v1.11   Oct, 2012

(a) Molecules without components may be treated as population variables
rather than individual agents. This feature is useful for reducing memory
requirements when a simple molecule has very large population. A molecule
type may be flagged for treatment as a population variable by using the
"population" keyword following the molecule type definition in the BNGL
model file. See section 8.i in the documentation for further detail.
(b) Added a new Reaction Class called "FunctionProduct" that permits local
functions defined on two reactants. The local rate law must have the form
f(x)\*g(y), where x and y are tags on two distinct reactants. See section
7.c in the documentation for complete information.
(c) Fixed some problems evaluating complex-scoped local functions (CSLF). 
CSLF were not updated properly after reactions that split complexes or 
deleted molecules. As in v1.10, complex-scoped local functions are
enabled by default. A new command-line switch, -nocslf, has been added
which disables complex-scoped evaluation. 
(d) Improved efficiency for matching patterns with connected-to syntax
when the connected-to component does not have reaction center. This may
be especially notable in models with large complexes, e.g. polymerization.
		
### v1.10   Aug, 2011
        
(a) Command line parser now detects arguments that are not properly preceeded
by a dash, and generates a warning. (b) Includes a check when creating
template molecules that throws an error when users attempt to use Null or Trash
in reactant or observable patterns (anything that requires the creation of a
Template molecule). (c) fixed a bug introduced in v1.09 whereby a site was
allowed to bind to itself, for instance, in a dimerization rxn. (d) Support
for creating a new molecule bound to an existing molecule, as in a rule like
A(a) -> A(a!1).A(a!1). Existing code that implemented this feature did not
function properly with the check for null conditions before reactions were
fired. (e) Fixed bug in template molecule when clearing molecules after a
connected-to syntax search. In some cases, not all molecules were being
cleared, giving rise to situations where adding one observable created dangling
matches which affected the results of other observables. (f) NFsim is now
packaged with Network3, an updated version of the run\_network code to execute
ODE and SSA simulations. Network3 allows global functions in BNGL models among
other release features given here:
http://bionetgen.org/index.php/Release\_Notes. Note that Network3 does not
support On-The-Fly Stochastic Simulation (you will have to recompile Network2
to use this feature). (g)  Mac 32bit is no longer supported by NFsim, but you
can make executables for older Macs by recompiling the code on your own
machine. See the manual for instructions.

### v1.09   Apr, 2011

(a) NFsim now allows the mixing of integers and strings as component labels,
although if numbers and strings are mixed, all labels are parsed as strings,
NOT integers. Therefore, PLUS and MINUS keywords cannot be used if mixing
integer states and string states, and a warning will be generated if a state is
set to PLUS. One can only use PLUS or MINUS when ALL states are an integer
value greater than zero. This new behavior was needed to handle BNGL files that
used the convention of ~P specifying phosphorylated, and ~0 specifying
unphosphorylated.  (b) Fixed bug whereby if verbose option was turned on
without specifying an output file location, no output would be generated. Now,
output to a gdat file will be generated in these cases.  (c) Use of 'ss' input
argument to 'saveSpecies', which prints a a list of all species at the end of a
simulation, is now handled. This feature was implemented to allow future
support in BNG by restarting an NFsim simulation after it ends, which can be
done by parsing the output species list together with a BNGL model file. This
feature still has to be tested in BNG, and will likely be fully documented in
v1.10. The 'ss' flag writes the file to either system_name_nf.species, or a
file designated by the user.  (d) fixed memory leak in TemplateMolecules that
caused memory / performance issues with molecules having multiple identical
sites and a high degree of aggregation. (e) fixed csv error, where when the csv
flag is used, the header line is not comma delimited. (f)  nfsim now supports
intra-molecular binding. Previously these events were rejected as null events.

### v1.08 Dec, 2010

With the new TotalRate keyword, users are now able to specify whether or not to
use the microscopic (default) interpretation or macroscopic (TotalRate)
interpretation of rate laws. Now, NFsim convention matches BNG. Previously,
NFsim interpreted all rates as microscopic except for global functions, which
were interpreted as macroscopic. This is now also explained in the user manual.
Example models for the flagellar motor and oscillating gene expression have
been updated correspondingly so that they still produce the same results as in
the NFsim paper.

Users also now have the option of outputting gdat files in a comma delimited
format (csv), which makes parsing the output file easier in some circumstances,
using the flag "-csv". Additionally, a bug in the parameter scanning script was
fixed that caused the script to crash when scanning a model that includes the
local function syntax. 

### v1.07   Nov, 2010

A series of updates to the code were made in this release. (1) RNF files that
are not found produce an error message. Previously, no error message was given
and execution proceeded as if the RNF flag was not given. (2) Input flags to
NFsim can be given in the original format with a single dash (as in ./NFsim
-logo), or with the more "linuxy" style double dash (as in ./NFsim --logo).
(3) The parameter scan script had problems when parsing BNGL files with local
functions due to the '%' character. This is now fixed. (4) The universal
traversal limit is automatically set to be the size of the largest pattern in
the system, which can be overridden by passing the -utl flag. This allows
users who are unfamiliar with this speedup to still take advantage of it to a
certain extent. (5) the -rtag flag was added that allows NFsim to produce
output whenever a particular reaction, given by the -rtag flag, is given. This
allows, for instance, users to track the fates of single particles exactly
without using the comprehensive molecule output feature. (6) The above changes
are documented in an updated user manual.

### v1.06   Sept 28, 2010

Added scripts for running NFsim from Matlab, parameter scanning, and basic
parameter estimation. The manual is also updated to reflect these changes.
However, the precompiled executables of NFsim remain unchanged from v1.05, so
running them will give the old version number unless you recompile the code on
your own computer. Also, models that were used to compare the performance of
NFsim to DYNSTOC, RuleMonkey, and Kappa are now included with a readme file
under: models/performance\_test\_models.

### v1.052  

First publicly released stable build
