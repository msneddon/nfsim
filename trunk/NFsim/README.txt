

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%     @@    @  @@@@@                %
%     @ @   @  @                    %
%     @  @  @  @@@@  ___            %
%     @   @ @  @    /__  | |\ /|    %
%     @    @@  @    ___\ | | v |    %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NFsim - the network free stochastic simulator, v1.07

michael w. sneddon
james r. faeder
thierry emonet

Yale University
University of Pittsburgh
funded by the National Science Foundation


################################################################################

NFsim is released under the GNU General Public License.  See LICENSE.txt for
more details about redistribution restrictions.
  
For help with running NFsim, see the user manual, NFsim_manual_[version].pdf,
and open the example model "simple_system.bngl".

Executable files for Windows, Mac and Linux are in the "bin" directory.  Source
code and makefiles for NFsim are in the NFcode directory.  Example models are
in the "models" directory, with README files.  BioNetGen and the ODE and SSA
solvers used by BioNetGen are in the BNG and Network2 directories.  Finally,
a suite of helpful analysis and other modeling tools can be found in the
"NFtools" directory.

Enjoy your new network-free world!



################################################################################

Release Notes

v1.08   Dec, 2010 - With the new TotalRate keyword, users are now able to
        specify whether or not to use the microscopic (default) interpretation
        or macroscopic (TotalRate) interpretation of rate laws.  Now, NFsim
        convention matches BNG.  Previously, NFsim interpreted all rates as
        microscopic except for global functions, which were interpreted as
        macroscopic.  This is now also explained in the user manual.  
        
        Although this was the main update, there were a few smaller items that
        are updated as well.  First, users now have the option of outputting
        gdat files in a comma delimited format (csv), which makes parsing the
        output file easier in some circumstances, using the flag "-csv".



v1.07   Nov, 2010 - A series of updates to the code were made in this
        release.  (1) RNF files that are not found produce an error message.
        Previously, no error message was given and execution proceeded as if
        the RNF flag was not given.  (2) Input flags to NFsim can be given
        in the original format with a single dash (as in ./NFsim -logo), or
        with the more "linuxy" style double dash (as in ./NFsim --logo).
        (3) The parameter scan script had problems when parsing BNGL files
        with local functions due to the '%' character.  This is now fixed.
        (4) The universal traversal limit is automatically set to be the size
        of the largest pattern in the system, which can be overridden by
        passing the -utl flag.  This allows users who are unfamiliar with
        this speedup to still take advantage of it to a certain extent.  (5)
        the -rtag flag was added that allows NFsim to produce output whenever
        a particular reaction, given by the -rtag flag, is given.  This allows,
        for instance, users to track the fates of single particles exactly
        without using the comprehensive molecule output feature.  (6) The above
        changes are documented in an updated user manual.


v1.06   Sept 28, 2010 - added scripts for running NFsim from Matlab, parameter
        scanning, and basic parameter estimation.  The manual is also updated
        to reflect these changes.  However, the precompiled executables of NFsim
        remain unchanged from v1.05, so running them will give the old version 
        number unless you recompile the code on your own computer.  Also, models
        that were used to compare the performance of NFsim to DYNSTOC, RuleMonkey,
        and Kappa are now included with a readme file under: 
        models/performance_test_models.

v1.052  First publically released stable build




