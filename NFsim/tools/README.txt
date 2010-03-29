

This directory contains other software tools that are useful for either running
or analyzing NFsim simulations and models.

#######################################################################################

runbng.m  -  This is a Matlab script that allows you to run NFsim simulations
             through Matlab, and automatically plot the results.  This script will
             be useful for users who are not comfortable with the command-line, but
             have experience with Matlab commands.  See the Matlab script help for
             more documentation on how to run it.

             
readNFsimBinary.m  - This Matlab function will read in NFsim binary output files that
                     were generated with the '-b' command-line option.  This function
                     reads the binary file and returns the data and headers in the same
                     way that the tblread command can read GDAT files.  See the function
                     help for more documentation on how to run it.
                     
                     
NFanalyzeDump      - This directory contains a set of Matlab tools for reading and
                     processing the output of NFsim dumps, which output the entire
                     state of the system.  To get started using these tools, first see
                     the help section of the "readNFdump.m" function.  This function
                     will read an entire directory of dump files and return a struct
                     that can be processed by the other methods.  These tools are
                     necessary, for instance, to extract the value of local functions
                     in NFsim, plot average aggregate sizes, or determine the structure
                     and configuration of polymers.