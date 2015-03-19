
#ifndef NFOUTPUT_HH_
#define NFOUTPUT_HH_



#include "../NFcore/NFcore.hh"
#include "../NFutil/NFutil.hh"

namespace NFcore
{
	class MoleculeType;
	class System;

	class Outputter {
		public:
			Outputter(string filename, System *s);
			virtual ~Outputter();

			virtual void output()=0;
			virtual void outputHeader()=0;
			virtual void tryToDump(double simTime)=0;
		protected:
			string filename;
			System *s;
			ofstream outputFileStream;
	};



	//! Class for outputting the state of the system.
	/*!
	    This class is useful for outputting and dumping all the states of the molecules
	    in the system.  This is ideal for then determining complex output values such
	    as average aggregate size.  DumpSystem objects are initialized at most once per
	    system by the createSystemDumper() function defined in the NFinput.hh file.

	    to invoke:  use -dump flag as:
	    -dump [t1;t2;t3;...tn]->/optional/path/to/directory
	    where t1 to tn are the output times you want dumped.  you can also list the
	    output times in a matlab style sequence as:
	    -dump [t_start:step:t_end]
	    if no step is given, the step is assumed to be at one second intervals.  the
	    path to the output directory is optional, and if omitted the files will be dumped
	    to the current directory.

	    For each output time, a number of files will be created.  There will be one
	    header file that contains the number of molecules and other MoleculeType information
	    needed to correctly read the binary files.  Then, there will be one binary file
	    for each MoleculeType giving the details of each molecule as

	    moleculeId,complexId,state_1,partner_1, ... state_n,partner_n,localFuncValue_1 ... localFuncValue_n
	    the moleculeId is the unique identifier of the molecule.  The complexId gives the unique
	    id of the complex the molecule is in.  For each component, the state value and the binding
	    partner if any is given.  Finally, each local function's value is given.

	    If you don't want to worry about all of that, take a look at the NFanalysis matlab library
	    for tools that let you easily read the binary files into a nice data structure that can
	    be used to analyze output.

	 */
	class DumpSystem {

		public:
			//! The only constructor of DumpSystem objects
			/*!
			   Takes a reference to the system, a sorted list of double valued number corresponding
			   to the simulation times that a dump is required.  The relative path to the output
			   directory of all your folders, and a flag that tells us whether or not to output
			   random things.
			 */
			DumpSystem(System *s, vector <double> dumpTimes, string pathToFolder, bool verbose);

			//! Deconstructor that does nothing - DumpSystems don't tie up any heap memory.
			/*!
			*/
			~DumpSystem();

			//! Key function of the DumpSystem class that is called at each simulation step
			/*!
			  The System will call this function at each step, and the DumpSystem will determine
			  whether or not to dump its files based on the output times originally provided.
			*/
			void tryToDump(double simTime);

		protected:

			//! Protected function for dumping the header file at the given time step
			/*!
			*/
			void dumpHeaderFile(double dumpTime);

			//! Protected function for dumping the MoleculeType files at the given time step
			/*!
			*/
			void dumpMoleculeTypeFiles(double dumpTime);


			System *s; /*!< reference to the system being dumped*/
			int currentDumpTimeIndex; /*!< keeps track of when we have dumped so far */
			vector <double> dumpTimes; /*!< the list of the times to dump */
			string pathToFolder; /*!< gives the path to the output directory */
			bool verbose; /*!< flag that tells us if we should output all sorts of things */

	};




}








#endif /* NFOUTPUT_HH_ */
