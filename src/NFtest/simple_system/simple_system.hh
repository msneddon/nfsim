#ifndef SIMPLE_SYSTEM_HH_
#define SIMPLE_SYSTEM_HH_



//The two header files listed here are all that we need to set up and
//run a simple simulation with basic reactions.
#include "../../NFcore/NFcore.hh"
#include "../../NFreactions/NFreactions.hh"
#include "../../NFreactions/reactions/reaction.hh"


//NFcore is the primary namespace of the NFsim program.  It contains
//all the basic classes and functions for creating, running, and
//outputting a given system.  Be sure to use this namespace!
using namespace NFcore;



//!  Namespace for the Simple System (abreviated ss) example.
/*!
	This namespace maintains functions that create, run, and output a simple
	enzymatic reaction type system.  If you want to learn the details of NFsim,
	this is the place to start looking.  Be sure to look at these files:
	simple_system.hh and simple_system.cpp .  Start with the run
	function in the cpp file, and follow the documentation to see how
	everything works.  Have fun!
	@author Michael Sneddon
*/
namespace NFtest_ss
{

	//!  Runs the simple enzymatic type reaction system as an example system.
	/*!
		Start here in simple_system.cpp to see how a simulation is created and run
		in the NFsim world.
		@author Michael Sneddon
	 */
	void run();



	/*!
		Creates Molecule of type X for the simple system.  In this system, X is the
		molecule that gets modified (in this case phosphorylated) by the enzyme.
		@author Michael Sneddon
	 */
	MoleculeType * createX(System *s);

	/*!
		Creates Molecule of type Y for the simple system.  In this system, Y is the
		enzyme that phosphorylates X.
		@author Michael Sneddon
	 */
	MoleculeType * createY(System *s);




	/*!
		Creates a simple dephosphorlyation reaction consisting of a single state change.  Look
		here first to get the basic idea of how reactions are defined and in particular, how
		simple unimolecular reactions are created.
		@author Michael Sneddon
	 */
	ReactionClass * createReactionXDephos(MoleculeType *molX, double rate);

	/*!
		Creates the binding reaction between X and Y.  Look here to see how bimolecular reactions
		such as binding reactions can be defined.
		@author Michael Sneddon
	 */
	ReactionClass * createReactionXYbind(MoleculeType *molX,MoleculeType *molY, double rate);

	/*!
		Creates an unbinding reaction between X and Y.  Look here to see how unbinding reactions
		can be defined.
		@author Michael Sneddon
	 */
	ReactionClass * createReactionXYunbind(MoleculeType *molX, MoleculeType *molY, double rate);

	/*!
		Creates the catalytic step of the enzymatic reaction in this simple system.  Look here
		to learn how to create reactions that include multiple transformations of the reactants.
		@author Michael Sneddon
	 */
	ReactionClass * createReactionYphosX(MoleculeType *molX, MoleculeType *molY, double rate);




	/*!
		Creates the observables used in the Simple System.  Look at this function to
		learn how basic observables can be created and added to the system.
		@author Michael Sneddon
	 */
	void addObs(System * s, MoleculeType *molX, MoleculeType *molY);
}


#endif /*SIMPLE_SYSTEM_HH_*/
