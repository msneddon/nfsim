#ifndef TRANSCRIPTION_HH_
#define TRANSCRIPTION_HH_



//The two header files listed here are all that we need to set up and
//run a simple simulation with basic reactions.
#include "../../NFcore/NFcore.hh"
#include "../../NFreactions/NFreactions.hh"
#include "../../NFreactions/reactions/reaction.hh"


//NFcore is the primary namespace of the NFsim program.  It contains
//all the basic classes and functions for creating, running, and
//outputting a given system.  Be sure to use this namespace!
using namespace NFcore;



//!  Namespace for the Transcription example.
/*!

	@author Michael Sneddon
*/
namespace NFtest_transcription
{

	//!
	/*!

		@author Michael Sneddon
	 */
	void run();



	/*!

		@author Michael Sneddon
	 */
	MoleculeType * createRNA(System *s);




	/*!

		@author Michael Sneddon
	 */
	ReactionClass * createReactionRNAdegrades(MoleculeType *molRNA, double rate);



	/*!
		@author Michael Sneddon
	 */
	ReactionClass * createReactionRNAtranscribed(MoleculeType *molRNA, double rate);




	/*!
		@author Michael Sneddon
	 */
	void addObs(System * s, MoleculeType *molRNA);
}


#endif /*TRANSCRIPTION_HH_*/
