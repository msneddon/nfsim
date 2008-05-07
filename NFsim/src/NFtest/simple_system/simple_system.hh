#ifndef SIMPLE_SYSTEM_HH_
#define SIMPLE_SYSTEM_HH_


#include "../../NFcore/NFcore.hh"
#include "../../NFreactions/basic/NFreactions.hh"
#include "../../NFreactions/transformation/transformation.hh"

using namespace NFcore;

namespace NFtest_simple_system
{
	void run();
	
	
	
	MoleculeType * createX(System *s);
	MoleculeType * createY(System *s);
	
	

	ReactionClass * createReactionXDephos(MoleculeType *molX, double rate);
	ReactionClass * createReactionYphosX(MoleculeType *molX, MoleculeType *molY, double rate);
	ReactionClass * createReactionXYbind(MoleculeType *molX,MoleculeType *molY, double rate);
	ReactionClass * createReactionXYunbind(MoleculeType *molX, MoleculeType *molY, double rate);
	void addObs(System * s, MoleculeType *molX, MoleculeType *molY);
}









#endif /*SIMPLE_SYSTEM_HH_*/
