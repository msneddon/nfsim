#ifndef TRANSFORMATIONS_HH_
#define TRANSFORMATIONS_HH_



#include "../../NFcore/NFcore.hh"

#include "../../NFreactions/transformation/transformation.hh"
#include "../../NFreactions/mapping/mapping.hh"

using namespace NFcore;

namespace NFtest_transformations
{
	void run();
	
	MoleculeType * tc_createX(System * s, int count);
	MoleculeType * tc_createY(System * s, int count);
}

#endif /*TRANSFORMATIONS_HH_*/
