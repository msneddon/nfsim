#ifndef TESTCOMPARE_HH_
#define TESTCOMPARE_HH_



#include "../../NFcore/NFcore.hh"
using namespace NFcore;

namespace NFtest_compare
{
	void run();
	
	MoleculeType * tc_createX(System * s, int count);
	MoleculeType * tc_createY(System * s, int count);
		
}

#endif /*TESTCOMPARE_HH_*/
