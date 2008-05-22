#ifndef NFREACTIONS_HH_
#define NFREACTIONS_HH_

//Standard Libraries that we need
#include <vector>
#include <string>

//Include the core NFsim objects
#include "../NFcore/NFcore.hh"

//Include the local header files needed for reactions
#include "reactantLists/reactantList.hh"
#include "transformations/transformationSet.hh"
#include "transformations/transformation.hh"
#include "transformations/speciesCreator.hh"
#include "mappings/mapping.hh"
#include "mappings/mappingSet.hh"
#include "mappings/mappingGenerator.hh"



using namespace std;

namespace NFcore
{

	/*!
	 	Function used for initial testing of the entire "mapping" construct including MappingSets,
	 	TransformationSets, ReactantLists, etc.
	    @author Michael Sneddon
	 */
	void test();
	
	/*!
	 	Some very simple initial testing and debugging of ReactantLists and the rest of the
	 	mapping functions unter NFreactions
	    @author Michael Sneddon
	 */
	void test_simple();
}

#endif /*NFREACTIONS_HH_*/
