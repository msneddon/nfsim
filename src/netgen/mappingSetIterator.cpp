/*
 * mappingSetIterator.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */


#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
MappingSetsIterator::MappingSetsIterator ( ReactionClass * rc_ )
{
	rc = rc_;
}


// Destructor
MappingSetsIterator::~MappingSetsIterator ( )
{

}


// load pointers to next MappingSets in the specified vector
bool MappingSetsIterator::getNextMappingSets ( vector <MappingSet *> * next_set )
{
	return true;
}

