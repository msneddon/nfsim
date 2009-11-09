/*
 * reactionList.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */



#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
ReactionList::ReactionList ( )
{

}


// Destructor
ReactionList::~ReactionList ( )
{

}


// add a reaction to the list
bool ReactionList::addReaction ( Reaction * rxn )
{
	reaction_list.push_back( rxn );
	return true;
}
