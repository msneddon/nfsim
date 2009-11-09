/*
 * reactionClassIterator.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
ReactionClassIterator::ReactionClassIterator ( )
{
}


// Destructor
ReactionClassIterator::~ReactionClassIterator ( )
{
}


void ReactionClassIterator::setReactionClassList ( vector <ReactionClass *> * rc_list )
{
	reactionClassList = rc_list;
}


bool ReactionClassIterator::getNextReactionClass ( ReactionClass * rc )
{
	return false;
}

