/*
 * reactionClassIter.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
ReactionClassIter::ReactionClassIter ( )
{

}


// Destructor
ReactionClassIter::~ReactionClassIter ( )
{

}


// set system method
void ReactionClassIter::setSystem ( vector <ReactionClass *> * _reactionClasses )
{
    reactionClasses = _reactionClasses;
    for ( rc_iter = reactionClasses->begin();
	     	rc_iter < reactionClasses->end(); rc_iter++ )
    {
	    MatchSetIter * ms_iter = new MatchSetIter( *rc_iter );
	    matchSetIters.push_back( ms_iter );
    }
    reset();
}


// reset iterator to first reaction class
void ReactionClassIter::reset ( )
{
	matchSetIters_iter = matchSetIters.begin();
}


// get pointer to matchSetIter for next reaction class
MatchSetIter * ReactionClassIter::nextReactionClass (  )
{
	MatchSetIter * match_set_iter = 0;
	if ( matchSetIters_iter < matchSetIters.end() )
	{
		match_set_iter = (*matchSetIters_iter);
		++matchSetIters_iter;
	}
	return match_set_iter;
}

