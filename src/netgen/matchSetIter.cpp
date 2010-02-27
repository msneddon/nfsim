/*
 * matchSetIter.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */


#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
MatchSetIter::MatchSetIter ( ReactionClass * _rc )
{
	// cast as a BasicRxnClass so we can access reactantLists
	rc = (BasicRxnClass*)_rc;
	n_reactants = rc->getNumOfReactants();

	// initialize curr_set vector and reactantLists vectors
	int ii = 0;
	while ( ii < n_reactants )
	{
		curr_set.push_back( 0 );
		reactantLists.push_back( (rc->reactantLists)[ii] );
		++ii;
    }

	// set curr_set to zero and determine if any match sets exist
	reset();
}



// Destructor
MatchSetIter::~MatchSetIter ( )
{

}


// reset iterator to zero
void MatchSetIter::reset ( )
{
	// reset curr_set to zero vector
	more_sets = true;
	reactantList_iter = reactantLists.begin();
	for ( curr_set_iter = curr_set.begin(); curr_set_iter < curr_set.end(); curr_set_iter++ )
	{
		// set to the first mappingSet in reactantList
		(*curr_set_iter) = 0;

		// check if the reactantList is empty
		if (  (*reactantList_iter)->size() == 0  )
	    {   more_sets = false;   }
	    reactantList_iter++;
	}
	return;
}

// advance iterator
// WARNING: be sure reactantList_iter and curr_set_iter are at the end
//   of their respective vectors.
void MatchSetIter::advance ( )
{
	--reactantList_iter;  --curr_set_iter;

	if ( curr_set_iter < curr_set.begin() )
	{   more_sets = false;  }
	else
	{
	    // advance index
	    (*curr_set_iter)++;
	    // if index is out of bounds ...
        if ( *curr_set_iter >= (unsigned int)(*reactantList_iter)->size()  )
        {
    	    (*curr_set_iter) = 0;
    	    advance( );
        }
	}

	++reactantList_iter;  ++curr_set_iter;
}




// load pointers to next match set in the specified vector
bool MatchSetIter::nextMatchSet ( vector <MappingSet *> & match_set )
{
	if ( !more_sets ) {  return false;  }

	// clear out match_set vector
	match_set.clear();

    // copy mappingSet pointers corresponding to curr_set into next_set vector
	cout << "curr_set: ";
	reactantList_iter = reactantLists.begin();
	for ( curr_set_iter = curr_set.begin(); curr_set_iter < curr_set.end(); curr_set_iter++ )
	{
		// get mappingSet from reactantList corresponding to the index at curr_set_iter
        match_set.push_back( (*reactantList_iter)->getMappingSet( *curr_set_iter ) );

		cout << *curr_set_iter << " ";
        // advance iterators
        ++reactantList_iter;
	}
	cout << endl;
	// advance iterator
	//  ..first make sure iterators are at the end
    reactantList_iter = reactantLists.end();
    curr_set_iter = curr_set.end();
    advance( );
    // return with true value
    return true;
}


