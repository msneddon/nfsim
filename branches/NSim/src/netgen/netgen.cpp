/*
 * netgen.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
Netgen::Netgen ( System * sys_ )
{
	sys = sys_;
	rc_iter.setReactionClassList ( &(sys->allReactions) );
}

// Destructor
Netgen::~Netgen ( )
{

}


// Network generation method
void Netgen::generate_network ( )
{
	cout << "generate_network!" << endl;
}

