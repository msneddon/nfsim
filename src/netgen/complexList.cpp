/*
 * complexList.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */


#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
ComplexList::ComplexList ( )
{

}


// Destructor
ComplexList::~ComplexList ( )
{
}


// add a complex to the complex list
bool ComplexList::addComplexToList ( Complex * c)
{
	string label;

	// get label for complex
	label = getComplexLabel ( c );

	// check to see if this is new?

	// add complex to list
	label_map[label] = c;
	complex_list.push_back(c);

	return true;
}


// use label to get complex pointer
Complex * ComplexList::getComplexByLabel ( string& label )
{
	return label_map[label];
}


// canonical labeling method
string ComplexList::getComplexLabel ( Complex * c )
{
	string label;

	// put HNauty here

	// make use of these pre-initialized iterators:
	// map <string, Complex *>::iterator labelMapIter;
	// vector <Complex *>::iterator      complexIter;
	// vector <Molecule *>::iterator     molIter;

	return label;
}


// compare two complexes by label
int ComplexList::compare ( Complex * c1, Complex * c2 )
{
	return 0;
}
