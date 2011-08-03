/*
 * ComplexList.cpp
 *
 *  Created on: Mar 2, 2010
 *      Author: justin
 */


#include "NFcore.hh"
#include <math.h>
//#include <fstream>
//#include "../NFscheduler/NFstream.h"
//#include "../NFscheduler/Scheduler.h"

using namespace std;
using namespace NFcore;


// Constructor
ComplexList::ComplexList( )
{
	sys = 0;
	useComplex = false;
}


// Destructor
ComplexList::~ComplexList()
{
    // delete complex objects on the list
	Complex * c;
	while( allComplexes.size() > 0 )
	{
		c = allComplexes.back();
		allComplexes.pop_back();
		delete c;
	}
}


// core methods:
//   createComplex
//   getNextAvailableComplex
//   notifyThatComplexIsAvailable

int ComplexList::createComplex(Molecule * m)
{
	if (!useComplex) return -1;  //Only create complexes if we intend on using them...
	int c_id = allComplexes.size();
	Complex * c = new Complex(sys, c_id, m);
	allComplexes.push_back(c);
	return c_id;
}


// NOTE: there should always be at least as many complex objects in existence as molecules,
//  since molecules are constructed as singleton complexes and the complex is not deleted (rather, unassigned)
//  if the molecule becomes bound in a larger complex. Therefore, we should never run out of complexes
//  on the queue (so long as we stick to the rule of: create new complex every time a new molecule is instantiated.
Complex * ComplexList::getNextAvailableComplex()
{
	Complex * c = allComplexes.at(nextAvailableComplex.front());
	nextAvailableComplex.pop();
	return c;
}



void ComplexList::notifyThatComplexIsAvailable(int ID_complex)
{
	nextAvailableComplex.push(ID_complex);
}






void ComplexList::purgeAndPrintAvailableComplexList()
{
	cout << "AvailableComplexes:";
	while(	!nextAvailableComplex.empty() )
	{
		cout << " -> " << nextAvailableComplex.front();
		nextAvailableComplex.pop();
	}
	cout << endl;
}



void ComplexList::printAllComplexes()
{
	cout<<"All System Complexes:"<<endl;
	// TODO: debug!
	for( complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
		(*complexIter)->printDetailsLong();
	cout<<endl;
}



void ComplexList::outputComplexSizes(double cSampleTime)
{
	int size = 0;
	(sys->getOutputFileStream())<<"\t"<<cSampleTime;
	for( complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getComplexSize();
		if (size!=0) (sys->getOutputFileStream())<<"\t"<<size;
	}
	(sys->getOutputFileStream())<<endl;
}



double ComplexList::outputMeanCount(MoleculeType *m)
{
	int count = 0;
	int sum = 0;
	int allSum = 0;
	int allCount=0;
	int size=0;
	(sys->getOutputFileStream())<<"\t"<<sys->getCurrentTime();
	for( complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);
		if(size>=2) { count++; sum+=size;}
		if(size>=1) { allSum += size; allCount++; }

	}
	//cout<<sum<<"/"<<count<<"   "<<allSum<<"/"<<allCount<<endl;
	if(count!=0)
	{
		(sys->getOutputFileStream())<<"\t"<<((double)sum/(double)count)<<endl;
		return ((double)sum/(double)count);
	}
	else
	{
		(sys->getOutputFileStream())<<"\t"<<0.0<<endl;
		return 0.0;
	}

	return ((double)sum/(double)count);
}



double ComplexList::calculateMeanCount(MoleculeType *m)
{
	int count = 0;
	int sum = 0;
	int allSum = 0;
	int allCount=0;
	int size=0;

	for( complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);
		if(size>=2) { count++; sum+=size; }
		if(size>=1) { allSum+=size; allCount++; }
	}
	return ((double)sum/(double)count);
}



void ComplexList::outputMoleculeTypeCountPerComplex(MoleculeType *m)
{
	int size = 0;
	(sys->getOutputFileStream())<<"\t"<<(sys->getCurrentTime());
	for( complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);

		if(size>=1) (sys->getOutputFileStream())<<"\t"<<size;
	}
	(sys->getOutputFileStream())<<endl;

}



// TODO: figure out how friend functions work!!
// friend functions
template<class T>
NFstream& operator<<(NFstream& nfstream, const T& value)
{
    if (nfstream.useFile_)
	nfstream.file_ << value;
    else
	nfstream.str_ << value;

    return nfstream;
}

