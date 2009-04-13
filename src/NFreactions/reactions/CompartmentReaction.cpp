/*
 * CompartmentReaction.cpp
 *
 *  Created on: Apr 13, 2009
 *      Author: kstanton
 */

#include "CompartmentReaction.h"
#include "reaction.hh"


using namespace std;
using namespace NFcore;

CompartmentReaction::CompartmentReaction() {
	// TODO Auto-generated constructor stub

}

CompartmentReaction::~CompartmentReaction() {
	// TODO Auto-generated destructor stub
}
inline bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking, that you should skip unless we are debugging...
//	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
//	{
//		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
//		exit(1);
//	}


	//Get the specified reactantList
	rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//cout<<"got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;


	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		if(!reactantTemplates[reactantPos]->compare(m)) {
		//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}

	} else {
		//Try to map it!
		ms = rl->pushNextAvailableMappingSet();
		if(!reactantTemplates[reactantPos]->compare(m,ms)) {
			rl->popLastMappingSet();
			//we just pushed, then popped, so we a has not changed...
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}

	return true;
}
