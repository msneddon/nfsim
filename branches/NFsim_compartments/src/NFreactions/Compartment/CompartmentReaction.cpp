/*
 * CompartmentReaction.cpp
 *
 *  Created on: Apr 13, 2009
 *      Author: kstanton
 */

#include "CompartmentReaction.h"

unsigned int CompartmentReaction::nCompartments = 0;
CompartmentReaction::CompartmentReaction(string name, double baseRate, TransformationSet *transformationSet) :
	ReactionClass(name,baseRate,transformationSet)
{
	this->reactionType = COMPARTMENT_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	if(nCompartments <1)
	{
		throw "Number of compartments not initialized";
	}

	for(unsigned int ii=0; ii<nCompartments; ii++)
	{
		// create the compartment
		m_mapCompartmentList.insert(
				pair<unsigned int, Compartment*>(
						ii,
						new Compartment(ii,baseRate,transformationSet,n_reactants,this)));
	}

}

CompartmentReaction::~CompartmentReaction()
{
	if(DEBUG) cout<<"Destorying rxn: "<<name<<endl;
	static map<unsigned int,Compartment*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		delete cmpIter->second;
	}
}


inline bool CompartmentReaction::tryToAdd(Molecule *m, unsigned int reactantPos)
{

	//First a bit of error checking, that you should skip unless we are debugging...
	if(DEBUG)
	{
		if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
		{
			cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
			exit(1);
		}
	}
	// gather a little background on this molecule and reaction pair

	//look up the index
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//Find the reaction ID of this molecule in this reaction if it is in this reaction at all if not -1
	int rxnId = m->getRxnListMappingId(rxnIndex);
	//Find the compartment that this Molecule is currently in in this reaction  if not -1
	int oldCompartmentId = m->getRxnListCompartmentMappingId(rxnIndex);
	// New Compartment is taken from the transformed molecule
	int newCompartmentId = m->getCompartmentId();

	// first off this is a compartment reaction and must have a compartment ID if its in the molecule
	// if its not in a compartment
	if(oldCompartmentId < 0)
	{
		if(rxnId<0) // and its not in this reaction
		{
			// then lets initialize it to the compartment that it is located in
			oldCompartmentId = newCompartmentId;
		}
		else // oops no compartment but its in this reaction
		{
			throw "No compartment RxnListCompartmentMappingId associated with this compartment reaction for which there is a RxnListMappingId";
		}
	}
	// if either the source or target compartment do not exist
	if(m_mapCompartmentList.find(oldCompartmentId) == m_mapCompartmentList.end()
			||m_mapCompartmentList.find(newCompartmentId) == m_mapCompartmentList.end()
			|| m_mapCompartmentList[oldCompartmentId] == NULL
			|| m_mapCompartmentList[newCompartmentId] == NULL)
	{
		throw "Dereferencing a null compartment!";
	}

	// remove clones first
	m_mapCompartmentList[oldCompartmentId]->removeClones(m,reactantPos,rxnIndex);
	m_mapCompartmentList[m->getCompartmentId()]->removeClones(m,reactantPos,rxnIndex);

	//Here we get the standard update...
	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		m_mapCompartmentList[oldCompartmentId]->tryToRemove(m, reactantPos, rxnIndex);
	} else {
		m_mapCompartmentList[m->getCompartmentId()]->tryToAdd(m, reactantPos, rxnIndex);
	}
	return true;
}


double CompartmentReaction::update_a()
{
	a = 0;
	static map<unsigned int,Compartment*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		a+=cmpIter->second->update_a();
	}
	return a;
}


void CompartmentReaction::init()
{
	static map<unsigned int,Compartment*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		cmpIter->second->init(this);
	}
}

void CompartmentReaction::prepareForSimulation()
{

}

void CompartmentReaction::remove(Molecule *m, unsigned int reactantPos)
{

	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}

	// get the reaction index
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	// check to make sure this molecule is in this reaction and that the compartment is also filled out
	if((m->getRxnListMappingId(rxnIndex)>=0) && (m->getRxnListCompartmentMappingId(rxnIndex)>=0))
	{
		m_mapCompartmentList[m->getRxnListCompartmentMappingId(rxnIndex)]->tryToRemove(m, reactantPos, rxnIndex);
	}
}

void CompartmentReaction::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	cerr<<"You are trying to notify a Basic Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}
unsigned int CompartmentReaction::getReactantCount(unsigned int reactantIndex) const
{
	unsigned int nReactants = 0;
	static map<unsigned int,Compartment*>::const_iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		nReactants = nReactants + cmpIter->second->getReactantCount(reactantIndex);
	}
	return nReactants;
}
void CompartmentReaction::printFullDetails() const
{
	cout<<"BasicRxnClass: "<<name<<endl;
	static map<unsigned int,Compartment*>::const_iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		cmpIter->second->printDetails();
	}
}
void CompartmentReaction::pickMappingSets(double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//So, we shall loop through the lists and extract out the MappingSets and
	//the molecule ID that the mappingSet was generated from
	//unsigned int *moleculeIDs = new unsigned int [n_reactants];

	// choose a compartment based on its propensity in this reaction
	double a_sum=0, last_a_sum=0;
	Compartment* chosenComp = 0;
	double leftover_random_A_number = 0;

	//WARNING - DO NOT USE THE DEFAULT C++ RANDOM NUMBER GENERATOR FOR THIS STEP
	// - IT INTRODUCES SMALL NUMERICAL ERRORS CAUSING THE ORDER OF RXNS TO
	//   AFFECT SIMULATION RESULTS
	static map<unsigned int,Compartment*>::const_iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		a_sum += cmpIter->second->a;
		if (random_A_number <= a_sum && chosenComp==0)
		{
			chosenComp = cmpIter->second;
			leftover_random_A_number = random_A_number-last_a_sum;
			break;
		}
		last_a_sum = a_sum;
	}
	if(chosenComp==0)
	{
		cerr<<"Error: random_A_number exceeds a_sum!!!"<<endl;
		cerr<<"randNum: "<<random_A_number<<"  a_sum: "<< a_sum<<" running a_tot:"<<a<<endl;
	}
	chosenComp->pickMappingSets(mappingSet, leftover_random_A_number);
}

void CompartmentReaction::restrictToCompartment(unsigned int compartmentId)
{
	static map<unsigned int,Compartment*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		cmpIter->second->active = false;
	}
	m_mapCompartmentList[compartmentId]->active = true;
}
