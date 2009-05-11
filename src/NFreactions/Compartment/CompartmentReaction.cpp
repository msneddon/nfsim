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

	// if the compartment doesnt exists
	if(m_mapCompartmentList.find(m->getCompartmentId()) == m_mapCompartmentList.end())
	{
		throw "Dereferencing a null compartment!";
	}
	// Try to add the molecule to the compartment.
	return m_mapCompartmentList[m->getCompartmentId()]->tryToAdd(m,reactantPos);
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


	//Get the specified reactantList
	ReactantList *rl = m_mapCompartmentList[m->getCompartmentId()]->reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);


	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
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

