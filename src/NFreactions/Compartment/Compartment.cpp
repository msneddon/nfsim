/*
 * Compartment.cpp
 *
 *  Created on: Apr 16, 2009
 *      Author: kstanton
 */

#include "Compartment.h"

Compartment::Compartment(
		unsigned int compId,
		double dbBaseRate,
		TransformationSet* transformationSet,
		unsigned int nReactants,
		CompartmentReaction* parentReaction)
{
	this->name = "No name assigned to compartment";
	reactantLists = new ReactantList *[nReactants];
	n_reactants = nReactants;
	m_pParentReaction = parentReaction;
	m_iCompId = compId;
		//Set up the reactantLists
	for(unsigned int r=0; r<nReactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
}
Compartment::~Compartment()
{
	if(DEBUG) cout<<"Destorying rxn Compartment: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		delete reactantLists[r];
	}
	delete [] reactantLists;
}
unsigned int Compartment::getReactantCount(unsigned int reactantIndex) const
{
	return reactantLists[reactantIndex]->size();
}
void Compartment::printDetails() const
{
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists[i]->printDetails();
}
void Compartment::pickMappingSets(MappingSet** mappingSet, double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//So, we shall loop through the lists and extract out the MappingSets and
	//the molecule ID that the mappingSet was generated from
	//unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int ii=0; ii<n_reactants; ii++)
	{
		reactantLists[ii]->pickRandom(mappingSet[ii]);
	}
}

double Compartment::update_a()
{
	for(unsigned int i=0; i<n_reactants; i++)
		a*=reactantLists[i]->size();
	a*=m_pParentReaction->baseRate;
	return a;
}

void Compartment::init(CompartmentReaction* thisReactionClass)
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		thisReactionClass->reactantTemplates[r]->getMoleculeType()->addReactionClass(thisReactionClass,r);
	}
}

inline bool Compartment::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//Get the specified reactantList
	rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(m_pParentReaction,reactantPos);

	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		if(!m_pParentReaction->reactantTemplates[reactantPos]->compare(m)) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}

	} else {
		//Try to map it!
		if(!m_pParentReaction->reactantTemplates[reactantPos]->compare(m,rl,ms)) {
			rl->popLastMappingSet();
			//we just pushed, then popped, so we a has not changed...
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}
	return true;
}
