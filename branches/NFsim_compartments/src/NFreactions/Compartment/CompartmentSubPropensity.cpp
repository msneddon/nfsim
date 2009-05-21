/*
 * CompartmentInteraction.cpp
 *
 *  Created on: May 18, 2009
 *      Author: Kelly Stanton
 */

#include "CompartmentSubPropensity.h"

CompartmentSubPropensity::CompartmentSubPropensity(
		CompartmentReaction* parentCompartmentReaction,
		unsigned int n_reactants,
		unsigned int CompartmentIdList[]) {

	for(unsigned int ii = 0; ii < n_reactants; ii++)
	{
		if(CompartmentIdList[ii]==(unsigned int)-1)
		{
			throw "CompartmentSubPropensity instantiated with bad CompartmentIdList:  Cannot contain (unsigned int)-1";
		}
	}

		// take the list of compartments in this interaction
	this->m_parentCompartmentReaction = parentCompartmentReaction;
	this->n_reactants = n_reactants;
	this->m_CompartmentIdList = CompartmentIdList;

}

CompartmentSubPropensity::~CompartmentSubPropensity() {

}

double CompartmentSubPropensity::update_a()
{

	a=1;
	// ii is both the position of the compartment in the ID list and the
	// reactant in the reactant list
	for(unsigned int ii=0; ii<n_reactants; ii++)
	{
		//get the compartment
		CompartmentReactantList* tmp = m_parentCompartmentReaction->m_mapCompartmentList[m_CompartmentIdList[ii]];
		//get the reactantlist
		a*=tmp->reactantLists[ii]->size();
	}
	a*=m_parentCompartmentReaction->baseRate;
	return a;
}
void CompartmentSubPropensity::pickMappingSets(MappingSet** mappingSet, double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//So, we shall loop through the lists and extract out the MappingSets and
	//the molecule ID that the mappingSet was generated from
	//unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int ii=0; ii<n_reactants; ii++)
	{
		// get the compartment that this list is in
		CompartmentReactantList* tmp = m_parentCompartmentReaction->m_mapCompartmentList[m_CompartmentIdList[ii]];
		// pick the mapping set
		tmp->reactantLists[ii]->pickRandom(mappingSet[ii]);
	}
}
