/*
 * CompartmentInteraction.cpp
 *
 *  Created on: May 18, 2009
 *      Author: Kelly Stanton
 */

#include "CompartmentInteraction.h"

CompartmentInteraction::CompartmentInteraction(
		CompartmentReaction* m_parentCompartmentReaction) {

		// take the list of compartments in this interaction
		static vector<unsigned int>::iterator compIter;
		/*
		for(compIter = CompartmentIdList.begin(); compIter != CompartmentIdList.end(); compIter++)
		{
			// add them to the compartments subset that this interaction takes care of
			unsigned int tmpCompId = *compIter;
			m_mapInteractionSubsetCompartmentsList.insert(
							pair<unsigned int, Compartment*>(
									tmpCompId,
									m_parentCompartmentReaction->m_mapCompartmentList[tmpCompId]));
		}
		*/
}

CompartmentInteraction::~CompartmentInteraction() {

}

double CompartmentInteraction::update_a()
{
	/*
	// We return 0 propensity for inactive compartments to fire
	if(!active) return 0;
	a=1;
	for(unsigned int i=0; i<n_reactants; i++)
		a*=reactantLists[i]->size();
	a*=m_pParentReaction->baseRate;
	return a;
	*/
}
