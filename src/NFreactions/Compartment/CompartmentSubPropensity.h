/*
 * CompartmentInteraction.h
 *
 *  Created on: May 18, 2009
 *      Author: Kelly Stanton
 *
 *      This Class represents the interaction between one or two compartments for
 *      a compartment reaction class.  These interaction compartnents can
 *      all contribute reactants to the reaction so the propensities need
 *      to consider this.
 */

#ifndef COMPARTMENTSUBPROPENSITY_H_
#define COMPARTMENTSUBPROPENSITY_H_

#include <map>
#include <vector>
#include "CompartmentReaction.h"
#include "Compartment.h"
#include "../NFreactions.hh"

using namespace NFcore;
using namespace std;

class CompartmentReaction;
class CompartmentReactantList;
class CompartmentSubPropensity{
	friend class CompartmentReactantList;
	friend class CompartmentReaction;
public:
	CompartmentSubPropensity(
			CompartmentReaction* parentCompartmentReaction,
			unsigned int n_reactants,
			unsigned int CompartmentIdList[]);
	virtual ~CompartmentSubPropensity();

	// pick a mapping set from the reactant lists located in the two compartments this interaction spans
	void pickMappingSets(MappingSet** mappingSet, double random_A_number) const;
	//Get and update the propensity
	double update_a();

	unsigned int* m_CompartmentIdList;
	CompartmentReaction* m_parentCompartmentReaction;
	unsigned int n_reactants;
	//propensity
	double a;
};

#endif /* COMPARTMENTSUBPROPENSITY_H_ */
