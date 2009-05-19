/*
 * CompartmentInteraction.h
 *
 *  Created on: May 18, 2009
 *      Author: Kelly Stanton
 *
 *      This Class represents the interaction between two compartments for
 *      a compartment reaction class.  These interaction compartmnents can
 *      all contribute reactants to the reaction so the propensities need
 *      to consider this.
 */

#ifndef COMPARTMENTINTERACTION_H_
#define COMPARTMENTINTERACTION_H_

#include <map>
#include <vector>
#include "CompartmentReaction.h"
#include "Compartment.h"

class CompartmentReaction;
class Compartment;
class CompartmentInteraction {
	friend class Compartment;
	friend class CompartmentReaction;
public:
	CompartmentInteraction(
			CompartmentReaction* m_parentCompartmentReaction);
	virtual ~CompartmentInteraction();

	//Get and update the propensity
	double update_a();

	map <unsigned int, Compartment*> m_mapInteractionSubsetCompartmentsList;
	CompartmentReaction* m_parentCompartmentReaction;
};

#endif /* COMPARTMENTINTERACTION_H_ */
