/*
 * CompartmentReaction.h
 *
 *  Created on: Apr 13, 2009
 *      Author: kstanton
 *
 *      The CompartmentReaction class specifies an extended version of the ReactionClass
 *      it adds functionality for storing molecules and the molecule lists in a list
 *      of compartments.  This is still an abstract class and will need to be derived
 *      for specific reaction types
 */


#ifndef COMPARTMENTREACTION_H_
#define COMPARTMENTREACTION_H_

#include "../NFreactions.hh"

#include "Compartment.h"
#include <vector>
using namespace NFcore;
using namespace std;

class Compartment;
class CompartmentReaction: public NFcore::ReactionClass {
	friend class Compartment;
public:
	//Constructor
	CompartmentReaction(string name, double baseRate, TransformationSet *transformationSet);
	//Destructor
	virtual ~CompartmentReaction();
	//Initilization
	virtual void init();
	//Finalization
	virtual void prepareForSimulation();
	//Adds a molecule to this reaction rule
	virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
	//Removese a molecule from this reaction rule
	virtual void remove(Molecule *m, unsigned int reactantPos);
	//Gets the propensity for this reaction rule to fire
	virtual double update_a();

	virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
	//Gets the number of reactants
	virtual unsigned int getReactantCount(unsigned int reactantIndex) const;
	//prints information about this reaction.  Mainly used for debugging
	virtual void printFullDetails() const;

protected:
	//Choose a mappingSet at random
	virtual void pickMappingSets(double randNumber) const;

	//note maps are access in log(n) time
	map<unsigned int,Compartment*> m_mapCompartmentList;
};

#endif /* COMPARTMENTREACTION_H_ */




