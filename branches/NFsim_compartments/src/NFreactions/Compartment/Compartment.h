/*
 * Compartment.h
 *
 *  Created on: Apr 16, 2009
 *      Author: kstanton
 *
 *      The compartment class represents a compartment that an individual reaction
 *      can occur in.  The compartment reaction contains a list of these objects
 *      and will choose reactions to occur based on the number of molecules present
 *      in the same compartment
 */

#ifndef COMPARTMENT_H_
#define COMPARTMENT_H_

#include "../NFreactions.hh"
#include "CompartmentReaction.h"
#include<string>

using namespace NFcore;
using namespace std;
class CompartmentReaction;
class Compartment {
	friend class CompartmentReaction;
public:
	//Constructor
	Compartment(
			unsigned int compId,
			double dbBaseRate,
			TransformationSet* transformationSet,
			unsigned int nReactants,
			CompartmentReaction*);
	//Destructor
	virtual ~Compartment();
	//initialization
	void init(CompartmentReaction* thisReactionClass);
	//Gets the number of reactants
	unsigned int getReactantCount(unsigned int reactantIndex) const;
	//Selects a random mapping set
	void pickMappingSets(MappingSet** mappingSet, double random_A_number) const;
	//Prints details of the reaction mainly used for debugging
	void printDetails() const;
	//Get and update the propensity
	double update_a();
	//Try to add the molecule to this reaction rule
	virtual bool tryToAdd(Molecule *m, unsigned int reactantPos, int rxnIndex);
	//Removese a molecule from this reaction rule
	virtual void tryToRemove(Molecule *m, unsigned int reactantPos, int rxnIndex);
	//remove clone mappings
	void removeClones(Molecule *m, unsigned int reactantPos, int rxnIndex);

protected:
	// number of reactants in this compartment
	unsigned int n_reactants;
	// name of this compartment
	string name;
	//List of reactant lists
	ReactantList **reactantLists;
	//Reactant List of all agents that satisfy the template for one species
	//required to fire this reaction
	ReactantList *rl;
	//mapping set
	MappingSet *ms;
	//Reaction this compartment belongs to
	CompartmentReaction* m_pParentReaction;
	//The compartment id
	unsigned int m_iCompId;
	//propensity
	double a;
	//Whether this compartment is used or not
	bool active;
	//Volume of this compartment
	double volume;
};

#endif /* COMPARTMENT_H_ */
