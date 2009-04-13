/*
 * CompartmentReaction.h
 *
 *  Created on: Apr 13, 2009
 *      Author: kstanton
 */

#ifndef COMPARTMENTREACTION_H_
#define COMPARTMENTREACTION_H_

#include "../NFreactions.hh"

class CompartmentReaction: public NFcore::ReactionClass {

public:
	CompartmentReaction(string name, double baseRate, TransformationSet *transformationSet);
	virtual ~CompartmentReaction();
	virtual void init();
	virtual void prepareForSimulation();
	virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
	virtual void remove(Molecule *m, unsigned int reactantPos);
	virtual double update_a();
	virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
	virtual unsigned int getReactantCount(unsigned int reactantIndex) const;

	virtual void printFullDetails() const;
	bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos);

protected:
	virtual void pickMappingSets(double randNumber) const;

	CompartmentList *compartmentLists;
};

class Compartment{
public:
	ReactantList **reactantLists;
	ReactantList *rl;
	MappingSet *ms;
};

#endif /* COMPARTMENTREACTION_H_ */




