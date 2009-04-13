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
	virtual double update_a() {
			a = 1;

			for(unsigned int i=0; i<n_reactants; i++)
				a*=reactantLists[i]->size();

			a*=baseRate;
			return a;
	}


	virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
	virtual unsigned int getReactantCount(unsigned int reactantIndex) const;

	virtual void printFullDetails() const;
	bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos);

protected:
	virtual void pickMappingSets(double randNumber) const;

	ReactantList **reactantLists;

	ReactantList *rl;
	MappingSet *ms;
};

#endif /* COMPARTMENTREACTION_H_ */




