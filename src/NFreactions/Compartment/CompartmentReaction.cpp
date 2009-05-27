/*
 * CompartmentReaction.cpp
 *
 *  Created on: Apr 13, 2009
 *      Author: kstanton
 */

#include "CompartmentReaction.h"
#include <algorithm> // check stl_algo.h for details

unsigned int CompartmentReaction::nCompartments = 0;
map<unsigned int,vector<unsigned int> > CompartmentReaction::m_mapCompartmentConnectivity;
CompartmentReaction::CompartmentReaction(string name, double baseRate, TransformationSet *transformationSet) :
	ReactionClass(name,baseRate,transformationSet)
{
	this->reactionType = COMPARTMENT_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	if(nCompartments <1)
	{
		throw "Number of compartments not initialized";
	}

	// create the correct number of compartments
	for(unsigned int ii=0; ii<nCompartments; ii++)
	{
		// create the compartment
		m_mapCompartmentList.insert(
				pair<unsigned int, CompartmentReactantList*>(
						ii,
						new CompartmentReactantList(ii,baseRate,transformationSet,n_reactants,this)));
	}

	// if this is a unimolecular reaction
	if(n_reactants ==1)
	{
		// find the compartment its in
		unsigned int cpt1 = reactantTemplates[0]->getCompartmentConstraint();

		// if the reactant is universal
		if(cpt1 == (unsigned int)-1)
		{
			// put the interaction in all compartments
			static map <unsigned int, CompartmentReactantList*>::iterator cmpIter;
			for(cmpIter =m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
			{
				this->addCompartmentInteraction((unsigned int[]){cmpIter->first});
			}
		}
		// just put it in the compartment it is constrained to
		else
		{
			this->addCompartmentInteraction((unsigned int[]){cpt1});
		}
	}
	//if we have exactly 2 reactants set up the interactions
	if(n_reactants ==2)
	{
		unsigned int cpt1 = reactantTemplates[0]->getCompartmentConstraint();
		unsigned int cpt2 = reactantTemplates[1]->getCompartmentConstraint();

		// if both reactants are universal
		if(cpt1 == (unsigned int)-1 && cpt2 == (unsigned int)-1)
		{
			// add all combinations of adjacent compartments
			static map<unsigned int,vector<unsigned int> >::iterator cmpVectIter;
			for(cmpVectIter = m_mapCompartmentConnectivity.begin(); cmpVectIter != m_mapCompartmentConnectivity.end(); cmpVectIter++)
			{
				static vector<unsigned int>::iterator cmpConnectedVectIter;
				for(cmpConnectedVectIter = cmpVectIter->second.begin(); cmpConnectedVectIter != cmpVectIter->second.end(); cmpConnectedVectIter++)
				{
					this->addCompartmentInteraction((unsigned int[]){cmpVectIter->first,*cmpConnectedVectIter});
				}
			}
			// Add interactions drawing reactants from the same compartments
			static map<unsigned int,CompartmentReactantList*>::iterator cmpIter;
			for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
			{
				this->addCompartmentInteraction((unsigned int[]){cmpIter->first,cmpIter->first});
			}

		}
		// if the first reactant is universal
		else if(cpt1 == (unsigned int)-1)
		{
			// all compartments for 1 that are connected to 2
			static vector<unsigned int>::iterator cmpConnectedVectIter;
			for(cmpConnectedVectIter =m_mapCompartmentConnectivity[cpt2].begin(); cmpConnectedVectIter != m_mapCompartmentConnectivity[cpt2].end(); cmpConnectedVectIter++)
			{
				this->addCompartmentInteraction((unsigned int[]){*cmpConnectedVectIter,cpt2});
			}
			// add interaction for drawing reactants from the same compartment
			this->addCompartmentInteraction((unsigned int[]){cpt2,cpt2});

		}
		// if the second reactant is universal
		else if(cpt2 == (unsigned int)-1)
		{
			// add compartments connected to 1
			static vector<unsigned int>::iterator cmpConnectedVectIter;
			for(cmpConnectedVectIter =m_mapCompartmentConnectivity[cpt1].begin(); cmpConnectedVectIter != m_mapCompartmentConnectivity[cpt1].end(); cmpConnectedVectIter++)
			{
				this->addCompartmentInteraction((unsigned int[]){cpt1,*cmpConnectedVectIter});
			}
			// add interaction for drawing reactants from the same compartment
			this->addCompartmentInteraction((unsigned int[]){cpt1,cpt1});
		}
		// if neither reactant is universal
		else
		{
			// if the compartments are adjacent
			vector <unsigned int> a = m_mapCompartmentConnectivity[cpt1];
			if(find(a.begin(), a.end(), cpt2) != a.end())
			{
				this->addCompartmentInteraction((unsigned int[]){cpt1,cpt2});
			}
		}
	}
	else if(n_reactants > 2)
	{
		throw "Compartment reactions do not currently support reactions with 3 or more reactants.  Update Compartment interaction initialization to add this.";
	}

	/* the general case is very complex and would only occur if we have more than 2
	 * template molecules
	 */

}

CompartmentReaction::~CompartmentReaction()
{
	if(DEBUG) cout<<"Destorying rxn: "<<name<<endl;
	static map<unsigned int,CompartmentReactantList*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		delete cmpIter->second;
	}
	static vector<CompartmentSubPropensity*>::iterator interactionIter;
	for(interactionIter = m_vectInteractionList.begin(); interactionIter != m_vectInteractionList.end(); interactionIter++)
	{
		delete *interactionIter;
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
	// gather a little background on this molecule and reaction pair

	//look up the index
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//Find the reaction ID of this molecule in this reaction if it is in this reaction at all if not -1
	int rxnId = m->getRxnListMappingId(rxnIndex);
	//Find the compartment that this Molecule is currently in in this reaction  if not -1
	int oldCompartmentId = m->getRxnListCompartmentMappingId(rxnIndex);
	// New Compartment is taken from the transformed molecule
	int newCompartmentId = m->getCompartmentId();

	// first off this is a compartment reaction and must have a compartment ID if its in the molecule
	// if its not in a compartment
	if(oldCompartmentId < 0)
	{
		if(rxnId<0) // and its not in this reaction
		{
			// then lets initialize it to the compartment that it is located in
			oldCompartmentId = newCompartmentId;
		}
		else // oops no compartment but its in this reaction
		{
			throw "No compartment RxnListCompartmentMappingId associated with this compartment reaction for which there is a RxnListMappingId";
		}
	}
	// if either the source or target compartment do not exist
	if(m_mapCompartmentList.find(oldCompartmentId) == m_mapCompartmentList.end()
			||m_mapCompartmentList.find(newCompartmentId) == m_mapCompartmentList.end()
			|| m_mapCompartmentList[oldCompartmentId] == NULL
			|| m_mapCompartmentList[newCompartmentId] == NULL)
	{
		throw "Dereferencing a null compartment!";
	}

	// remove clones first
	m_mapCompartmentList[oldCompartmentId]->removeClones(m,reactantPos,rxnIndex);
	m_mapCompartmentList[m->getCompartmentId()]->removeClones(m,reactantPos,rxnIndex);

	//Here we get the standard update...
	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		m_mapCompartmentList[oldCompartmentId]->tryToRemove(m, reactantPos, rxnIndex);
	} else {
		m_mapCompartmentList[m->getCompartmentId()]->tryToAdd(m, reactantPos, rxnIndex);
	}
	return true;
}


double CompartmentReaction::update_a()
{
	// we no longer check propensity of the compartments alone
	// the CompartmentSubPropensity takes care of this now
	/*
	static map<unsigned int,CompartmentReactantList*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		a+=cmpIter->second->update_a();
	}
	*/

	static vector<CompartmentSubPropensity*>::iterator subPropIter;
	a = 0;
	for(subPropIter = m_vectInteractionList.begin(); subPropIter != m_vectInteractionList.end(); subPropIter++)
	{
		a+=(*subPropIter)->update_a();
	}
	return a;
}


void CompartmentReaction::init()
{
	static map<unsigned int,CompartmentReactantList*>::iterator cmpIter;
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

	// get the reaction index
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	// check to make sure this molecule is in this reaction and that the compartment is also filled out
	if((m->getRxnListMappingId(rxnIndex)>=0) && (m->getRxnListCompartmentMappingId(rxnIndex)>=0))
	{
		m_mapCompartmentList[m->getRxnListCompartmentMappingId(rxnIndex)]->tryToRemove(m, reactantPos, rxnIndex);
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
	static map<unsigned int,CompartmentReactantList*>::const_iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		nReactants = nReactants + cmpIter->second->getReactantCount(reactantIndex);
	}
	return nReactants;
}
void CompartmentReaction::printFullDetails() const
{
	cout<<"BasicRxnClass: "<<name<<endl;
	static map<unsigned int,CompartmentReactantList*>::const_iterator cmpIter;
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
	//CompartmentReactantList* chosenComp = 0;
	double leftover_random_A_number = 0;

	//WARNING - DO NOT USE THE DEFAULT C++ RANDOM NUMBER GENERATOR FOR THIS STEP
	// - IT INTRODUCES SMALL NUMERICAL ERRORS CAUSING THE ORDER OF RXNS TO
	//   AFFECT SIMULATION RESULTS
	// we now choose a subpropensity instead of a compartment
	/*
	static map<unsigned int,CompartmentReactantList*>::const_iterator cmpIter;
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
	// if we chose a compartment get the mapping set from it
	if(chosenComp!=0)
	{
		chosenComp->pickMappingSets(mappingSet, leftover_random_A_number);
	}
	*/

	CompartmentSubPropensity* chosenSubProp = 0;
	// or choose an interaction
	static vector<CompartmentSubPropensity*>::const_iterator interactIter;
	for(interactIter = m_vectInteractionList.begin(); interactIter != m_vectInteractionList.end(); interactIter++)
	{
		a_sum += (*interactIter)->a;
		if (random_A_number <= a_sum && chosenSubProp==0)
		{
			chosenSubProp = *interactIter;
			leftover_random_A_number = random_A_number-last_a_sum;
			break;
		}
		last_a_sum = a_sum;
	}
	if(chosenSubProp==0)
	{
		cerr<<"Error: random_A_number exceeds a_sum!!!"<<endl;
		cerr<<"randNum: "<<random_A_number<<"  a_sum: "<< a_sum<<" running a_tot:"<<a<<endl;
	}
	chosenSubProp->pickMappingSets(mappingSet, leftover_random_A_number);

}

void CompartmentReaction::restrictToCompartment(unsigned int compartmentId)
{
	static map<unsigned int,CompartmentReactantList*>::iterator cmpIter;
	for(cmpIter = m_mapCompartmentList.begin(); cmpIter != m_mapCompartmentList.end(); cmpIter++)
	{
		cmpIter->second->active = false;
	}
	m_mapCompartmentList[compartmentId]->active = true;
}
void CompartmentReaction::addCompartmentInteraction(unsigned int CompartmentIdList[])
{
	// not the position of the compartment in the compartmentIdList corresponds to the reactant pos in the reactantList
	CompartmentSubPropensity* tmp = new CompartmentSubPropensity(this, n_reactants, CompartmentIdList);
	m_vectInteractionList.push_back(tmp);
}
void CompartmentReaction::SetNumCompartments(unsigned int nCompartments)
{
	CompartmentReaction::nCompartments = nCompartments;
}
void CompartmentReaction::addConnectivity(unsigned int compartmentId1, unsigned int compartmentId2)
{
	m_mapCompartmentConnectivity[compartmentId1].push_back(compartmentId2);
	m_mapCompartmentConnectivity[compartmentId2].push_back(compartmentId1);
}
