


#include "reaction.hh"

#define DEBUG_MESSAGE 0


using namespace std;
using namespace NFcore;


//should also accept list of local functions and list of PointerNames for each of the functions...
PopulationRxnClass::PopulationRxnClass(
		string name,
		double baseRate,
		string baseRateName,
		TransformationSet *transformationSet,
		System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	//Remember that we are a Population ReactionClass
	this->reactionType = ReactionClass::POP_RXN;

	//Set up the reactantTrees
	reactantTrees = new ReactantTrees *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		reactantTrees[r]=(new ReactantTree(r,transformationSet,32));
	}

	//Initialize a to zero
	this->a=0;
}

PopulationRxnClass::~PopulationRxnClass() {

	if(DEBUG) cout<<"Destroying rxn: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++) {
			delete reactantTrees[r];
	}

	delete [] reactantTrees;
}

void PopulationRxnClass::init() {

	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}



void PopulationRxnClass::remove(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking...
	if( reactantPos<0 || reactantPos>=n_reactants || m==NULL )
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	if(m->getRxnListMappingId(rxnIndex)>=0) {
		// was in the tree, so we should remove
		reactantTrees[reactantPos]->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}


bool PopulationRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {

	// Get the specified reactantTree
	rt = reactantTrees[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	//If this reaction has multiple instances, we always remove them all!
	// then we remap because other mappings may have changed.  Yes, this may
	// be more ineffecient, but it is the fast implementation
	if( rt->getHasClonedMappings() ) {
		if( m->getRxnListMappingId(rxnIndex)>=0 ) {
			rt->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex, Molecule::NOT_IN_RXN);
		}
	}

	// Here we get the standard update...
	if( m->getRxnListMappingId(rxnIndex)>=0 ) {
		// was in the tree, so checking if we should remove"
		if(!reactantTemplates[reactantPos]->compare(m)) {
			//if(DEBUG_MESSAGE)cout<<"removing..."<<endl;
			rt->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex, Molecule::NOT_IN_RXN);
		} else {}
	} else {
		// wasn't in the tree, so trying to push and compare
		ms = rt->pushNextAvailableMappingSet();
		if(!reactantTemplates[reactantPos]->compare(m,rt,ms)) {
			//cout<<"shouldn't be in the tree, so we pop"<<endl;
			rt->removeMappingSet(ms->getId());
		} else {
			// should be in the tree, so confirm push
			rt->confirmPush(ms->getId(), ms->getComplexPopulation());
			m->setRxnListMappingId(rxnIndex, ms->getId());
		}
	}
	return true;
}


int PopulationRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return reactantTrees[reactantIndex]->getRateFactorSum();
}


int PopulationRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	return reactantTrees[reactantIndex]->getRateFactorSum();
}


double PopulationRxnClass::update_a() {
	a = baseRate;
	for(unsigned int i=0; i<n_reactants; i++) {
		a*=reactantTrees[i]->getRateFactorSum();
	}
	return a;
}

void PopulationRxnClass::pickMappingSets(double randNumber) const
{
	// we need a random number for each reactant tree, so we will generate our own,
	//  rather than using the provided number
	double randnum;
	for(unsigned int i=0; i<n_reactants; i++) {
		randnum = NFutil::RANDOM(reactantTrees[i]->getRateFactorSum());
		reactantTrees[i]->pickReactantFromValue(mappingSet[i],randnum,1.0);
	}
}

void PopulationRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) {
	cerr<<"You are trying to notify a Population Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}


void PopulationRxnClass::printDetails() const
{
	cout<<"PopulationRxnClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout << "      -|"<< this->reactantTrees[r]->size() << " mappings|\t";
		cout << this->reactantTemplates[r]->getPatternString() << "\n";
		cout << "             (rateFactorSum="<<reactantTrees[r]->getRateFactorSum() << ")."<<endl;
	}

	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}





