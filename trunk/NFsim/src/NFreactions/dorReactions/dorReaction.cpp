

#include "dorReaction.hh"




using namespace NFcore;
using namespace std;



DORrxnClass::DORrxnClass(
	string name, 
	double baseRate, 
	TransformationSet *transformationSet, 
	unsigned int DORreactantIndex, 
	string DORgroupName, 
	int DORgroupValueIndex) :
	ReactionClass(name,baseRate,transformationSet)
{
	this->reactionType = DOR_RXN;
	this->DORreactantIndex = DORreactantIndex;
	this->DORgroupName = DORgroupName;
	this->DORgroupValueIndex = DORgroupValueIndex;
	this->rateFactorSum = new double [n_reactants];
	
	reactantTree = new ReactantTree(DORreactantIndex, transformationSet, 4);
	
	reactantLists = new ReactantList *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if(r==DORreactantIndex) {
			reactantLists[r]=NULL;
		} else {
			reactantLists[r]=new ReactantList(r,transformationSet,25);
		}
	}
	
	for(unsigned int s=0; s<n_reactants; s++) rateFactorSum[s] = 0;
}
	
DORrxnClass::~DORrxnClass()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if((int)r!=DORreactantIndex) delete reactantLists[r];
		delete [] rateFactors[r];
	}
	
	delete [] rateFactors;
	delete [] rateFactorSum;
	delete reactantTree;
	delete [] reactantLists;
}
				
				
void DORrxnClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++) {
		if((int)r==DORreactantIndex) {
			reactantTemplates[r]->getMoleculeType()->addDORrxnClass(this,r);
		} else {
			reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
		}
	}
}


void DORrxnClass::prepareForSimulation() 
{
	
}


bool DORrxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL) 
	{
		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}
	
	
	//Get the specified reactantList
	//ReactantList *rl = reactantLists.at(position);
	
	//Check if the molecule is in the list or tree at this position
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);
	
	
	if(isInRxn)
	{
		if(!reactantTemplates[reactantPos]->compare(m)) {
			if(DORreactantIndex==(int)reactantPos) {
				reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			} else {
				reactantLists[reactantPos]->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			}
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
		
	} else {
		//Try to map it!
		MappingSet *ms;
		if(DORreactantIndex==(int)reactantPos) {
			double rFactor = m->getDORvalueFromGroup(DORgroupName, DORgroupValueIndex);
			ms = reactantTree->pushNextAvailableMappingSet(rFactor);
		} else {
			ms = reactantLists[reactantPos]->pushNextAvailableMappingSet();
		}
		if(!reactantTemplates[reactantPos]->compare(m,ms)) {
			if(DORreactantIndex==(int)reactantPos) {
				reactantTree->removeMappingSet(ms->getId());
			} else {
				reactantLists[reactantPos]->popLastMappingSet();
			}
		}
		m->setRxnListMappingId(rxnIndex,ms->getId());
	}
		
	return true;
}


void DORrxnClass::remove(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL) 
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}
		
		
	//Get the specified reactantList
	//ReactantList *rl = reactantLists.at(position);
		
	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);
		
		
	if(isInRxn)
	{
		if(DORreactantIndex==(int)reactantPos) {
			reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		} else {
			reactantLists[reactantPos]->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		}
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}	
}
				
double DORrxnClass::update_a()
{
	bool CHECK_PRECISION = false;
	double MAX_ERROR = 0.0001;
	
	//When we calculate the new reaction propensity, it is just 
	a = baseRate;  //Start at 1 because rate factor sum already has rate multiplied in
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if((int)r == DORreactantIndex)
		{
			//Here we can check the precision to make sure we don't have any
			//numerical errors from adding and subtracting from the rateFactorSum
			if(CHECK_PRECISION)
			{
				double cumRateFactorSum = 0;
//			for(int m=0; m<lastIndex[r]+1; m++) cumRateFactorSum +=  rateFactors[r][m];
//				if(fabs(cumRateFactorSum-rateFactorSum[r]) >= MAX_ERROR)
				{
					cerr<<"In DOR RXN: "<<name<<" NUMERICAL ERROR IN RATE FACTOR >= " << MAX_ERROR << ".  CORRECTING. "<<endl;
					rateFactorSum[r] = cumRateFactorSum;
				}
			}
			//cout<<" dor Rate Factor Sum: " << dorReactants->getRateFactorSum() << endl;
			double rfSum = reactantTree->getRateFactorSum();
			if(reactantTree->size()==0 || rfSum <=0 ) {a=0; return a; }
			a*= rfSum;
		}
		else
		{
			a*= reactantLists[r]->size();
		}
	}
	return a;
	
}
				
				
void DORrxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	//First make sure m is in fact a molecule here, if it is not, return
	if(rxnListIndex<0) return;
	if(reactantIndex!=DORreactantIndex) {
		cerr<<"Trying to update a molecule that is not a DOR reactant!!"<<endl;
		return;
	}
	
	//Get the new rFactor for this molecule
	double rFactor = m->getDORvalueFromGroup(DORgroupName, DORgroupValueIndex);

	//Add in this rFactor to the rFactorSum
	//rateFactorSum[reactantIndex] -= rateFactors[reactantIndex][rxnListIndex];
	//rateFactorSum[reactantIndex] += rFactor;
	
	//Remember the new rFactor
	//rateFactors[reactantIndex][rxnListIndex] = rFactor;
	
	reactantTree->updateValue(rxnListIndex, rFactor);
	
}
				
				
				
unsigned int DORrxnClass::getReactantCount(unsigned int reactantIndex) const
{
	if(DORreactantIndex==(int)reactantIndex) {
		reactantTree->size();
	} else {
		reactantLists[reactantIndex]->size();
	}
}

void DORrxnClass::printFullDetails() const
{
	cout<<"** DORReaction: " << name <<"  ( rate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if((int)r==DORreactantIndex)
		{
			cout<<"      -(DOR) "<< this->reactantTemplates[r]->getMoleculeTypeName();
			cout<<"	(count="<< reactantTree->size() <<", RF Sum="<<reactantTree->getRateFactorSum()<<")."<<endl;		
		}
		else
		{
			cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
			cout<<"	(count="<< reactantLists[r]->size() <<")."<<endl;
		}
	}
}



void DORrxnClass::pickMappingSets(double random_A_number) const
{
	
	//First, loop through all the non-DOR reactants to figure out the base rate factor.
	//This 'total' base rate factor is needed because the random number given will be
	//between 0 and the total propensity of the entire reaction.  To get that total propensity,
	//we have to add up the propensity based on the other rates, and then (and ONLY then)
	//can we pick the reactant from value from the tree using that value.
	double baseRateFactor = this->baseRate;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if((int)r == DORreactantIndex) {
			continue;
		} else {
			reactantLists[r]->pickRandom(mappingSet[r]);
			baseRateFactor*=reactantLists[r]->size();
		}
	}
	reactantTree->pickReactantFromValue(mappingSet[DORreactantIndex],random_A_number,baseRateFactor);
}

