

#include "NFrateChangeReactions.hh"
#include <math.h>

#define CHECK_PRECISION 0
#define MAX_ERROR 0.00000001

using namespace NFcore;

////////////////////////////////////////////////////////////////////
///////////////  Generalized DOR reaction

DORrxn::DORrxn( char * name, 
				int n_reactants, 
				TemplateMolecule ** reactantTemplates,
				int DORreactantIndex,
				char * DORgroupName,
				int DORgroupValueIndex,
				double baseRate) : ReactionClass(name, n_reactants, reactantTemplates, baseRate)
{
	
	if(DEBUG) cout<<"Creating generalized DOR reaction..."<<endl;	
	this->reactionType = DOR_RXN;
	this->DORreactantIndex = DORreactantIndex;
	this->DORgroupName = DORgroupName;
	this->DORgroupValueIndex = DORgroupValueIndex;
	this->rateFactorSum = new double [n_reactants];
	
	for(int s=0; s<n_reactants; s++) rateFactorSum[s] = 0;
}


DORrxn::~DORrxn()
{
	for(int r=0; r<n_reactants; r++)
		delete [] rateFactors[r];
	
	delete [] rateFactors;
	delete [] rateFactorSum;
	delete dorReactants;
}


void DORrxn::init()
{ 
	if(DEBUG) cout<<"Init DOR rxn"<<endl;
	for(int r=0; r<n_reactants; r++)
	{
		if(r==DORreactantIndex)
		{
			dorReactants = new NFReactantTree(reactantTemplates[r]->getMoleculeType()->getMoleculeCount());
			reactantTemplates[r]->getMoleculeType()->addDORrxnClass(this,r);
		}
		else
			reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}

void DORrxn::prepareForSimulation()
{
	reactantArray = new Molecule ** [n_reactants];
	maxSize = new int[n_reactants];
	lastIndex = new int[n_reactants];
	rxnIndexInMoleculeType = new int[n_reactants];
	
	rateFactors = new double * [n_reactants];
	
	for(int r=0; r<n_reactants; r++)
	{
		maxSize[r] = reactantTemplates[r]->getMoleculeType()->getMoleculeCount();
		lastIndex[r] = -1;
		reactantArray[r] = new Molecule * [maxSize[r]];
		rateFactors[r] = new double[maxSize[r]];
		for(int k=0; k<maxSize[r]; k++)
			rateFactors[r][k]=0;
		rxnIndexInMoleculeType[r] = reactantTemplates[r]->getMoleculeType()->getRxnIndex(this,r);
	}
	
	prepared = true;
}

int DORrxn::addToReactantList(Molecule *m, int reactantIndex)
{
	
	//If we added a DOR reactant, then we also have to add a rate factor
	//at the cooresponding position in the array
	if(DORreactantIndex==reactantIndex)
	{
		double rFactor = m->getDORvalueFromGroup(DORgroupName, DORgroupValueIndex);
		//rateFactors[reactantIndex][lastIndex[reactantIndex]] = rFactor;
		//rateFactorSum[reactantIndex] += rFactor;
		
		//Add it to the tree
		return dorReactants->insert(m, rFactor);
	}
	
	//Always add to the end of the array (which works even if array is empty)
	reactantArray[reactantIndex][++lastIndex[reactantIndex]] = m;
	
	
	return lastIndex[reactantIndex];;
}

void DORrxn::removeFromReactantList(Molecule *m, int reactantIndex, int rxnListIndex)
{
	//If we want to remove the last reactant, no swapping is necessary
	if(lastIndex[reactantIndex]== rxnListIndex)
	{
		
		
		//If this is a DOR reactant, we have to subtract from the rateFactorSum
		if(reactantIndex == DORreactantIndex)
		{
			//rateFactorSum[reactantIndex] -= rateFactors[reactantIndex][rxnListIndex];
			//if(0==rxnListIndex)
			//	cout<<"here "<< rateFactorSum[reactantIndex] << endl;
			dorReactants->remove(m,rxnListIndex);
			return;
		}
		//cout<<"here in last index?"<<endl;
		
		reactantArray[reactantIndex][rxnListIndex] = 0;
		lastIndex[reactantIndex]--;
		
		return;
	}
	
	
	//And make sure if it was a DOR reactant, we update the rateFactorSum
	if(reactantIndex == DORreactantIndex)
	{
		//First decrement the rate factor sum with the value we last had
		//rateFactorSum[reactantIndex] -= rateFactors[reactantIndex][rxnListIndex];
		
		//then swap the positions in the array
		//rateFactors[reactantIndex][rxnListIndex] = rateFactors[reactantIndex][lastIndex[reactantIndex]];
		//rateFactors[reactantIndex][lastIndex[reactantIndex]] = 0;
		dorReactants->remove(m,rxnListIndex);
		return;
	}
	
	
	//If we get here, then we're going to have to swap some positions
	reactantArray[reactantIndex][rxnListIndex] = reactantArray[reactantIndex][lastIndex[reactantIndex]];
	reactantArray[reactantIndex][lastIndex[reactantIndex]]=0;
	reactantArray[reactantIndex][rxnListIndex]->setRxnListIndex(rxnIndexInMoleculeType[reactantIndex], rxnListIndex); 
	
	
	//Finally, remember that we need to decrease the last index
	lastIndex[reactantIndex]--;
	
	return;
}

void DORrxn::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	//First make sure m is in fact a molecule here, if it is not, return
	if(rxnListIndex==Molecule::NOT_IN_RXN) return;
	
	//Get the new rFactor for this molecule
	double rFactor = m->getDORvalueFromGroup(DORgroupName, DORgroupValueIndex);

	//Add in this rFactor to the rFactorSum
	//rateFactorSum[reactantIndex] -= rateFactors[reactantIndex][rxnListIndex];
	//rateFactorSum[reactantIndex] += rFactor;
	
	//Remember the new rFactor
	//rateFactors[reactantIndex][rxnListIndex] = rFactor;
	dorReactants->updateValue(m,rxnListIndex, rFactor);
}

Molecule ** DORrxn::pickReactants(double randA_Value)
{
	Molecule ** reactantSelections = new Molecule * [n_reactants];
	
	//Loop through the reactants
	for(int r=0; r<this->n_reactants; r++)
	{
		//If the reactant is a DOR reactant, we will have to select
		//which guy fires based on the rate factor
		if(r == DORreactantIndex)
		{
			double baseRateFactor = rate;
			for(int r2 =0; r2<this->n_reactants; r2++)
			{
				if(r2!=DORreactantIndex)
					baseRateFactor*=(lastIndex[r2]+1);
			}
			
			//cout<<" picking reactants for " << name;
			//dorReactants->printDetails();
			
			reactantSelections[r] = dorReactants->getReactantFromValue(randA_Value, baseRateFactor);
			
			
//			//First calculate contribution from other molecules for the rate of this guy
//			double a_value = rate;
//			for(int r2 =0; r2<this->n_reactants; r2++)
//			{
//				if(r2!=DORreactantIndex)
//					a_value*=(lastIndex[r2]+1);
//			}
//			
//			//Now, loop through all the molecules to get the molecule that should fire
//			//Start the loop from the beginning if we get a small randA_Value or the
//			//end if we get a large randA_Value
//			double cumRateFactorSum = 0;
//			if(randA_Value<=(rateFactorSum[r]*a_value/2))
//			{
//				//bool uninit = true;
//				for(int m=0; m<lastIndex[r]+1; m++)
//				{
//					cumRateFactorSum += a_value * rateFactors[r][m];
//					if(randA_Value<=cumRateFactorSum)
//					{
//						reactantSelections[r] = reactantArray[r][m];
//						break;
//					}
//				}
//			}
//			else
//			{
//				cumRateFactorSum = rateFactorSum[r]*a_value;
//				for(int m=lastIndex[r]; m>=0; m--)
//				{
//					cumRateFactorSum -= a_value * rateFactors[r][m];
//					if(randA_Value>=cumRateFactorSum)
//					{
//						reactantSelections[r] = reactantArray[r][m];
//						break;
//					}
//				}
//			}
		}
		
		//If it wasn't a DOR reactant, select the next molecule randomly
		else
		{
			int rand = NFutil::RANDOM_INT(0,(lastIndex[r]+1));
			reactantSelections[r] = reactantArray[r][rand];
		}
		
	}
	return reactantSelections;
}


double DORrxn::update_a()
{
	//When we calculate the new reaction propensity, it is just 
	a = rate;  //Start at 1 because rate factor sum already has rate multiplied in
	for(int r=0; r<n_reactants; r++)
	{
		if(r == DORreactantIndex)
		{
			//Here we can check the precision to make sure we don't have any
			//numerical errors from adding and subtracting from the rateFactorSum
			if(CHECK_PRECISION)
			{
				double cumRateFactorSum = 0;
				for(int m=0; m<lastIndex[r]+1; m++) cumRateFactorSum +=  rateFactors[r][m];
				if(fabs(cumRateFactorSum-rateFactorSum[r]) >= MAX_ERROR)
				{
					cerr<<"In DOR RXN: "<<name<<" NUMERICAL ERROR IN RATE FACTOR >= " << MAX_ERROR << ".  CORRECTING. "<<endl;
					rateFactorSum[r] = cumRateFactorSum;
				}
			}
			//cout<<" dor Rate Factor Sum: " << dorReactants->getRateFactorSum() << endl;
			double rfSum = dorReactants->getRateFactorSum();
			if(dorReactants->getNumOfMolecules()==0 || rfSum <=0 ) {a=0; return a; }
			a*= rfSum;
		}
		else
		{
			a*= (lastIndex[r]+1);
		}
	}
	return a;
}

void DORrxn::transformReactants(Molecule ** reactants, int nReactants)
{
	cerr<<"Calling transform function from a generalized DOR reaction!"<<endl;
}



void DORrxn::printDetails() const 
{
	cout<<"** DORReaction: " << name <<"  ( rate="<<rate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(int r=0; r<n_reactants; r++)
	{
		if(r==DORreactantIndex)
		{
			cout<<"      -(DOR) "<< this->reactantTemplates[r]->getMoleculeTypeName();
			cout<<"	(count="<< dorReactants->getNumOfMolecules() <<", RF Sum="<<dorReactants->getRateFactorSum()<<")."<<endl;		
		}
		else
		{
			cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
			cout<<"	(count="<< lastIndex[r]+1 <<")."<<endl;
		}
	}
}


/////////////////////////////////////////////////////





ObsDrxn::ObsDrxn( char * name, 
				int n_reactants, 
				TemplateMolecule ** reactantTemplates,
				int n_observables,
				Observable ** obs) : ReactionClass(name, n_reactants, reactantTemplates, 0)
{
	if(DEBUG) cout<<"Creating generalized Observable Dependent Reaction..."<<endl;
	this->obs = obs;
}


ObsDrxn::~ObsDrxn()
{
	delete [] obs;
}

double ObsDrxn::update_a()
{
	cerr<<" !!! Calling update_a from a generalized ObsDrxn!!!  ("<<name<<")"<<endl;
	
	a = 1;
	for(int r=0; r<n_reactants; r++)
		a*= (lastIndex[r]+1);
		
	a*=rate;
	return a;
}


void ObsDrxn::transformReactants(Molecule ** reactants, int nReactants)
{
	cerr<<" !!! Calling transformReactants from a generalized ObsDrxn!!!  ("<<name<<")"<<endl;
}


void ObsDrxn::printDetails() const
{
	cout<<"** ObsDrxn: " << name <<"  ( rate="<<rate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< lastIndex[r]+1 <<")."<<endl;
	}
	
}


















