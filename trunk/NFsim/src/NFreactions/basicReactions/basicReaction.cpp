



#include "basicReaction.hh"


using namespace std;
using namespace NFcore;




FunctionalRxnClass::FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet) :
	BasicRxnClass(name,1, transformationSet)
{
	this->gf=gf;
}
FunctionalRxnClass::~FunctionalRxnClass() {};
			
double FunctionalRxnClass::update_a() {
	a = 1;
	for(unsigned int i=0; i<n_reactants; i++)
		a*=reactantLists.at(i)->size();
	
	a*=FuncFactory::Eval(gf->p);
	return a;
}

void FunctionalRxnClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( baseFunction="<<gf->getNiceName()<<"="<<FuncFactory::Eval(gf->p)<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}



BasicRxnClass::BasicRxnClass(string name, double baseRate, TransformationSet *transformationSet) : 
	ReactionClass(name,baseRate,transformationSet)
{
	this->reactionType = BASIC_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	
	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists.push_back(new ReactantList(r,transformationSet,25));
}


BasicRxnClass::~BasicRxnClass()
{
	
//	this->reactantLists.at(0)->printDetails();
	
	if(DEBUG) cout<<"Destorying rxn: "<<name<<endl;
	
	for(unsigned int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		reactantTemplates[r] = 0;
	}
	delete [] reactantTemplates;
	delete [] mappingSet;
	
/*
	Transformation *tr;
	while(transformations.size()>0)
	{
		tr = transformations.back();
		transformations.pop_back();
		delete tr;
	}
	*/
	
	ReactantList *rl;
	while(reactantLists.size()>0)
	{
		rl = reactantLists.back();
		reactantLists.pop_back();
		delete rl;
	}
	
	delete transformationSet;
}
			
void BasicRxnClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void BasicRxnClass::prepareForSimulation() 
{
	
}



bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL) 
	{
		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}
	
	
	//Get the specified reactantList
	ReactantList *rl = reactantLists.at(reactantPos);
	
	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//cout<<"got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;
	
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);
	
	
	if(isInRxn)
	{
		if(!reactantTemplates[reactantPos]->compare(m)) {
		//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
		
	} else {
		//Try to map it!
		MappingSet *ms = rl->pushNextAvailableMappingSet();
		if(!reactantTemplates[reactantPos]->compare(m,ms)) {
			rl->popLastMappingSet();
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}
		
	return true;
}


void BasicRxnClass::remove(Molecule *m, unsigned int reactantPos)
{
	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL) 
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}
		
		
	//Get the specified reactantList
	ReactantList *rl = reactantLists.at(reactantPos);
		
	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);
		
		
	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}





						
						
void BasicRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) 
{
	cerr<<"You are trying to notify a Basic Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}


unsigned int BasicRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return reactantLists.at(reactantIndex)->size();
}			



void BasicRxnClass::printFullDetails() const
{
	cout<<"BasicRxnClass: "<<name<<endl;
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists.at(i)->printDetails();
}

			
void BasicRxnClass::pickMappingSets(double random_A_number) const
{
	//Note here that we completely ignore the arguement.  The arguement is only
	//used for DOR reactions because we need that number to select the reactant to fire
	
	//So, we shall loop through the lists and extract out the MappingSets and 
	//the molecule ID that the mappingSet was generated from
	//unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int i=0; i<n_reactants; i++)
	{
		reactantLists.at(i)->pickRandom(mappingSet[i]);
		//mappingSet[i]->get(0)
	}
	
	//delete [] moleculeIDs;
}
































