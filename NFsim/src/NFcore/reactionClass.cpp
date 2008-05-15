#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;





ReactionClass::~ReactionClass()
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
}
	

void ReactionClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if(reactionType==NORMAL_RXN)
			reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void ReactionClass::prepareForSimulation() {}

		


		
		

		

		
void ReactionClass::printDetails() const {
	cout<<"** Reaction: " << name <<"  ( rate="<<rate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< getReactantCount(r) <<")."<<endl;
	}
}
	












////////////////////
/////////////////////
////////////////////
/////////////////////
////////////////////
/////////////////////
////////////////////
/////////////////////
////////////////////
/////////////////////
////////////////////
/////////////////////
////////////////////
/////////////////////



unsigned int ReactionClass::getReactantCount(unsigned int reactantIndex) const
{
	return reactantLists.at(reactantIndex)->size();
}


double ReactionClass::update_a()
{
	a = 1;
	
	for(unsigned int i=0; i<n_reactants; i++)
		a*=reactantLists.at(i)->size();
	
	
	a*=rate;
	return a;
}



bool ReactionClass::tryToAdd(Molecule *m, unsigned int position)
{
	//First a bit of error checking...
	if(position<0 || position>=n_reactants || m==NULL) 
	{
		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}
	
	
	//Get the specified reactantList
	ReactantList *rl = reactantLists.at(position);
	
	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,position);
	//cout<<"got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;
	
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);
	
	
	if(isInRxn)
	{
		if(!reactantTemplates[position]->compare(m)) {
		//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,-1);
		}
		
	} else {
		//Try to map it!
		MappingSet *ms = rl->pushNextAvailableMappingSet();
		if(!reactantTemplates[position]->compare(m,ms)) {
			rl->popLastMappingSet();
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}
		
	return false;
}
		

void ReactionClass::pickMappingSets(double random_A_number)
{
	//Note here that we completely ignore the arguement.  The arguement is only
	//used for DOR reactions because we need that number to select the reactant to fire
	
	//So, we shall loop through the lists and extract out the MappingSets and 
	//the molecule ID that the mappingSet was generated from
	unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int i=0; i<n_reactants; i++)
	{
		reactantLists.at(i)->pickRandom(mappingSet[i]);
		//mappingSet[i]->get(0)
	}
	
	delete [] moleculeIDs;
}



void ReactionClass::fire2(double random_A_number)
{
	
	fireCounter++; //Remember that we fired
//	cout<<"firing: "<<name<<endl;
	
	//First randomly pick the reactants to fire by selecting the MappingSets
	pickMappingSets(random_A_number);

	//Generate the set of possible products that we need to update
	list <Molecule *> products;
	this->transformationSet->getListOfProducts(mappingSet,products,traversalLimit);
	
	//Loop through the products and remove them from thier observables
	list <Molecule *>::iterator molIter;
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
		(*molIter)->removeFromObservables();
		
	//Through the MappingSet, transform all the molecules as neccessary
	this->transformationSet->transform(this->mappingSet);

	//Tell each molecule in the list of products to add itself back into 
	//the counts of observables and update its class lists, and update any DOR Groups
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	{
		(*molIter)->addToObservables();
	  	(*molIter)->updateRxnMembership();
	  	(*molIter)->notifyGroupsThatRateMayChange();
	}
//	Molecule::printMoleculeList(products);
	
	//Tidy up
	products.clear();
}



void ReactionClass::printFullDetails()
{
	cout<<"Reaction: "<<name<<endl;
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists.at(i)->printDetails();
}




ReactionClass::ReactionClass(string name, double rate, TransformationSet *transformationSet)
{
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->rate = rate;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;
	
	
	//Set up the template molecules from the transformationSet
	this->n_reactants = transformationSet->getNreactants();
	this->reactantTemplates = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		reactantTemplates[r] = transformationSet->getTemplateMolecule(r);
	

	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists.push_back(new ReactantList(r,transformationSet,50));
	
	
	mappingSet = new MappingSet *[n_reactants];
	
	/*
	
	//Check here to see if we have molecule types that are the same across different reactants
	//Because if so, we will give a warning
	if(n_reactants>2) cerr<<"Warning!! You created a reaction ("<< name <<") that has more than 2 reactants.  This has not been tested!"<<endl;
	if(n_reactants==2)
	{
		//If the reactants are of the same type, then we have to make a few special considerations
		if(reactantTemplates[0]->getMoleculeType()==reactantTemplates[1]->getMoleculeType())
		{
			cout<<endl;
			cout<<"Warning! You have a rxn (" << name << ") that allows a moleculeType to bind another of the same type."<<endl;
			cout<<"Make sure that is correct."<<endl;
			cout<<endl;
			rate = rate*0.5;  //We have to correct the rate to get the proper factor
		}	
	}
	
	
	*/
	
	this->reactionType = NORMAL_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	
	
}
















