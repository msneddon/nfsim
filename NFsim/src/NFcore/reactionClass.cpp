


#include "NFcore.hh"


using namespace std;
using namespace NFcore;



ReactionClass::ReactionClass(string name, double baseRate, TransformationSet *transformationSet)
{
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->baseRate = baseRate;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;
	
	
	//Set up the template molecules from the transformationSet
	this->n_reactants = transformationSet->getNreactants();
	this->reactantTemplates = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		reactantTemplates[r] = transformationSet->getTemplateMolecule(r);
	mappingSet = new MappingSet *[n_reactants];
	
	
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
			baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
		}	
	}
}


ReactionClass::~ReactionClass()
{
	
}



void ReactionClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}


void ReactionClass::fire(double random_A_number)
{
	
	fireCounter++; //Remember that we fired
//	cout<<"firing: "<<name<<endl;
	
	//First randomly pick the reactants to fire by selecting the MappingSets
	pickMappingSets(random_A_number);

	//Generate the set of possible products that we need to update
	this->transformationSet->getListOfProducts(mappingSet,products,traversalLimit);
	
	//Loop through the products and remove them from thier observables
	//cout<<"------------------------------------------"<<endl;
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	{
		//cout<<"Removing: "<<(*molIter)->getMoleculeTypeName()<<"_"<<(*molIter)->getMoleculeID()<<endl;
		(*molIter)->removeFromObservables();
	}
	//cout<<"++++++++++++++++++++++++++++++++++++++++++"<<endl;
		
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











