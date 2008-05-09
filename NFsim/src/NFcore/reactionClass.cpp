#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;





ReactionClass::~ReactionClass()
{
	if(DEBUG) cout<<"Destorying rxn: "<<name<<endl;
	
	for(int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		reactantTemplates[r] = 0;
	}
	delete [] reactantTemplates;
	
	

	Transformation *tr;
	while(transformations.size()>0)
	{
		tr = transformations.back();
		transformations.pop_back();
		delete tr;
	}
	
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
	for(int r=0; r<n_reactants; r++)
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
	
	//No matter the situation, we should pop this molecule first
	rl->pop(m);
	
	//Try to match this molecule to the template at this position
	MappingSet *ms = new MappingSet();
	bool match = reactantTemplates[position]->compare(m,ms);
	
	//cout<<"comparing, resulting mappingSet:"<<endl;
	//ms->printDetails();
	
//	cout<<"in rxn "<<name<<", comparing: " << m->getMoleculeTypeName() << "_" <<m->getUniqueID()<<endl;
//	m->printDetails();
	
	//If there was a match, add the mappingSet
	if(match)
	{
	//	cout<<"Match!!"<<endl;
		rl->push(m,ms);
		return true;
	}
	
	
//	cout<<"Nope!!"<<endl;
	//Otherwise, do nothing and let the superiors know we didn't add
	delete ms;
	return false;
}
		

void ReactionClass::pickMappingSets(double random_A_number, vector <MappingSet *> &mappingSets)
{
	//Note here that we completely ignore the arguement.  The arguement is only
	//used for DOR reactions because we need that number to select the reactant to fire
	
	//So, we shall loop through the lists and extract out the MappingSets and 
	//the molecule ID that the mappingSet was generated from
	unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int i=0; i<n_reactants; i++)
	{
		MappingSet *ms= NULL;
		reactantLists.at(i)->pickRandom(ms,moleculeIDs[i]);
		//cout<<"here: "<<endl;
		//ms->printDetails();
		mappingSets.push_back(ms);
		
	//	cout<<"selected: \t"<<moleculeIDs[i]<<endl;
	}
	
	delete [] moleculeIDs;
}



void ReactionClass::fire2(double random_A_number)
{
	fireCounter++; //Remember that we fired
//	cout<<"firing: "<<name<<endl;
	//First randomly pick the reactants to fire by selecting the MappingSets
	vector <MappingSet *> mappingSets;
	pickMappingSets(random_A_number, mappingSets);
	
	//Only continue if we found some MappingSets.  PickMappingSets may return
	//nothing if certain constraints are not met, so we have to check here
	if(mappingSets.size()==0)
	{
		return;
	}
	
	//Generate the set of possible products that we need to update
	list <Molecule *> products;
	MappingSet::getPossibleProducts(mappingSets,products,traversalLimit);
	
	//Loop through the products and remove them from thier observables
	list <Molecule *>::iterator molIter;
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
		(*molIter)->removeFromObservables();
		
	//Through the MappingSet, transform all the molecules as neccessary
	MappingSet::transform(mappingSets);
		

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
	mappingSets.clear();
}



void ReactionClass::printFullDetails()
{
	cout<<"Reaction: "<<name<<endl;
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists.at(i)->printDetails();
	
	
	
}





ReactionClass::ReactionClass(string name, vector <TemplateMolecule *> templateMolecules, double rate)
{
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->rate = rate;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	
	//Set up the template molecules
	this->n_reactants = templateMolecules.size();
	this->reactantTemplates = new TemplateMolecule *[n_reactants];
	for(unsigned int i=0; i<n_reactants; i++)
		reactantTemplates[i] = templateMolecules.at(i);
	
	//Set up the reactantLists
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists.push_back(new ReactantList(50));
	
	
	
	
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
	
	
	
	
	this->reactionType = NORMAL_RXN;  //set as normal reaction here, but deriving reaction classes can change this
}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


ReactantList::ReactantList(unsigned int init_capacity)
{
	this->n_reactants = 0;
	this->capacity = init_capacity;
	this->reactants = new MappingSet *[init_capacity];
	this->reactantMoleculeIDs = new unsigned int [init_capacity];
	for(unsigned int i=0; i<capacity; i++)
	{
		reactants[i]=NULL;
		reactantMoleculeIDs[i]=0;
	}
}

ReactantList::~ReactantList()
{
	moleculeLookupTable.clear();
	
	for(unsigned int i=0; i<capacity; i++)
		if(reactants[i]!=NULL)
		{
			reactants[i]->printDetails();
			delete reactants[i];
		}
	
	delete [] reactants;
	delete [] reactantMoleculeIDs;
}


void ReactantList::pickRandom(MappingSet *&ms, unsigned int &moleculeID)
{
	unsigned int rand = NFutil::RANDOM_INT(0,n_reactants);
	ms = reactants[rand];
	moleculeID = reactantMoleculeIDs[rand];
}



unsigned int ReactantList::size()
{
	return n_reactants;
}



void ReactantList::push(Molecule *reactant, MappingSet *mappingSet)
{
	if(n_reactants<capacity)
	{
		reactantMoleculeIDs[n_reactants] = reactant->getUniqueID();
		moleculeLookupTable.insert(pair<unsigned int, unsigned int>(reactantMoleculeIDs[n_reactants],n_reactants));
		reactants[n_reactants] = mappingSet;
		
		n_reactants++;
	}
	else
	{
		MappingSet ** newMapArray = new MappingSet *[capacity*2];
		unsigned int * newMolArray = new unsigned int [capacity*2];
		for(unsigned int i=0; i<capacity; i++)
		{
			newMapArray[i] = reactants[i];
			newMolArray[i] = reactantMoleculeIDs[i];
		}
		for(unsigned int i=capacity; i<capacity*2; i++)
		{
			newMapArray[i] = NULL;
			newMolArray[i]=0;
		}
		
		delete [] reactants;
		delete [] reactantMoleculeIDs;
		reactants = newMapArray;
		reactantMoleculeIDs = newMolArray;
		capacity*=2;
		
		reactantMoleculeIDs[n_reactants] = reactant->getUniqueID();
		moleculeLookupTable.insert(pair<unsigned int, unsigned int>(reactantMoleculeIDs[n_reactants],n_reactants));
		reactants[n_reactants] = mappingSet;
		n_reactants++;
	}
	
}


void ReactantList::pop(Molecule *m)
{
	multimap <unsigned int, unsigned int>::iterator it;
	for(it = moleculeLookupTable.find(m->getUniqueID());
		it != moleculeLookupTable.end();
		it = moleculeLookupTable.find(m->getUniqueID()) )
	{
		this->removeReactant(it->second);
		moleculeLookupTable.erase(it);
	}
}


void ReactantList::removeReactant(unsigned int index)
{
	
	//If the array has only one element, or we happened to select the last element,
	//Then just remove the last element
	if(index == (n_reactants-1) )
	{
		delete reactants[index];
		reactants[index] = NULL;
		reactantMoleculeIDs[index]=0;
	}
	
	//Otherwise we have to swap with the last element in the list
	else
	{
		MappingSet *ms = reactants[index];
		reactants[index] = reactants[n_reactants-1];
		reactantMoleculeIDs[index] = reactantMoleculeIDs[n_reactants-1] ;
		
		pair<multimap<unsigned int,unsigned int>::iterator,multimap<unsigned int,unsigned int>::iterator> ret;
		multimap<unsigned int,unsigned int>::iterator it;
		
		ret = moleculeLookupTable.equal_range(reactantMoleculeIDs[n_reactants-1]);
		for(it = ret.first;it != ret.second;it++)
			if(it->second == (n_reactants-1))
			{
				it->second = index;
				break;
			}
		
		reactants[n_reactants-1] = NULL;
		reactantMoleculeIDs[n_reactants-1] = 0;
		delete ms;
	}
	
	n_reactants--;
}



void ReactantList::printDetails()
{
	cout<<"ReactantList that contains: "<<size()<<" MappingSets and has a capacity for "<<capacity<<" total sets."<<endl;
	cout<<"Members: [molecule unique Id, private mappingSet index]"<<endl;
	
	
	for(unsigned int i=0; i<capacity; i++)
		cout<<i<<" ";
	cout<<endl;
	for(unsigned i=0; i<capacity; i++)
		if(reactants[i]!=NULL)
			cout<<1<<" ";
		else
			cout<<0<<" ";
	cout<<endl;
	for(unsigned i=0; i<capacity; i++)
			cout<<reactantMoleculeIDs[i]<<" ";
		cout<<endl;
	
	multimap <unsigned int, unsigned int>::iterator it;
	for ( it=moleculeLookupTable.begin(); it != moleculeLookupTable.end(); it++ )
	   cout << "  [" << (*it).first << ", " << (*it).second << "]" << endl;
	
	cout<<endl;
}









