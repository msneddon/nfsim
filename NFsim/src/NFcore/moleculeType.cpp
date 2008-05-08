#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;

MoleculeType::MoleculeType(
	string name, 
	string * stateNames, 
	int *defaultStateValues, 
	int numOfStates,
	string * bindingSiteNames,
	int numOfBindingSites,
	System * system )
{
	if(DEBUG) cout << "Creating MoleculeType " << name;
	
	this->name = name.c_str();
	
	//Set the state names and default states
	this->numOfStates = numOfStates;
	this->stateNames = stateNames;
	this->defaultStateValues = defaultStateValues;
	
	//Set the binding site names
	this->numOfBindingSites = numOfBindingSites;
	this->bindingSiteNames = bindingSiteNames;
	
	//Register myself with the system, and get an ID number
	this->system = system;
	this->type_id = this->system->addMoleculeType(this);
	
	
	if(DEBUG) {
		cout << "Creating MoleculeType " << name;
		cout << ": ("<< type_id <<") with " << numOfStates << " states (";
		for(int s=0; s<numOfStates; s++)
			cout<<" "<<this->stateNames[s]<<" ";
		cout << ") and " << numOfBindingSites << " binding sites (";
		for(int b=0; b<numOfBindingSites; b++)
			cout<<" "<<bindingSiteNames[b]<<" ";
		cout <<")."<<endl; }
}




MoleculeType::MoleculeType(
	const char * name, 
	char ** stateNames, 
	int *defaultStateValues, 
	int numOfStates,
	char ** bindingSiteNames,
	int numOfBindingSites,
	System * system )
{
/*	if(DEBUG) cout << "Creating MoleculeType " << name;
	
	string str(name);
	this->name = str.c_str();
	
	//Set the state names and default states
	this->numOfStates = numOfStates;
	this->stateNames = stateNames;
	this->defaultStateValues = defaultStateValues;
	
	//Set the binding site names
	this->numOfBindingSites = numOfBindingSites;
	this->bindingSiteNames = bindingSiteNames;
	
	//Register myself with the system, and get an ID number
	this->system = system;
	this->type_id = this->system->addMoleculeType(this);
	
	
	if(DEBUG) {
		cout << "("<< type_id <<") with " << numOfStates << " states (";
		for(int s=0; s<numOfStates; s++)
			cout<<" "<<this->stateNames[s]<<" ";
		cout << ") and " << numOfBindingSites << " binding sites (";
		for(int b=0; b<numOfBindingSites; b++)
			cout<<" "<<bindingSiteNames[b]<<" ";
		cout <<")."<<endl; }
*/
}

MoleculeType::~MoleculeType()
{
	if(DEBUG) cout << "Destroying MoleculeType " << name << endl;
	
	delete [] stateNames;
	delete [] defaultStateValues;
	delete [] bindingSiteNames;
	
	
	
	//Delete all molecules of this type that exist
	Molecule *m;
	while(mInstances.size()>0)
	{
		m = mInstances.back();
		mInstances.pop_back();
		delete m;
	}
	
	//Delete all template molecules of this type that exist
	TemplateMolecule *t;
	while(allTemplates.size()>0)
	{
		t = allTemplates.back();
		allTemplates.pop_back();
		delete t;
	}
	
	//Delete all observables of this type that exist
	Observable *o;
	while(observables.size()>0)
	{
		o = observables.back();
		observables.pop_back();
		delete o;
	}
}

int MoleculeType::addMolecule(Molecule *m)
{
	if(m->getMoleculeType()==this)
		mInstances.push_back(m);
	else {
		cout<<"!!!!Error: trying to add molecule of type " << m->getMoleculeTypeName() << " to MoleculeType " << name << endl;
		return -1;
	}
	return mInstances.size()-1;

}

void MoleculeType::addTemplateMolecule(TemplateMolecule *t)
{
	if(t->getMoleculeType()==this)
		allTemplates.push_back(t);
	else
		cout<<"!!!!Error: trying to add molecule of type " << t->getMoleculeTypeName() << " to MoleculeType " << name << endl;
}


const char * MoleculeType::getObservableAlias(int obsIndex) const 
{ 
	return observables.at(obsIndex)->getAliasName();
}

unsigned long int MoleculeType::getObservableCount(int obsIndex) const
{
	return observables.at(obsIndex)->getCount();
}

int MoleculeType::getBindingSiteIndex(const char * siteName ) const
{
	for(int b=0; b<numOfBindingSites; b++)
		if(siteName==bindingSiteNames[b]) return b;
	cerr<<"!!! warning !!! cannot find site name "<< siteName << " in MoleculeType: "<<name<<endl;
	this->printDetails();
	exit(1);
}

int MoleculeType::getStateIndex(const char * stateName ) const
{
	for(int s=0; s<numOfStates; s++)
		if(stateName==stateNames[s]) return s;
	cerr<<"!!! warning !!! cannot find state name "<< stateName << " in MoleculeType: "<<name<<endl;
	exit(1);
}


void MoleculeType::addReactionClass(ReactionClass * r, int rPosition)
{
	this->reactions.push_back(r);
	this->reactionPositions.push_back(rPosition);
}

void MoleculeType::addDORrxnClass(ReactionClass * r, int rPosition)
{
	//cout<<"adding dor rxn: "<<r->getName()<<" to moleculeType: "<<name<<endl;
	this->reactions.push_back(r);
	this->reactionPositions.push_back(rPosition);
	this->indexOfDORrxns.push_back(reactions.size()-1);
	//cout<<"index is "<< reactions.size()-1 << endl;
}

void MoleculeType::populateWithDefaultMolecules(int moleculeCount)
{
	if(DEBUG) cout<< " Populating "<< this->name << " with " << moleculeCount << " molecule(s)";
	if(DEBUG) cout<< " for a total of " << mInstances.size()+moleculeCount << " molecule(s)."<<endl;
	mInstances.reserve(mInstances.size()+moleculeCount);
	for(int m=0; m<moleculeCount; m++)
	{	
		if(DEBUG) cout<<" ("<<m+1<<") ";
		
		//Create the molecule (which knows how many states and binding sites to make)
		new Molecule(this);
		
		//Add the molecule to the list of molecules so we save it (does this automatically now!!!! )
		//mInstances.push_back(mol);
	}
}

void MoleculeType::prepareForSimulation()
{
	//Our iterators that we will use to loop through every molecule
	int r=0;
	
  	for( molIter = mInstances.begin(); molIter != mInstances.end(); molIter++ )
  	{
  		//First prepare the molecule for simulation
  		(*molIter)->prepareForSimulation();
  		
  		//Check each observable and see if this molecule should be counted
  		for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
  		{
			if((*obsIter)->isObservable(*molIter))
				(*obsIter)->add();
  		}
  		
  		//Check each reaction and add this molecule as a reactant if we have to
		for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
		{
			(*rxnIter)->tryToAdd((*molIter), reactionPositions.at(r));
  		}
	}
}

void MoleculeType::updateRxnMembership(Molecule * m)
{
	int r=0;
	
//	cout<<"update membership of "<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<endl;
	
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		(*rxnIter)->tryToAdd(m, reactionPositions.at(r));
  	}
  	
}


int MoleculeType::getRxnIndex(ReactionClass * rxn, int rxnPosition)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		if((*rxnIter)==rxn)
			if(reactionPositions.at(r) == rxnPosition)
				return r;
	}
	cerr<<"Could not find this rxn: " << rxn->getName() << " in molecule Type: "<<name<<endl;
	exit(1);
}




void MoleculeType::removeFromObservables(Molecule *m)
{
	//Check each observable and see if this molecule was counted, and if so, remove
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		if((*obsIter)->isObservable(m))
			(*obsIter)->subtract();
}

void MoleculeType::addToObservables(Molecule *m)
{
	//Check each observable and see if this molecule should be counted
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		if((*obsIter)->isObservable(m))
			(*obsIter)->add();
}

void MoleculeType::outputObservableNames(ofstream &fout)
{
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		fout<<"\t"<<(*obsIter)->getAliasName();
}

void MoleculeType::outputObservableCounts(ofstream &fout)
{
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		fout<<"\t"<<(*obsIter)->getCount();
}

void MoleculeType::printObservableCounts()
{
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		cout<<"\t"<<(*obsIter)->getCount();
}


void MoleculeType::printDetails() const
{
	cout<<"Molecule Type: "<< name << " type ID: " << type_id <<endl;
	cout<<"   -has "<< mInstances.size() <<" molecules."<<endl;
	cout<<"   -has "<< reactions.size() <<" reactions"<<endl;
	cout<<"        of which "<< indexOfDORrxns.size() <<" are DOR rxns. "<<endl;
	cout<<"   -has "<< observables.size() <<" observables " <<endl;
	
}





