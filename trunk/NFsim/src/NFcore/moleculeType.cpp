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
	this->isStateEq = new int [numOfStates];
	for(int s=0; s<numOfStates; s++) {
		isStateEq[s] = -1;
	}
	
	//Set the binding site names
	this->numOfBindingSites = numOfBindingSites;
	this->bindingSiteNames = bindingSiteNames;
	this->isSiteEq = new int [numOfBindingSites];
	for(int b=0; b<numOfBindingSites; b++) {
		isSiteEq[b] = -1;
	}
	
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
	
	
	
	mList = new MoleculeList(this,2);
	n_eqSites = 0;
	n_eqStates = 0;
}



MoleculeType::~MoleculeType()
{
	if(DEBUG) cout << "Destroying MoleculeType " << name << endl;
	
	delete [] stateNames;
	delete [] defaultStateValues;
	delete [] bindingSiteNames;
	
	delete [] isSiteEq;
	delete [] eqSiteSizes;
	for(unsigned int i=0; i<n_eqSites; i++) {
		delete [] eqSiteName[i];
		delete [] eqSiteIndex[i];
	}
	delete [] eqSiteName;
	delete [] eqSiteIndex;
	
	
	delete [] isStateEq;
	delete [] eqStateSizes;
	for(unsigned int i=0; i<n_eqStates; i++) {
		delete [] eqStateName[i];
		delete [] eqStateIndex[i];
	}
	delete [] eqStateName;
	delete [] eqStateIndex;
		
	
	//Delete all molecules of this type that exist
	//Molecule *m;
	//while(.size()>0)
	//{
	//	m = mInstances.back();
	//	mInstances.pop_back();
	//	delete m;
	//}
	
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
	
	
	delete mList;
}


void MoleculeType::addEquivalentSites(vector <vector <string> > &identicalSites)
{
	this->n_eqSites = identicalSites.size();
	eqSiteName=new string * [n_eqSites];
	eqSiteIndex=new int *[n_eqSites];
	eqSiteSizes=new unsigned int [n_eqSites];
	
	for(unsigned int i=0; i<n_eqSites; i++) {
		eqSiteSizes[i]=identicalSites.at(i).size();
		eqSiteName[i] = new string [eqSiteSizes[i]];
		eqSiteIndex[i] = new int [eqSiteSizes[i]];
		for(unsigned int k=0; k<eqSiteSizes[i]; k++) {
			eqSiteName[i][k] = identicalSites.at(i).at(k);
			eqSiteIndex[i][k] = getBindingSiteIndex(eqSiteName[i][k]);
			isSiteEq[eqSiteIndex[i][k]]=i;
		}
	}
}


bool MoleculeType::isAnEquivalentSite(string siteName) {
	for(unsigned int i=0; i<n_eqSites; i++) {
		for(unsigned int k=0; k<eqSiteSizes[i]; k++) {
			if(eqSiteName[i][k]==siteName)
				return true;
		}
	}
	return false;
}
bool MoleculeType::isAnEquivalentSite(int siteIndex) {
	for(unsigned int i=0; i<n_eqSites; i++) {
		for(unsigned int k=0; k<eqSiteSizes[i]; k++) {
			if(eqSiteIndex[i][k]==siteIndex)
				return true;
		}
	}
	return false;
}
void MoleculeType::getEquivalencyClass(int *&sites, int &n_sites, string siteName) {
	for(unsigned int i=0; i<n_eqSites; i++) {
		for(unsigned int k=0; k<eqSiteSizes[i]; k++) {
			if(eqSiteName[i][k]==siteName) {
				sites = eqSiteIndex[i];
				n_sites=eqSiteSizes[i];
			}
		}
	}
}


void MoleculeType::addEquivalentStates(vector <vector <string> > &identicalStates)
{
	this->n_eqStates = identicalStates.size();
	eqStateName=new string * [n_eqStates];
	eqStateIndex=new int *[n_eqStates];
	eqStateSizes=new unsigned int [n_eqStates];
	
	for(unsigned int i=0; i<n_eqStates; i++) {
		eqStateSizes[i]=identicalStates.at(i).size();
		eqStateName[i] = new string [eqStateSizes[i]];
		eqStateIndex[i] = new int [eqStateSizes[i]];

		for(unsigned int k=0; k<eqStateSizes[i]; k++) {
			eqStateName[i][k] = identicalStates.at(i).at(k);
			eqStateIndex[i][k] = getStateIndex(eqStateName[i][k]);
			isStateEq[eqStateIndex[i][k]]=i;
		}
	}
}



bool MoleculeType::isAnEquivalentState(string stateName) {
	for(unsigned int i=0; i<n_eqStates; i++) {
		for(unsigned int k=0; k<eqStateSizes[i]; k++) {
			if(eqStateName[i][k]==stateName)
				return true;
		}
	}
	return false;
}
bool MoleculeType::isAnEquivalentState(int stateIndex) {
	for(unsigned int i=0; i<n_eqStates; i++) {
		for(unsigned int k=0; k<eqStateSizes[i]; k++) {
			if(eqStateIndex[i][k]==stateIndex)
				return true;
		}
	}
	return false;
}
void MoleculeType::getEquivalencyStateClass(int *&states, int &n_states, string stateName) {
	for(unsigned int i=0; i<n_eqStates; i++) {
		for(unsigned int k=0; k<eqStateSizes[i]; k++) {
			if(eqStateName[i][k]==stateName) {
				states = eqStateIndex[i];
				n_states=eqStateSizes[i];
			}
		}
	}
}


Molecule *MoleculeType::genDefaultMolecule()
{
	Molecule *m;
	mList->create(m);
	return m;
}


void MoleculeType::addMoleculeToRunningSystem(Molecule *&mol)
{
	//First prepare the molecule for simulation
	mol->prepareForSimulation();
	  		
	//Check each observable and see if this molecule should be counted
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
	{
		if((*obsIter)->isObservable(mol))
			(*obsIter)->add();
	}
	  		
	//Check each reaction and add this molecule as a reactant if we have to
	int r=0;
	for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
	{
		(*rxnIter)->tryToAdd(mol, reactionPositions.at(r));
	}
}



void MoleculeType::removeMoleculeFromRunningSystem(Molecule *&m)
{
	mList->remove(m->getMolListId());
	removeFromObservables(m);
	removeFromRxns(m);
	
}


Molecule * MoleculeType::getMolecule(int ID_molecule) const { return mList->at(ID_molecule); }
int MoleculeType::getMoleculeCount() const { return mList->size(); }


void MoleculeType::addTemplateMolecule(TemplateMolecule *t)
{
	if(t->getMoleculeType()==this)
		allTemplates.push_back(t);
	else
		cout<<"!!!!Error: trying to add molecule of type " << t->getMoleculeTypeName() << " to MoleculeType " << name << endl;
}


string MoleculeType::getObservableAlias(int obsIndex) const 
{ 
	return observables.at(obsIndex)->getAliasName();
}

unsigned long int MoleculeType::getObservableCount(int obsIndex) const
{
	return observables.at(obsIndex)->getCount();
}

unsigned int MoleculeType::getBindingSiteIndex(string siteName ) const
{
	for(int b=0; b<numOfBindingSites; b++)
		if(siteName==bindingSiteNames[b]) return b;
	cerr<<"!!! warning !!! cannot find site name "<< siteName << " in MoleculeType: "<<name<<endl;
	this->printDetails();
	exit(1);
}

unsigned int MoleculeType::getStateIndex(string stateName ) const
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
	if(DEBUG) cout<< " for a total of " << mList->size()+moleculeCount << " molecule(s)."<<endl;
	//mInstances.reserve(mInstances.size()+moleculeCount);
	for(int m=0; m<moleculeCount; m++)
	{	
		if(DEBUG) cout<<" ("<<m+1<<") ";
		
		//Create the molecule (which knows how many states and binding sites to make)
		this->genDefaultMolecule();
		//new Molecule(this);
		
		//Add the molecule to the list of molecules so we save it (does this automatically now!!!! )
		//mInstances.push_back(mol);
	}
}





void MoleculeType::prepareForSimulation()
{
	
	//Check each reaction and add this molecule as a reactant if we have to
	int r=0;
	for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
	{
		system->registerRxnIndex((*rxnIter)->getRxnId(), reactionPositions.at(r),r);
  	}
	
	//Our iterators that we will use to loop through every molecule
	
	Molecule *mol;
  	for( int m=0; m<mList->size(); m++ )
  	{
  		//First prepare the molecule for simulation
  		mol = mList->at(m);
  		mol->prepareForSimulation();
  		
  		//Check each observable and see if this molecule should be counted
  		for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
  		{
			if((*obsIter)->isObservable(mol))
				(*obsIter)->add();
  		}
  		
  		//Check each reaction and add this molecule as a reactant if we have to
		for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
		{
			(*rxnIter)->tryToAdd(mol, reactionPositions.at(r));
  		}
	}
}

void MoleculeType::updateRxnMembership(Molecule * m)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		(*rxnIter)->tryToAdd(m, reactionPositions.at(r));
  	}
}


int MoleculeType::getRxnIndex(ReactionClass * rxn, int rxnPosition)
{
	return system->getRxnIndex(rxn->getRxnId(),rxnPosition);
	
	//The old way!!  (that is slow if we have many rxns of course!)
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
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ ){
		//cout<<"Comparing (in subtract): "<<endl;
		//m->printDetails();
		if((*obsIter)->isObservable(m)) 
			(*obsIter)->subtract();
	}
}

void MoleculeType::removeFromRxns(Molecule * m)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		(*rxnIter)->remove(m, reactionPositions.at(r));
  	}
}

void MoleculeType::addToObservables(Molecule *m)
{
	//Check each observable and see if this molecule should be counted
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ ){
		//cout<<"Comparing(in add: "<<endl;
		//m->printDetails();
		if((*obsIter)->isObservable(m))
			(*obsIter)->add();
	}
			
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

void MoleculeType::printObservableNames()
{
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		cout<<"\t"<<(*obsIter)->getAliasName();
}

void MoleculeType::printObservableCounts()
{
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ )
		cout<<"\t"<<(*obsIter)->getCount();
}


void MoleculeType::printDetails() const
{
	cout<<"Molecule Type: "<< name << " type ID: " << type_id <<endl;
	
	cout<<"   -states ( ";
	for(int s=0; s<numOfStates; s++) cout<<stateNames[s]<<" ";
	cout<<")"<<endl;
	
	cout<<"   -sites ( ";
	for(int s=0; s<numOfBindingSites; s++) cout<<bindingSiteNames[s]<<" ";
	cout<<")"<<endl;
	
	for(unsigned int i=0; i<n_eqSites; i++) {
		cout<<"   -equivalent sites ( ";
		for(unsigned int k=0; k<eqSiteSizes[i]; k++)
			cout<<eqSiteName[i][k]<<" ";
		cout<<")"<<endl;
	}
	
	for(unsigned int i=0; i<n_eqStates; i++) {
		cout<<"   -equivalent states ( ";
		for(unsigned int k=0; k<eqStateSizes[i]; k++)
			cout<<eqStateName[i][k]<<" ";
		cout<<")"<<endl;
	}
	
	cout<<"   -has "<< mList->size() <<" molecules."<<endl;
	cout<<"   -has "<< reactions.size() <<" reactions"<<endl;
	cout<<"        of which "<< indexOfDORrxns.size() <<" are DOR rxns. "<<endl;
	cout<<"   -has "<< observables.size() <<" observables " <<endl;
}





