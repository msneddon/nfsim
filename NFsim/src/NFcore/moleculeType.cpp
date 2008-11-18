#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;




MoleculeType::MoleculeType(
	string name,
	vector <string> &compName,
	System *s)
{
	vector <string> defaultCompState;
	vector < vector <string> > possibleCompStates;
	vector <bool> isIntegerComponent;
	for(unsigned int i=0; i<compName.size(); i++) {
		vector <string> v;
		possibleCompStates.push_back(v);
		defaultCompState.push_back("NO_STATE");
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, s);

}

MoleculeType::MoleculeType(
	string name,
	vector <string> &compName,
	vector <string> &defaultCompState,
	System *s)
{
	vector < vector <string> > possibleCompStates;
	vector <bool> isIntegerComponent;
	for(unsigned int i=0; i<compName.size(); i++) {
		vector <string> v;
		possibleCompStates.push_back(v);
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, s);
}







MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		System *system)
{
	vector <bool> isIntegerComponent;
	for(unsigned int i=0; i<compName.size(); i++) {
		isIntegerComponent.push_back(false);
	}
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}

MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		vector <bool> isIntegerComponent,
		System *system)
{
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent,system);
}


void MoleculeType::init(
	string name,
	vector <string> &compName,
	vector <string> &defaultCompState,
	vector < vector<string> > &possibleCompStates,
	vector <bool> isIntegerComponent,
	System *system)
{
	//Basics...
	this->name=name;
	this->numOfComponents=compName.size();

	//First, some quick error checks
	if((int)defaultCompState.size()!=numOfComponents || (int)possibleCompStates.size()!=numOfComponents ||
			(int)isIntegerComponent.size()!=numOfComponents) {
		cout<<"Error creating MoleculeType: '"<<name<<"': The length of the input vectors\n";
		cout<<"do not match, so I can't initialize this object.\n";
		cout<<"quitting now."<<endl; exit(1);
	}

	//Now we can get on with initializing the MoleculeType information
	this->compName=new string [numOfComponents];
	this->defaultCompState = new int [numOfComponents];
	this->isIntegerCompState = new bool [numOfComponents];

	for(int c=0; c<numOfComponents; c++) {
		this->compName[c]=compName.at(c);
		this->isIntegerCompState[c]=isIntegerComponent.at(c);

		bool foundDefaultState=false;
		vector <string> p;
		for(unsigned int i=0; i<possibleCompStates.at(c).size(); i++) {
			p.push_back(possibleCompStates.at(c).at(i));
			if(possibleCompStates.at(c).at(i) == defaultCompState.at(c)) {
				this->defaultCompState[c]=i; foundDefaultState=true;
			}
		}
		if(!foundDefaultState) this->defaultCompState[c]=Molecule::NOSTATE;
		this->possibleCompStates.push_back(p);
	}


	//Register myself with the system, and get an ID number
	this->system = system;
	this->type_id = this->system->addMoleculeType(this);


	mList = new MoleculeList(this,2);
	n_eqComp = 0;
}






MoleculeType::~MoleculeType()
{
	if(DEBUG) cout << "Destroying MoleculeType " << name << endl;

	//Delete freestore component information
	delete [] compName;
	delete [] defaultCompState;
	delete [] isIntegerCompState;

	//Delete details about equivalent components
	delete [] eqCompSizes;
	for(int i=0; i<n_eqComp; i++) {
		delete [] eqCompName[i];
		delete [] eqCompIndex[i];
	}
	delete [] eqCompName;
	delete [] eqCompIndex;




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

void MoleculeType::addEquivalentComponents(vector <vector <string> > &identicalComponents)
{
	this->n_eqComp = identicalComponents.size();
	eqCompName=new string * [n_eqComp];
	eqCompIndex=new int *[n_eqComp];
	eqCompSizes=new int [n_eqComp];

	for(int i=0; i<n_eqComp; i++) {
		eqCompSizes[i]=identicalComponents.at(i).size();
		eqCompName[i] = new string [eqCompSizes[i]];
		eqCompIndex[i] = new int [eqCompSizes[i]];
		for(int k=0; k<eqCompSizes[i]; k++) {
			eqCompName[i][k] = identicalComponents.at(i).at(k);
			eqCompIndex[i][k] = getCompIndexFromName(eqCompName[i][k]);
		}
	}
}


bool MoleculeType::isIntegerComponent(string cName) const {
	for(int c=0; c<numOfComponents; c++)
			if(cName==compName[c]) {
				return this->isIntegerCompState[c];
			}
	cerr<<"!!! error !!! cannot find site name "<< cName << " in MoleculeType: "<<name;
	cerr<<"in function isIntegerComponent(string cName).  "<<endl;
	this->printDetails();
	exit(1);
}
bool MoleculeType::isIntegerComponent(int cIndex) const {
	if(cIndex>=0 && cIndex<numOfComponents) {
		return this->isIntegerCompState[cIndex];
	} else {
		cerr<<"!!! error !!! "<< cIndex << " is not a valid component index in MoleculeType: "<<name;
		cerr<<"in function isIntegerComponent(int cIndex).  "<<endl;
		this->printDetails();
		exit(1);
	}
}


bool MoleculeType::isEquivalentComponent(string cName) const {
	for(int i=0; i<n_eqComp; i++) {
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(eqCompName[i][k]==cName)
				return true;
		}
	}
	return false;
}
bool MoleculeType::isEquivalentComponent(int cIndex) const {
	for(int i=0; i<n_eqComp; i++) {
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(eqCompIndex[i][k]==cIndex)
				return true;
		}
	}
	return false;
}
void MoleculeType::getEquivalencyClass(int *&components, int &n_components, string cName) const {
	for(int i=0; i<n_eqComp; i++) {
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(eqCompName[i][k]==cName) {
				components = eqCompIndex[i];
				n_components=eqCompSizes[i];
			}
		}
	}
}





string MoleculeType::getComponentStateName(int cIndex, int cValue) {
	if(cValue==Molecule::NOSTATE) return "NO_STATE";
	return possibleCompStates.at(cIndex).at(cValue);
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
	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ ) {
		if((*obsIter)->isObservable(mol))
			(*obsIter)->add();
	}

	//Check each reaction and add this molecule as a reactant if we have to
	int r=0;
	for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ ) {
		(*rxnIter)->tryToAdd(mol, reactionPositions.at(r));
	}

}



void MoleculeType::removeMoleculeFromRunningSystem(Molecule *&m)
{
	mList->remove(m->getMolListId());
	removeFromObservables(m);
	removeFromRxns(m);
}


Molecule * MoleculeType::getMolecule(int ID_molecule) const {
	return mList->at(ID_molecule);
}
int MoleculeType::getMoleculeCount() const {
	return mList->size();
}


void MoleculeType::addTemplateMolecule(TemplateMolecule *t)
{
	if(t->getMoleculeType()==this)
		allTemplates.push_back(t);
	else
		cout<<"!!!!Error: trying to add molecule of type " << t->getMoleculeTypeName() << " to MoleculeType " << name << endl;
}


string MoleculeType::getObservableAlias(int obsIndex) const {
	return observables.at(obsIndex)->getAliasName();
}

unsigned long int MoleculeType::getObservableCount(int obsIndex) const {
	return observables.at(obsIndex)->getCount();
}


int MoleculeType::getCompIndexFromName(string cName) const
{
	for(int c=0; c<numOfComponents; c++)
		if(cName==compName[c]) return c;
	cerr<<"!!! warning !!! cannot find site name "<< cName << " in MoleculeType: "<<name<<endl;
	this->printDetails();
	exit(1);
}

int MoleculeType::getStateValueFromName(int cIndex, string stateName) const
{
	for(unsigned int s=0; s<possibleCompStates.at(cIndex).size(); s++) {
		if(possibleCompStates.at(cIndex).at(s)==stateName) {
			return s;
		}
	}
	cerr<<"Error!  '"<<stateName<<" is not a recognized possible state for '"<<compName[cIndex]<<"' in MoleculeType: '"<<name<<"'"<<endl;
	cerr<<"For that, I'm quitting!";
	printDetails();
	exit(1);
}




void MoleculeType::addReactionClass(ReactionClass * r, int rPosition)
{
	this->reactions.push_back(r);
	this->reactionPositions.push_back(rPosition);

	//We also have to check to make sure that if the reaction is a DOR reaction,
	//we remember it so we can updated it
	if(r->getRxnType()==ReactionClass::DOR_RXN) {
		if( r->getDORreactantPosition()==rPosition) {
			indexOfDORrxns.push_back(reactions.size()-1);
		}
	}
}



void MoleculeType::populateWithDefaultMolecules(int moleculeCount)
{
	if(DEBUG) cout<< " Populating "<< this->name << " with " << moleculeCount << " molecule(s)";
	if(DEBUG) cout<< " for a total of " << mList->size()+moleculeCount << " molecule(s)."<<endl;
	//mInstances.reserve(mInstances.size()+moleculeCount);
	for(int m=0; m<moleculeCount; m++)
	{
		if(DEBUG) cout<<" ("<<m+1<<") ";

		//Create the molecule (which knows how many components to make)
		this->genDefaultMolecule();
		//new Molecule(this);

		//Add the molecule to the list of molecules so we save it (does this automatically now!!!! )
		//mInstances.push_back(mol);
	}
}



void MoleculeType::setUpLocalFunctionListForMolecules()
{
	Molecule *mol;
	for(int m=0; m<mList->size(); m++ )
	{
	  	mol = mList->at(m);
	  	mol->setUpLocalFunctionList();
	}
}

void MoleculeType::prepareForSimulation()
{
	//cout<<"Preparing: "<<name<<endl;
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
  		//cout<<"start here"<<endl;
  		//First prepare the molecule for simulation
  		mol = mList->at(m);
  		mol->prepareForSimulation();

  		//Check each observable and see if this molecule should be counted
  		this->addToObservables(mol);
  		//cout<<"got here"<<endl;

  		//Check each reaction and add this molecule as a reactant if we have to
		for(rxnIter = reactions.begin(), r=0; rxnIter != reactions.end(); rxnIter++, r++ )
		{
			(*rxnIter)->tryToAdd(mol, reactionPositions.at(r));
  		}
	}
}

void MoleculeType::updateRxnMembership(Molecule * m)
{
	unsigned int r=0; ReactionClass *rxn;
	for(; r<reactions.size(); r++ )
	{
		rxn=reactions.at(r);
		double oldA = rxn->get_a();
		rxn->tryToAdd(m, reactionPositions.at(r));
		this->system->update_A_tot(oldA,rxn->update_a());
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
	int o=0;
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++ ){
  		//Only subtract if m happened to be an observable... this saves us a compare call
  		if(m->isObs(o)) { //(*obsIter)->isObservable(m)){
			(*obsIter)->subtract();
  		}
  		o++;
	}
}

void MoleculeType::removeFromRxns(Molecule * m)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		double oldA = (*rxnIter)->get_a();
		(*rxnIter)->remove(m, reactionPositions.at(r));
		this->system->update_A_tot(oldA,(*rxnIter)->update_a());
  	}
}




//TypeI local function: this molecule type depends on the value of this
//evaluated function
int MoleculeType::addLocalFunc_TypeI(LocalFunction *lf) {
	cout<<"adding type I local function to "<<name<<endl;
	locFuncs_typeI.push_back(lf);
	return locFuncs_typeI.size()-1;

}

//TypeII local function: this molecule type, when updated, changes the
//value of this function
int MoleculeType::addLocalFunc_TypeII(LocalFunction *lf) {
	cout<<"adding type II local functions to "<<name<<endl;
	locFuncs_typeII.push_back(lf);
	return locFuncs_typeII.size()-1;

}
















void MoleculeType::addAllToObservables()
{
	//Check each observable and see if this molecule should be counted
	Molecule *mol; int o=0;
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++){
  		(*obsIter)->clear();
  		for( int m=0; m<mList->size(); m++ )
  		{
  			mol = mList->at(m);
  		  	if((*obsIter)->isObservable(mol)) {
  		  		(*obsIter)->straightAdd();
  		  		mol->setIsObs(o,true);
  		  	} else {
  		  		mol->setIsObs(o,false);
  		  	}
  		}
  		o++;
	}

}



void MoleculeType::addToObservables(Molecule *m)
{
	//Check each observable and see if this molecule should be counted
	int o=0;
  	for(obsIter = observables.begin(); obsIter != observables.end(); obsIter++){
		//cout<<"Comparing(in add: "<<endl;
		//m->printDetails();
		if((*obsIter)->isObservable(m)) {
			(*obsIter)->add();
			m->setIsObs(o,true);
		} else {
			m->setIsObs(o,false);
		}
		o++;
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


void MoleculeType::printAllMolecules()
{
	for( int m=0; m<mList->size(); m++ ) {
		mList->at(m)->printDetails();
	}

}


void MoleculeType::printDetails() const
{
	cout<<"Molecule Type: "<< name << " type ID: " << type_id <<endl;

	cout<<"   -components ( ";
	for(int c=0; c<numOfComponents; c++) {

		cout<<compName[c];
		if(!isIntegerCompState[c]) {
			for(unsigned int s=0; s<possibleCompStates.at(c).size(); s++) {
				cout<<"~"<<possibleCompStates.at(c).at(s);
			}
		} else {
			cout<<"~integer[0-"<<possibleCompStates.at(c).at(possibleCompStates.at(c).size()-1)<<"]";
		}
		if(c<(numOfComponents-1)) cout<<", ";
	}
	cout<<" )"<<endl;

	//Output the local functions...
	cout<<"  Type I local functions include:";
	if(locFuncs_typeI.size()==0) cout<<"  none.";
	for(unsigned int ti=0; ti<locFuncs_typeI.size(); ti++) {
		cout<<"  "<<locFuncs_typeI.at(ti)->getNiceName();
	} cout<<endl;
	cout<<"  Type II local functions include:";
	if(locFuncs_typeII.size()==0) cout<<"  none.";
	for(unsigned int tii=0; tii<locFuncs_typeII.size(); tii++) {
		cout<<"  "<<locFuncs_typeII.at(tii)->getNiceName();
	} cout<<endl;


	cout<<"   -has "<< mList->size() <<" molecules."<<endl;
	cout<<"   -has "<< reactions.size() <<" reactions"<<endl;
//	cout<<"        of which "<< indexOfDORrxns.size() <<" are DOR rxns. "<<endl;
	cout<<"   -has "<< observables.size() <<" observables " <<endl;
}





