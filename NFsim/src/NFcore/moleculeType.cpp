#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;




MoleculeType::MoleculeType(
	string name,
	vector <string> &compName,
	System *s)
 : population_type( false )
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
 : population_type( false )
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
 : population_type( false )
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
 : population_type( false )
{
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}


MoleculeType::MoleculeType(
		string name,
		vector <string> &compName,
		vector <string> &defaultCompState,
		vector < vector<string> > &possibleCompStates,
		vector <bool> isIntegerComponent,
		bool pop_type,
		System *system)
 : population_type( pop_type )
{
	init(name, compName, defaultCompState, possibleCompStates, isIntegerComponent, system);
}


string MoleculeType::getComponentName(int cIndex) const {
	if(cIndex>=this->numOfComponents) {
		cerr<<"in MoleculeType:getComponentName, can't get that component!"<<endl;
		cerr<<"looking for cIndex: "<<cIndex<<endl;
		cerr<<"I only have: "<<numOfComponents<<" components."<<endl;
		exit(2);
	}
	return this->compName[cIndex];
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


	mList = new MoleculeList(this,2,system->getGlobalMoleculeLimit());
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
	delete [] eqCompOriginalName;




	//Delete all template molecules of this type that exist
	TemplateMolecule *t;
	while(allTemplates.size()>0)
	{
		t = allTemplates.back();
		allTemplates.pop_back();
		delete t;
	}




	delete mList;
}

void MoleculeType::addEquivalentComponents(vector <vector <string> > &identicalComponents)
{
	this->n_eqComp = identicalComponents.size();
	eqCompOriginalName = new string [n_eqComp];
	eqCompName=new string * [n_eqComp];
	eqCompIndex=new int *[n_eqComp];
	eqCompSizes=new int [n_eqComp];

	for(int i=0; i<n_eqComp; i++) {
		eqCompSizes[i]=identicalComponents.at(i).size();
		eqCompName[i] = new string [eqCompSizes[i]];
		eqCompIndex[i] = new int [eqCompSizes[i]];
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(k==0) {
				string tempString = identicalComponents.at(i).at(k);
				eqCompOriginalName[i] = tempString.substr(0,tempString.size()-1);
			}
			eqCompName[i][k] = identicalComponents.at(i).at(k);
			eqCompIndex[i][k] = getCompIndexFromName(eqCompName[i][k]);
		}
	}
}


bool MoleculeType::isIntegerComponent(string cName) const {
	for(int c=0; c<numOfComponents; c++)
			if(compName[c].compare(cName)==0) {
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
		if(eqCompOriginalName[i].compare(cName)==0) {
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

int MoleculeType::getEquivalenceClassNumber(int cIndex) const {
	for(int i=0; i<n_eqComp; i++) {
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(eqCompIndex[i][k]==cIndex) {
				return i;
			}
		}
	}
	return -1;
}

string MoleculeType::getEquivalenceClassComponentNameFromComponentIndex(int cIndex) const {
	for(int i=0; i<n_eqComp; i++) {
		for(int k=0; k<eqCompSizes[i]; k++) {
			if(eqCompIndex[i][k]==cIndex) {
				return eqCompOriginalName[i];
			}
		}
	}
	cerr<<"Could not find equivalency class component string for component number: "<<cIndex<<"!!!"<<endl;
	exit(1);
}

void MoleculeType::getEquivalencyClass(int *&components, int &n_components, string cName) const {
	for(int i=0; i<n_eqComp; i++) {
		if(eqCompOriginalName[i].compare(cName)==0) {
			components = eqCompIndex[i];
			n_components=eqCompSizes[i];
			return;
		}
	}
}
int MoleculeType::getEquivalencyClassNumber(string cName) const {
	for(int i=0; i<n_eqComp; i++) {
		if(eqCompOriginalName[i].compare(cName)==0) {
				return i;
		}
	}
	cerr<<"Could not find equivalency class number for component named: "<<cName<<"!!!"<<endl;
	exit(1);
}




string MoleculeType::getComponentStateName(int cIndex, int cValue) {
	if(cValue==Molecule::NOSTATE) return "NO_STATE";
	if(cValue>possibleCompStates.size() || cValue<0)
		return "?-value was: "+NFutil::toString(cValue);
	return possibleCompStates.at(cIndex).at(cValue);
}




Molecule *MoleculeType::genDefaultMolecule()
{
	Molecule *m;
	mList->create(m);
	m->setAlive(true);
	//cout<<"adding molecule: "<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<endl;

	return m;
}


void MoleculeType::addMoleculeToRunningSystem(Molecule *&mol)
{
	//cout<<"adding molecule: "<<mol->getMoleculeTypeName()<<"_"<<mol->getUniqueID()<<endl;
	//First prepare the molecule for simulation
	mol->setUpLocalFunctionList();
	mol->prepareForSimulation();
	mol->setAlive(true);

	mol->addToObservables();
	this->updateRxnMembership(mol);
}


void MoleculeType::addMoleculeToRunningSystemButDontUpdate(Molecule *&mol)
{
	//First prepare the molecule for simulation

	mol->setUpLocalFunctionList();
	mol->prepareForSimulation();
	mol->setAlive(true);

	//We assume observables and reaction membership will be updated later
	// (this is now the case for reaction firing)
}


void MoleculeType::removeMoleculeFromRunningSystem(Molecule *&m)
{
	if (system->isUsingComplex())
		// Need to manually unset canonical flag since we're not calling a Complex method
		m->getComplex()->unsetCanonical();

	//Remove this guy from the list, the observables list, and from all rxns
	mList->remove(m->getMolListId(), m);
	removeFromObservables(m);
	removeFromRxns(m);


	//We also have to remove all bonds
	for(int c=0; c<getNumOfComponents(); c++) {
		if(m->isBindingSiteBonded(c)) {
			Molecule::unbind(m,c);
		}
	}

	m->setAlive(false);

}


void MoleculeType::removeMoleculeFromRunningSystemButDontUpdate(Molecule *&m)
{
	//Remove this guy from the list, the observables list, and from all rxns
	mList->remove(m->getMolListId(), m);
	//removeFromObservables(m);
	//removeFromRxns(m);

	//We also have to remove all bonds
	for(int c=0; c<getNumOfComponents(); c++) {
		if(m->isBindingSiteBonded(c)) {
			Molecule::unbind(m,c);
		}
	}

	m->setAlive(false);
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

string MoleculeType::getMolObsName(int obsIndex) const {
	return molObs.at(obsIndex)->getName();
}

int MoleculeType::getMolObsCount(int obsIndex) const {
	return molObs.at(obsIndex)->getCount();
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
	else if(r->getRxnType()==ReactionClass::DOR2_RXN) {
		if( r->getDORreactantPosition()==rPosition) {
			indexOfDORrxns.push_back(reactions.size()-1);
		}
		else if( r->getDORreactantPosition2()==rPosition) {
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
  		//First prepare the molecule for simulation
  		mol = mList->at(m);
  		mol->prepareForSimulation();

  		//Check each observable and see if this molecule should be counted
  		this->addToObservables(mol);

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
		//cout<<"\n\n\n\n*****************"<<endl;
		//cout<<"looking at rxn: "<<rxn->getName()<<endl;
		rxn->tryToAdd(m, reactionPositions.at(r));
		this->system->update_A_tot(rxn,oldA,rxn->update_a());
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
	//cout<<"removing from observables:"<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<endl;
	//m->printDetails();

	//Check each observable and see if this molecule was counted, and if so, remove
	int ind=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
  	{
  		//Only subtract if m happened to be an observable... this saves us a compare call
  		//int matches = (*molObsIter)->isObservable(m);

  		// How many times does this observable match the molecule?
  		int matches = m->isObs(ind);
  		// subtract matches from observable
  		(*molObsIter)->subtract(matches);
		// set IsObs to zero, so if remove is called twice, we don't get negative counts.
  		m->setIsObs(ind,0);

  		ind++;
	}
}

void MoleculeType::removeFromRxns(Molecule * m)
{
	int r=0;
	for(rxnIter = reactions.begin(); rxnIter != reactions.end(); rxnIter++, r++ )
	{
		double oldA = (*rxnIter)->get_a();
		(*rxnIter)->remove(m, reactionPositions.at(r));
		this->system->update_A_tot((*rxnIter),oldA,(*rxnIter)->update_a());
  	}
}




//TypeI local function: this molecule type depends on the value of this
//evaluated function
int MoleculeType::addLocalFunc_TypeI(LocalFunction *lf) {
	locFuncs_typeI.push_back(lf);
	return locFuncs_typeI.size()-1;

}

//TypeII local function: this molecule type, when updated, changes the
//value of this function
int MoleculeType::addLocalFunc_TypeII(LocalFunction *lf) {
	locFuncs_typeII.push_back(lf);
	return locFuncs_typeII.size()-1;

}
















void MoleculeType::addAllToObservables()
{
	/////  WARNING:: when calling this function, be sure to clear all observables
	/////  first, because this function will not clear observables.

//	cout<<"+++++++++ "<<this->getName()<<endl;

	//Check each observable and see if this molecule should be counted
	Molecule *mol;  int o=0;  int matches=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
  	{
  		//cout<<"comparing to obs: "<<(*molObsIter)->getName()<<endl;

  		for( int m=0; m<mList->size(); m++ )
  		{
  			mol = mList->at(m);
  			matches = (*molObsIter)->isObservable(mol);
  			(*molObsIter)->add(matches);
  			mol->setIsObs(o,matches);
  			//cout<<"matches:"<<matches<<endl;
  		}
  		o++;
	}

}




void MoleculeType::addToObservables(Molecule *m)
{
	//Check each observable and see if this molecule should be counted
	int o=0;
  	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
  	{
		//cout<<"Comparing(in add: ";
		//cout<<m->getUniqueID()<<")"<<endl;

		int matches = (*molObsIter)->isObservable(m);
		m->setIsObs(o,matches);

		(*molObsIter)->add(matches);
		o++;
	}
	
}


void MoleculeType::outputMolObsNames(NFstream &fout)
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		fout<<"\t"<<(*molObsIter)->getName();
}

void MoleculeType::outputMolObsCounts(NFstream &fout)
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		fout<<"\t"<<(*molObsIter)->getCount();
}

void MoleculeType::printMolObsNames()
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++)
		cout<<"\t"<<(*molObsIter)->getName();
}

void MoleculeType::printMolObsCounts()
{
	for(molObsIter = molObs.begin(); molObsIter != molObs.end(); molObsIter++ )
		cout<<"\t"<<(*molObsIter)->getCount();
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
	cout<<"   -has "<< molObs.size() <<" molecules observables " <<endl;
}


// friend functions
template<class T>
NFstream& operator<<(NFstream& nfstream, const T& value)
{
    if (nfstream.useFile_)
	nfstream.file_ << value;
    else
	nfstream.str_ << value;

    return nfstream;
}



