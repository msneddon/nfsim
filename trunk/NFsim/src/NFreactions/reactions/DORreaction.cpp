


#include "reaction.hh"


using namespace std;
using namespace NFcore;


//should also accept list of local functions and list of PointerNames for each of the functions...
DORRxnClass::DORRxnClass(
		string name,
		double baseRate,
		TransformationSet *transformationSet,
		vector <LocalFunction *> &lfList,
		vector <string> &lfArgumentPointerNameList) :
	ReactionClass(name,baseRate,transformationSet)
{
	cout<<"ok, here we go..."<<endl;
	vector <TemplateMolecule *> dorMolecules;

	//Step 1: Find the DOR reactant, and make sure there is only one.  DOR reactants
	//can be found because they have a LocalFunctionPointer Transformation that keeps
	//information about the pointer onto either a reactant species or a particular molecule
	//in the pattern.
	this->DORreactantIndex = -1;
	for(int r=0; (unsigned)r<n_reactants; r++) {
		for(int i=0; i<transformationSet->getNumOfTransformations(r); i++) {
			Transformation *transform = transformationSet->getTransformation(r,i);
			cout<<"found transformation of type: "<<transform->getType()<<" for reactant: "<<r<<endl;
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				if(DORreactantIndex==-1) { DORreactantIndex=r; }
				else if(DORreactantIndex!=r) {
					cout<<"Error when creating DORRxnClass: "<<name<<endl;
					cout<<"DOR reactions currently only support one DOR reactant.  This means that you can"<<endl;
					cout<<"only have a pointer to one or the other of the two reactants, but not both."<<endl;
					exit(1);
				}
			}
		}
	}
	if(DORreactantIndex==-1) {
		cout<<"Error when creating DORRxnClass: "<<name<<endl;
		cout<<"You don't have any pointers onto the Molecules or Species, so you can't have a local function!"<<endl;
		cout<<"That means that this is not a DOR reaction at all!"<<endl;
		exit(1);
	}

	cout<<"I determined that the DOR reactant is in fact: "<<DORreactantIndex<<endl;

	cout<<"N_reactants: "<<transformationSet->getNreactants();



	cout<<endl<<endl<<endl;
	//Step 2: Some bookkeeping so that we can quickly get the function values from a mapping set
	// Now that we have found the DOR reactant, which can potentially have multiple functions, lets
	// figure out which functions apply to which
	// vector <int> indexIntoMappingSet;  //list of the index into the transformations for each of the local functions
	//vector <double> localFunctionValue;  //list of the value of each of the local functions needed to evaluate
	                                    //the rate law
	//Array to double check that we have used all pointer references we have created
	bool *hasMatched = new bool [transformationSet->getNumOfTransformations(DORreactantIndex)];
	for(int i=0; i<transformationSet->getNumOfTransformations(DORreactantIndex); i++) hasMatched[i]=false;

	//make sure that we have the right number of functions and argument names
	if(lfList.size()!=lfArgumentPointerNameList.size()) {
		cout<<"Error when creating DORRxnClass: "<<name<<endl;
		cout<<"LocalFunctionList size and LocalFunctionArguementPointerNameList size do not match!"<<endl;
		exit(1);
	}
	for(int i=0; i<(signed)lfList.size(); i++) {
		cout<<"Received local function: "<< lfList.at(i)->getNiceName()<<endl;
		cout<<" Takes as argument this thang: "<< lfArgumentPointerNameList.at(i)<<endl;

		//Now search for the function argument...
		bool match = false;
		for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex); k++) {
			Transformation *transform = transformationSet->getTransformation(DORreactantIndex,k);
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				LocalFunctionReference *lfr = static_cast<LocalFunctionReference*>(transform);
				if(lfr->getPointerName()==lfArgumentPointerNameList.at(i)) {
					cout<<"Found a match here!"<<endl;
					//If we got here, we found a match, so remember the index of the transformation
					//so we can quickly get the value of the function for any mapping object we try
					//to push on the reactant Tree.
					this->lfList.push_back(lfList.at(i));
					localFunctionValue.push_back(0);
					indexIntoMappingSet.push_back(k);
					hasMatched[k]=true;
					match=true;
				}
			}
		}
		if(!match){  //If there was no match found, then we've got issues...
			cout<<"Error when creating DOR reaction: "<<name<<endl;
			cout<<"Could not find a match in the templateMolecules for the pointer reference to species/molecule: ";
			cout<<lfArgumentPointerNameList.at(i)<<endl;
			exit(1);
		}
	}

	//Just send out a warning if we didn't use one of the pointer references we were given
	for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex); k++) {
		Transformation *transform = transformationSet->getTransformation(DORreactantIndex,k);
		if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
			if(!hasMatched[k]) {
				cout<<endl<<"Warning!  when creating DORrxnClass: "<<name<<endl;
				cout<<"Pointer reference: "<<  static_cast<LocalFunctionReference*>(transform)->getPointerName();
				cout<<" that was provided is not used in the local function definition."<<endl;
	}	}	}


	///  Step 3: Wheh! now we can finally get on the business of creating the reactant lists
	///  and the reactant tree and setting the usual reactionClass parameters

	//Remeber that we are a DOR ReactionClass
	this->reactionType = ReactionClass::DOR_RXN;

	//Set up the reactant tree
	reactantTree = new ReactantTree(this->DORreactantIndex,transformationSet,32);
	//reactantTree = new ReactantTree(this->DORreactantIndex,transformationSet,32);

	//Set up the reactantLists
	reactantLists = new ReactantList *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if((signed)r!=this->DORreactantIndex)
			reactantLists[r]=(new ReactantList(r,transformationSet,25));
	}

	//Initialize a to zero
	this->a=0;


	// i think we be done now
}
DORRxnClass::~DORRxnClass() {


}

void DORRxnClass::init() {

	//cout<<"called init."<<endl; exit(1);
	//Here we have to tell the molecules that they are part of this function
	//and for single molecule functions, we have to tell them also that they are in
	//this function
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}


}


bool DORRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {


	cout<<endl<<endl<<"adding molecule to DORRxnClass"<<endl;
	m->printDetails();
	if(reactantPos==(unsigned)this->DORreactantIndex) {
		cout<<" ... as a DOR"<<endl;


		cout<<"RxnListMappingId: "<<m->getRxnListMappingId(m->getMoleculeType()->getRxnIndex(this,reactantPos))<<endl;

		// handle the DOR reactant
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			cout<<"was in the tree, so checking if we should remove"<<endl;
			if(!reactantTemplates[reactantPos]->compare(m)) {
				cout<<"removing..."<<endl;
				reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			} else { cout<<"not removing"<<endl; }
		} else {
			cout<<"wasn't in the tree, so trying to push and compare"<<endl;
			ms=reactantTree->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,ms)) {
				cout<<"shouldn't be in the tree, so we pop"<<endl;
				reactantTree->popLastMappingSet();
			} else {
				cout<<"should be in the tree, so confirm push."<<endl;
				//we are keeping it, so evaluate the function and confirm the push
				double localFunctionValue = this->evaluateLocalFunctions(ms);
				reactantTree->confirmPush(ms->getId(),localFunctionValue);
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}


		reactantTree->printDetails();
//		ReactantList * rl = new ReactantList(DORreactantIndex,transformationSet,25);
//		MappingSet *ms = rl->pushNextAvailableMappingSet();
//		bool match = this->reactantTemplates[DORreactantIndex]->compare(m,ms);
//		cout<<"first, did we match?:"<<match<<endl;
//		cout<<"now, let us evaluate!"<<endl;
//		evaluateLocalFunctions(ms);
	} else {

		// handle it normally...
		cout<<" ... as a normal reactant"<<endl;
		ReactantList *rl = reactantLists[reactantPos];
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			if(!reactantTemplates[reactantPos]->compare(m)) {
				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		} else {
			//try to map it.
			ms = rl->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,ms)) {
				rl->popLastMappingSet();
				//we just pushed, then popped, so molecule has not changed...
			} else {
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
	}
	cout<<"finished adding"<<endl;
	return true;
}


unsigned int DORRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	if(reactantIndex==(unsigned)this->DORreactantIndex) {
		return reactantTree->size();
	}
	reactantLists[reactantIndex]->size();
}

//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double DORRxnClass::evaluateLocalFunctions(MappingSet *ms)
{
	//Go through each function, and set the value of the function
	for(int i=0; i<(signed)lfList.size(); i++) {
		Molecule *molObject = ms->get(this->indexIntoMappingSet.at(i))->getMolecule();
		int index = lfList.at(i)->getIndexOfTypeIFunctionValue(molObject);
		this->localFunctionValue.at(i)=molObject->getLocalFunctionValue(index);
		//cout<<"found that local function: "<<getName()<<" evaluates to: " <<localFunctionValue.at(i)<<endl;
	}

}
