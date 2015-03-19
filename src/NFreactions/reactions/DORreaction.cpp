


#include "reaction.hh"

#define DEBUG_MESSAGE 0


using namespace std;
using namespace NFcore;


//should also accept list of local functions and list of PointerNames for each of the functions...
DORRxnClass::DORRxnClass(
		string name,
		double baseRate,
		string baseRateName,
		TransformationSet *transformationSet,
		CompositeFunction *function,
		vector <string> &lfArgumentPointerNameList, System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
//	cout<<"ok, here we go..."<<endl;
	vector <TemplateMolecule *> dorMolecules;

	//////////////////////////////////////////////////////////////////////////////////////////
	//Step 1: Find the DOR reactant, and make sure there is only one.  DOR reactants
	//can be found because they have a LocalFunctionPointer Transformation that keeps
	//information about the pointer onto either a reactant species or a particular molecule
	//in the pattern.
	this->DORreactantIndex = -1;
	for(int r=0; (unsigned)r<n_reactants; r++) {
		for(int i=0; i<transformationSet->getNumOfTransformations(r); i++) {
			Transformation *transform = transformationSet->getTransformation(r,i);
//			cout<<"found transformation of type: "<<transform->getType()<<" for reactant: "<<r<<endl;
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {

				if(DORreactantIndex==-1)
				{
					if ( transformationSet->getTemplateMolecule(r)->getMoleculeType()->isPopulationType() )
					{   // DOR reactant is a population!
						cout<<"Error when creating DORRxnClass: "<<name<<endl;
						cout<<"DOR reactant cannot be a population type."<<endl;
						exit(1);
					}

					DORreactantIndex=r;
				}
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

	if(DEBUG_MESSAGE)cout<<"I determined that the DOR reactant is in fact: "<<DORreactantIndex<<endl;
	if(DEBUG_MESSAGE)cout<<"N_reactants: "<<transformationSet->getNreactants()<<endl;

	//////////////////////////////////////////////////////////////////////////////////////////
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
	if((unsigned)function->getNumOfArgs()!=lfArgumentPointerNameList.size()) {
		cout<<"Error when creating DORRxnClass: "<<name<<endl;
		cout<<"Number of arguments and LocalFunctionArgumentPointerNameList size do not match!"<<endl;
		exit(1);
	}


	//
	this->n_argMolecules=lfArgumentPointerNameList.size();
	argIndexIntoMappingSet =  new int [n_argMolecules];
	argMappedMolecule = new Molecule *[n_argMolecules];
	argScope = new int [n_argMolecules];


	for(int i=0; i<(int)lfArgumentPointerNameList.size(); i++) {
//		cout<<"Received local function arg: "<< lfArgumentPointerNameList.at(i)<<endl;
//		cout<<" Takes as argument this thang: "<< lfArgumentPointerNameList.at(i)<<endl;

		//Now search for the function argument...
		bool match = false;
		for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex); k++) {
			Transformation *transform = transformationSet->getTransformation(DORreactantIndex,k);
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				LocalFunctionReference *lfr = static_cast<LocalFunctionReference*>(transform);
				if(lfr->getPointerName()==lfArgumentPointerNameList.at(i)) {
//					cout<<"Found a match here!"<<endl;
					//cout<<"found scope should be: "<<lfr->getFunctionScope()<<endl;
					//If we got here, we found a match, so remember the index of the transformation
					//so we can quickly get the value of the function for any mapping object we try
					//to push on the reactant Tree.

					argIndexIntoMappingSet[i] =  k;
					argMappedMolecule[i] = 0;
					argScope[i] = lfr->getFunctionScope();


					//this->lfList.push_back(lfList.at(i));
					//localFunctionValue.push_back(0);
					//indexIntoMappingSet.push_back(k);
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


	delete [] hasMatched;

	//////////////////////////////////////////////////////////////////////////////////////////
	///  Step 3: Wheh! now we can finally get on the business of creating the reactant lists
	///  and the reactant tree and setting the usual reactionClass parameters

	//Remember that we are a DOR ReactionClass
	this->reactionType = ReactionClass::DOR_RXN;

	//Set up the reactant tree
	//reactantTree = new ReactantTree(this->DORreactantIndex,transformationSet,4);
	reactantTree = new ReactantTree(this->DORreactantIndex,transformationSet,32);

	//Set up the reactantLists
	reactantLists = new ReactantList *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if((signed)r!=this->DORreactantIndex)
			reactantLists[r]=(new ReactantList(r,transformationSet,25));
	}

	//Initialize a to zero
	this->a=0;


	//Set the actual function
	this->cf = function;

	//Add type I molecule dependencies, so that when this function
	//is reevaluated on a molecule, the molecule knows to update this reaction
	//for(unsigned int r=0; r<n_reactants; r++) {
	//	cf->addTypeIMoleculeDependency(this->reactantTemplates[r]->getMoleculeType());
	//}
	// TODO: determine if it's sufficient to only add the DORreactantIndex
	cf->addTypeIMoleculeDependency( reactantTemplates[DORreactantIndex]->getMoleculeType() );

}
DORRxnClass::~DORRxnClass() {

	for(unsigned int r=0; r<n_reactants; r++) {
		if(this->DORreactantIndex!=r)
			delete reactantLists[r];
	}

	delete [] reactantLists;
	delete reactantTree;

	delete [] argIndexIntoMappingSet;
	delete [] argMappedMolecule;
	delete [] argScope;

}

void DORRxnClass::init() {

	//Here we have to tell the molecules that they are part of this function
	//and for single molecule functions, we have to tell them also that they are in
	//this function, so they need to update thier value should they be transformed
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}



void DORRxnClass::remove(Molecule *m, unsigned int reactantPos)
{
	//cout<<"removing from a DOR!!"<<endl;
	if(reactantPos==(unsigned)this->DORreactantIndex) {
		//if(DEBUG_MESSAGE)cout<<" ... as a DOR"<<endl;

		// handle the DOR reactant
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			//cout<<"was in the tree, so we should remove"<<endl;
			reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	} else {

		// handle it normally...
		//if(DEBUG_MESSAGE)cout<<" ... as a normal reactant"<<endl;
		ReactantList *rl = reactantLists[reactantPos];
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}
	//if(DEBUG_MESSAGE)cout<<"finished removing"<<endl;
}


bool DORRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {

	//if(DEBUG_MESSAGE)cout<<endl<<endl<<"adding molecule to DORRxnClass"<<endl;
	//if(DEBUG_MESSAGE)m->printDetails();
	if(reactantPos==(unsigned)this->DORreactantIndex) {
	//	if(DEBUG_MESSAGE)cout<<" ... as a DOR"<<endl;
		//cout<<"RxnListMappingId: "<<m->getRxnListMappingId(m->getMoleculeType()->getRxnIndex(this,reactantPos))<<endl;

		// handle the DOR reactant
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

		//cout<<"trying to add to the tree:"<<endl;
		//m->printDetails();

		if(reactantTree->getHasClonedMappings()) {
			if(m->getRxnListMappingId(rxnIndex)>=0) {
				//cout<<"removing"<<endl;
				reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		}

		if(m->getRxnListMappingId(rxnIndex)>=0) {
			//if(DEBUG_MESSAGE)cout<<"was in the tree, so checking if we should remove"<<endl;
			if(!reactantTemplates[reactantPos]->compare(m)) {
				//if(DEBUG_MESSAGE)cout<<"removing..."<<endl;
				reactantTree->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			} else {}
		} else {
			//if(DEBUG_MESSAGE)cout<<"wasn't in the tree, so trying to push and compare"<<endl;
			ms=reactantTree->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,reactantTree,ms)) {
				//cout<<"shouldn't be in the tree, so we pop"<<endl;
				reactantTree->removeMappingSet(ms->getId());
			} else {
				//cout<<"should be in the tree, so confirm push."<<endl;
				//m->printDetails();
				//ms->printDetails();
				//we are keeping it, so evaluate the function and confirm the push
				double localFunctionValue = this->evaluateLocalFunctions(ms);
				//if(DEBUG_MESSAGE)cout<<"local function value is: "<<localFunctionValue<<endl;
				reactantTree->confirmPush(ms->getId(),localFunctionValue);
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
	} else {

		//Get the specified reactantList
		ReactantList *rl = reactantLists[reactantPos];
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

		if(rl->getHasClonedMappings()) {
			if(m->getRxnListMappingId(rxnIndex)>=0) {
				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		}

		//Here we get the standard update...
		if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
		{
			if(!reactantTemplates[reactantPos]->compare(m)) {
				//cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}

		} else {
			//Try to map it!
			ms = rl->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,rl,ms)) {
				//we must remove, if we did not match.  This will also remove
				//everything that was cloned off of the mapping set
				rl->removeMappingSet(ms->getId());
			} else {
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}



//		// handle it normally...
//		//if(DEBUG_MESSAGE)cout<<" ... as a normal reactant"<<endl;
//		ReactantList *rl = reactantLists[reactantPos];
//		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
//		if(m->getRxnListMappingId(rxnIndex)>=0) {
//			if(!reactantTemplates[reactantPos]->compare(m)) {
//				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
//				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
//			}
//		} else {
//			//try to map it.
//			ms = rl->pushNextAvailableMappingSet();
//			if(!reactantTemplates[reactantPos]->compare(m,rl,ms)) {
//				rl->popLastMappingSet();
//				//we just pushed, then popped, so molecule has not changed...
//			} else {
//				m->setRxnListMappingId(rxnIndex,ms->getId());
//			}
//		}
	}
	//if(DEBUG_MESSAGE)cout<<"finished adding"<<endl;
	return true;
}


int DORRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	if(reactantIndex==(unsigned)this->DORreactantIndex) {
		return reactantTree->size();
	}
	return isPopulationType[reactantIndex] ?
		       reactantLists[reactantIndex]->getPopulation()
	         : reactantLists[reactantIndex]->size();
}


int DORRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	if(reactantIndex==(unsigned)this->DORreactantIndex) {
		return reactantTree->size();
	}
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}


//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double DORRxnClass::evaluateLocalFunctions(MappingSet *ms)
{
	//Go through each function, and set the value of the function
	//this->argMappedMolecule
	//cout<<"\t\t\t\tDORRxnClass::evaluateLocalFunctions()"<<endl;
	//cout<<"dor is reevaluating its function."<<endl;

	//Grab the molecules needed for the local function to evaluate
	for(int i=0; i<this->n_argMolecules; i++) {
		//cout<<"here."<<endl;
		//cout<<"\t\t\t\t\t"<<i<<": argMappedMolecule="<<argMappedMolecule[i]<<" argIndexIntoMappingSet="<<argIndexIntoMappingSet[i]<<endl;
		this->argMappedMolecule[i] = ms->get(this->argIndexIntoMappingSet[i])->getMolecule();
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeType="<<argMappedMolecule[i]->getMoleculeTypeName()<<endl;
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeScope="<<argScope[i]<<endl;
	}

	//cout<<"done setting molecules, so know calling the composite function evaluate method."<<endl;
	int * reactantCounts = new int[this->n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if(r==this->DORreactantIndex) {
			reactantCounts[r]= reactantTree->size();
		}
		else {
			reactantCounts[r]=reactantLists[r]->size();
		}
	}

	double value = this->cf->evaluateOn(argMappedMolecule,argScope, reactantCounts, n_reactants);
	delete [] reactantCounts;
	//cout<<"\t\t\t\t\t"<<"composite function value="<<value<<endl;

	return value;

	/*Molecule

	for(int i=0; i<(signed)lfList.size(); i++) {
		Molecule *molObject = ms->get(this->indexIntoMappingSet.at(i))->getMolecule();
		int index = lfList.at(i)->getIndexOfTypeIFunctionValue(molObject);
		this->localFunctionValue.at(i)=molObject->getLocalFunctionValue(index);
		//cout<<"found that local function: "<<getName()<<" evaluates to: " <<localFunctionValue.at(i)<<endl;
	}
	return this->localFunctionValue.at(0);
	*/
}


double DORRxnClass::update_a() {
	a = baseRate;
	for(unsigned int i=0; i<n_reactants; i++) {
		if(i!=DORreactantIndex) {
			a*=(double)getCorrectedReactantCount(i);
		} else {
			a*=reactantTree->getRateFactorSum();
		}
	}
	return a;
}

void DORRxnClass::pickMappingSets(double randNumber) const
{
	//here we cannot just select a random molecule.  This is where all of our hard
	//work pays off.  We can use the tree to correctly select the next DOR reactant
	//(as well as all the other reactants.  So here we go...
	double rateFactorMultiplier = baseRate;
	for(unsigned int i=0; i<n_reactants; i++) {
		if(i!=(unsigned)DORreactantIndex) {
			if ( isPopulationType[i] ) {
				reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
			} else {
				reactantLists[i]->pickRandom(mappingSet[i]);
			}
			rateFactorMultiplier*=getReactantCount(i);
		}
	}

	if(randNumber<0) randNumber = NFutil::RANDOM(this->a);
	reactantTree->pickReactantFromValue(mappingSet[DORreactantIndex],randNumber,rateFactorMultiplier);

	//cout<<"tree size:        "<<reactantTree->size()<<endl;
	//cout<<"Choosing at tree: "<<mappingSet[DORreactantIndex]->getId()<<endl;
	//mappingSet[DORreactantIndex]->printDetails();
	//reactantTree->printDetails();
}

void DORRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) {
	if(reactantIndex==DORreactantIndex) {
		double newValue = evaluateLocalFunctions(reactantTree->getMappingSet(rxnListIndex));
		reactantTree->updateValue(rxnListIndex,newValue);
	} else {
		cout<<"Internal Error in DORRxnClass::notifyRateFactorChange!!  : trying to change a rate\n";
		cout<<"factor of a non-DOR reactant.  That means this function was called in error!\n";
		exit(1);
	}
}


void DORRxnClass::printDetails() const
{
	cout<<"DORRxnClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if(r!=(unsigned)DORreactantIndex) {
			cout<<"      -|"<< this->getReactantCount(r)<<" mappings|\t";
			cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
		} else {

			cout<<"      -(DOR) |"<< this->getReactantCount(r)<<" mappings|\t";
			cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
			cout<<"             (rateFactorSum="<<reactantTree->getRateFactorSum();
			cout<<")."<<endl;
		    //reactantTree->printDetails();
		}
	}

	//this->printFullDetails();

	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}







/*
 * DOR2RxnClass
 */

DOR2RxnClass::DOR2RxnClass(
		string name,
		double baseRate,
		string baseRateName,
		TransformationSet *transformationSet,
		CompositeFunction *function1,
		CompositeFunction *function2,
		vector <string> &lfArgumentPointerNameList1,
		vector <string> &lfArgumentPointerNameList2,
		System *s
	) : ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	// TODO: figure out if there are used for anything
	//vector <TemplateMolecule *> dorMolecules1;
	//vector <TemplateMolecule *> dorMolecules2;

	//////////////////////////////////////////////////////////////////////////////////////////
	//Step 1: Find the DOR reactants, and make sure there are exactly 2.  DOR reactants
	//can be found because they have a LocalFunctionPointer Transformation that keeps
	//information about the pointer onto either a reactant species or a particular molecule
	//in the pattern.
	DORreactantIndex1 = -1;
	DORreactantIndex2 = -1;
	for (int r=0; (unsigned)r<n_reactants; r++) {

		for (int i=0; i < transformationSet->getNumOfTransformations(r); i++) {

			Transformation *transform = transformationSet->getTransformation(r,i);
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {

				if (DORreactantIndex1 ==-1){
					if ( transformationSet->getTemplateMolecule(r)->getMoleculeType()->isPopulationType() )
					{   // DOR reactant is a population!
						cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
						cout<<"DOR reactant1 cannot be a population type."<<endl;
						exit(1);
					}

					DORreactantIndex1 = r;
				}
				else if (DORreactantIndex1 == r) {
					// this is ok
				}
				else if (DORreactantIndex2 ==-1) {
					if ( transformationSet->getTemplateMolecule(r)->getMoleculeType()->isPopulationType() ) {
						// DOR reactant is a population!
						cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
						cout<<"DOR reactant2 cannot be a population type."<<endl;
						exit(1);
					}

					DORreactantIndex2 = r;
				}
				else if (DORreactantIndex2 == r) {
					// this is ok
				}
				else {
					cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
					cout<<"DOR2 reactions only support two DOR reactants."<<endl;
					exit(1);
				}
			}
		}
	}

	if (DORreactantIndex1==-1 || DORreactantIndex2==-1) {
		cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
		cout<<"DOR2RxnClass requires pointers into two different reactant patterns, but fewer than 2 were found!"<<endl;
		exit(1);
	}

	if(DEBUG_MESSAGE)cout<<"I determined that the DOR reactant1 is in fact: "<<DORreactantIndex1<<endl;
	if(DEBUG_MESSAGE)cout<<"I determined that the DOR reactant2 is in fact: "<<DORreactantIndex2<<endl;
	if(DEBUG_MESSAGE)cout<<"N_reactants: "<<transformationSet->getNreactants()<<endl;


	//////////////////////////////////////////////////////////////////////////////////////////
	//Step 2: Some bookkeeping so that we can quickly get the function values from a mapping set
	// Now that we have found the DOR reactant, which can potentially have multiple functions, lets
	// figure out which functions apply to which
	// vector <int> indexIntoMappingSet;    //list of the index into the transformations for each of the local functions
	// vector <double> localFunctionValue;  //list of the value of each of the local functions needed to evaluate the rate law

	// DOR reactant1
	//Array to double check that we have used all pointer references we have created
	bool *hasMatched1 = new bool [transformationSet->getNumOfTransformations(DORreactantIndex1)];
	for (int i=0; i<transformationSet->getNumOfTransformations(DORreactantIndex1); i++) hasMatched1[i]=false;

	//make sure that we have the right number of functions and argument names
	if((unsigned)function1->getNumOfArgs()!=lfArgumentPointerNameList1.size()) {
		cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
		cout<<"Number of arguments in function1 and LocalFunctionArgumentPointerList1 size do not match!"<<endl;
		exit(1);
	}

	n_argMolecules1=lfArgumentPointerNameList1.size();
	argIndexIntoMappingSet1 =  new int [n_argMolecules1];
	argMappedMolecule1 = new Molecule *[n_argMolecules1];
	argScope1 = new int [n_argMolecules1];

	for(int i=0; i<(int)lfArgumentPointerNameList1.size(); i++) {
		//Now search for the function argument...
		bool match = false;
		for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex1); k++) {
			Transformation *transform = transformationSet->getTransformation(DORreactantIndex1,k);
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				LocalFunctionReference *lfr = static_cast<LocalFunctionReference*>(transform);
				if(lfr->getPointerName()==lfArgumentPointerNameList1.at(i)) {
					//If we got here, we found a match, so remember the index of the transformation
					//so we can quickly get the value of the function for any mapping object we try
					//to push on the reactant Tree.

					argIndexIntoMappingSet1[i] =  k;
					argMappedMolecule1[i] = 0;
					argScope1[i] = lfr->getFunctionScope();

					hasMatched1[k]=true;
					match=true;
				}
			}
		}
		if(!match){  //If there was no match found, then we've got issues...
			cout<<"Error when creating DOR2 reaction: "<<name<<endl;
			cout<<"Could not find a match in the templateMolecules for a pointer reference to species/molecule: ";
			cout<<lfArgumentPointerNameList1.at(i)<<endl;
			exit(1);
		}
	}

	//Just send out a warning if we didn't use one of the pointer references we were given
	for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex1); k++) {
		Transformation *transform = transformationSet->getTransformation(DORreactantIndex1,k);
		if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
			if(!hasMatched1[k]) {
				cout<<endl<<"Warning!  when creating DOR2RxnClass: "<<name<<endl;
				cout<<"Pointer reference: "<<  static_cast<LocalFunctionReference*>(transform)->getPointerName();
				cout<<" that was provided is not used in the local function definition."<<endl;
			}
		}
	}
    // done with DOR reactant 1
	delete [] hasMatched1;


	// DOR reactant2
	//Array to double check that we have used all pointer references we have created
	bool *hasMatched2 = new bool [transformationSet->getNumOfTransformations(DORreactantIndex2)];
	for (int i=0; i<transformationSet->getNumOfTransformations(DORreactantIndex1); i++) hasMatched2[i]=false;

	//make sure that we have the right number of functions and argument names
	if((unsigned)function2->getNumOfArgs()!=lfArgumentPointerNameList2.size()) {
		cout<<"Error when creating DOR2RxnClass: "<<name<<endl;
		cout<<"Number of arguments in function2 and LocalFunctionArgumentPointerList2 size do not match!"<<endl;
		exit(1);
	}

	n_argMolecules2=lfArgumentPointerNameList2.size();
	argIndexIntoMappingSet2 =  new int [n_argMolecules2];
	argMappedMolecule2 = new Molecule *[n_argMolecules2];
	argScope2 = new int [n_argMolecules2];

	for(int i=0; i<(int)lfArgumentPointerNameList2.size(); i++) {
		//Now search for the function argument...
		bool match = false;
		for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex2); k++) {
			Transformation *transform = transformationSet->getTransformation(DORreactantIndex2,k);
			if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
				LocalFunctionReference *lfr = static_cast<LocalFunctionReference*>(transform);
				if(lfr->getPointerName()==lfArgumentPointerNameList2.at(i)) {
					//If we got here, we found a match, so remember the index of the transformation
					//so we can quickly get the value of the function for any mapping object we try
					//to push on the reactant Tree.

					argIndexIntoMappingSet2[i] =  k;
					argMappedMolecule2[i] = 0;
					argScope2[i] = lfr->getFunctionScope();

					hasMatched2[k]=true;
					match=true;
				}
			}
		}
		if(!match){  //If there was no match found, then we've got issues...
			cout<<"Error when creating DOR2 reaction: "<<name<<endl;
			cout<<"Could not find a match in the templateMolecules for a pointer reference to species/molecule: ";
			cout<<lfArgumentPointerNameList2.at(i)<<endl;
			exit(1);
		}
	}



	//Just send out a warning if we didn't use one of the pointer references we were given
	for(int k=0; k<transformationSet->getNumOfTransformations(DORreactantIndex2); k++) {
		Transformation *transform = transformationSet->getTransformation(DORreactantIndex2,k);
		if((unsigned)transform->getType()==TransformationFactory::LOCAL_FUNCTION_REFERENCE) {
			if(!hasMatched2[k]) {
				cout<<endl<<"Warning!  when creating DOR2RxnClass: "<<name<<endl;
				cout<<"Pointer reference: "<<  static_cast<LocalFunctionReference*>(transform)->getPointerName();
				cout<<" that was provided is not used in the local function definition."<<endl;
			}
		}
	}
    // done with DOR reactant 2
	delete [] hasMatched2;


	//////////////////////////////////////////////////////////////////////////////////////////
	///  Step 3: Wheh! now we can finally get on the business of creating the reactant lists
	///  and the reactant tree and setting the usual reactionClass parameters

	//Remember that we are a DOR ReactionClass
	this->reactionType = ReactionClass::DOR2_RXN;

	//Set up the reactant trees
	reactantTree1 = new ReactantTree(this->DORreactantIndex1,transformationSet,32);
	reactantTree2 = new ReactantTree(this->DORreactantIndex2,transformationSet,32);

	//Set up the reactantLists
	reactantLists = new ReactantList *[n_reactants];
	for (unsigned int r=0; r<n_reactants; r++) {
		if( (signed)r!=this->DORreactantIndex1  &&  (signed)r!=this->DORreactantIndex2 )
			reactantLists[r]=(new ReactantList(r,transformationSet,25));
	}

	//Initialize a to zero
	this->a=0;

	//Set the actual functions
	this->cf1 = function1;
	this->cf2 = function2;

	// TODO: figure out if we really need to add all these molecules as TypeI dependencies.
	//  It seems like we really only need to do this for the DOR reactant

	//Add type I molecule dependencies, so that when this function
	//is reevaluated on a molecule, the molecule knows to update this reaction
	//for (unsigned int r=0; r<n_reactants; r++) {
	//	cf1->addTypeIMoleculeDependency( reactantTemplates[r]->getMoleculeType() );
	//	cf2->addTypeIMoleculeDependency( reactantTemplates[r]->getMoleculeType() );
	//}
	cf1->addTypeIMoleculeDependency( reactantTemplates[DORreactantIndex1]->getMoleculeType() );
	cf2->addTypeIMoleculeDependency( reactantTemplates[DORreactantIndex2]->getMoleculeType() );

}


DOR2RxnClass::~DOR2RxnClass() {

	for(unsigned int r=0; r<n_reactants; r++) {
		if( r != DORreactantIndex1  && 	r != DORreactantIndex2 )
			delete reactantLists[r];
	}

	delete [] reactantLists;

	delete reactantTree1;
	delete reactantTree2;

	delete [] argIndexIntoMappingSet1;
	delete [] argIndexIntoMappingSet2;
	delete [] argMappedMolecule1;
	delete [] argMappedMolecule2;
	delete [] argScope1;
	delete [] argScope2;
}


void DOR2RxnClass::init() {

	//Here we have to tell the molecules that they are part of this function
	//and for single molecule functions, we have to tell them also that they are in
	//this function, so they need to update thier value should they be transformed
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void DOR2RxnClass::remove(Molecule *m, unsigned int reactantPos)
{
	// removing molecule from a DOR!!
	if(reactantPos==(unsigned)this->DORreactantIndex1){
		// handle the DOR reactant1
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			//cout<<"was in the tree, so we should remove"<<endl;
			reactantTree1->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}
	else if (reactantPos==(unsigned)this->DORreactantIndex2){
		// handle the DOR reactant2
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			//cout<<"was in the tree, so we should remove"<<endl;
			reactantTree2->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}
	else {
		// handle it normally...
		ReactantList *rl = reactantLists[reactantPos];
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
		if(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}
}


bool DOR2RxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {

	// adding molecule to DOR2RxnClass
	//if(DEBUG_MESSAGE)m->printDetails();
	if (reactantPos==(unsigned)this->DORreactantIndex1) {

		// handle the DOR reactant
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

		if(reactantTree1->getHasClonedMappings()) {
			if(m->getRxnListMappingId(rxnIndex)>=0) {
				reactantTree1->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		}

		if(m->getRxnListMappingId(rxnIndex)>=0) {
			// was in the tree, so checking if we should remove
			if(!reactantTemplates[reactantPos]->compare(m)) {
				// removing
				reactantTree1->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			} else {}
		} else {
			// wasn't in the tree, so trying to push and compare
			ms=reactantTree1->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,reactantTree1,ms)) {
				//cout<<"shouldn't be in the tree, so we pop"<<endl;
				reactantTree1->removeMappingSet(ms->getId());
			} else {
				//we are keeping it, so evaluate the function and confirm the push
				double localFunctionValue = evaluateLocalFunctions1(ms);
				//if(DEBUG_MESSAGE)cout<<"local function value is: "<<localFunctionValue<<endl;
				reactantTree1->confirmPush(ms->getId(),localFunctionValue);
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
	}
	else if (reactantPos==(unsigned)this->DORreactantIndex2) {

		// handle the DOR reactant
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

		if(reactantTree2->getHasClonedMappings()) {
			if(m->getRxnListMappingId(rxnIndex)>=0) {
				reactantTree2->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		}

		if(m->getRxnListMappingId(rxnIndex)>=0) {
			// was in the tree, so checking if we should remove
			if(!reactantTemplates[reactantPos]->compare(m)) {
				// removing
				reactantTree2->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			} else {}
		} else {
			// wasn't in the tree, so trying to push and compare
			ms=reactantTree2->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,reactantTree2,ms)) {
				//cout<<"shouldn't be in the tree, so we pop"<<endl;
				reactantTree2->removeMappingSet(ms->getId());
			} else {
				//we are keeping it, so evaluate the function and confirm the push
				double localFunctionValue = this->evaluateLocalFunctions2(ms);
				//if(DEBUG_MESSAGE)cout<<"local function value is: "<<localFunctionValue<<endl;
				reactantTree2->confirmPush(ms->getId(),localFunctionValue);
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
	}
	else {
		//Get the specified reactantList
		ReactantList *rl = reactantLists[reactantPos];
		int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

		if(rl->getHasClonedMappings()) {
			if(m->getRxnListMappingId(rxnIndex)>=0) {
				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}
		}

		//Here we get the standard update...
		if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
		{
			if(!reactantTemplates[reactantPos]->compare(m)) {
				//cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
				rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
				m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
			}

		} else {
			//Try to map it!
			ms = rl->pushNextAvailableMappingSet();
			if(!reactantTemplates[reactantPos]->compare(m,rl,ms)) {
				//we must remove, if we did not match.  This will also remove
				//everything that was cloned off of the mapping set
				rl->removeMappingSet(ms->getId());
			} else {
				m->setRxnListMappingId(rxnIndex,ms->getId());
			}
		}
	}
	return true;
}


int DOR2RxnClass::getReactantCount(unsigned int reactantIndex) const
{
	if (reactantIndex==(unsigned)this->DORreactantIndex1) {
		return reactantTree1->size();
	}
	if (reactantIndex==(unsigned)this->DORreactantIndex2) {
		return reactantTree2->size();
	}
	return isPopulationType[reactantIndex] ?
		       reactantLists[reactantIndex]->getPopulation()
	         : reactantLists[reactantIndex]->size();
}


int DOR2RxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	if (reactantIndex==(unsigned)DORreactantIndex1) {
		return reactantTree1->size();
	}
	else if (reactantIndex==(unsigned)DORreactantIndex2) {
		return reactantTree1->size();
	}
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}



//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double DOR2RxnClass::evaluateLocalFunctions1(MappingSet *ms)
{
	//cout << "DOR2RxnClass::evaluateLocalFunctions1(" << ms << ")" << endl;
	//cout << "n_argMolecules1: " << n_argMolecules1 << endl;
	//cout << "argIndexIntoMappingSet1: " << argIndexIntoMappingSet1[0] << endl;
	//Go through each function, and set the value of the function

	//Grab the molecules needed for the local function to evaluate
	for(int i=0; i < n_argMolecules1; i++) {
		argMappedMolecule1[i] = ms->get(argIndexIntoMappingSet1[i])->getMolecule();
	}
	//cout << "argMappedMolecule1: " << argMappedMolecule1[0]->getMoleculeTypeName() << endl;
	//cout << "argScope1: " << argScope1[0] << endl;

	// done setting molecules, so now calling the composite function evaluate method
	int * reactantCounts = new int[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if(r==(unsigned int)DORreactantIndex1) {
			reactantCounts[r] = reactantTree1->size();
		}
		else if(r==(unsigned int)DORreactantIndex2) {
			reactantCounts[r] = reactantTree2->size();
		}
		else {
			reactantCounts[r] = reactantLists[r]->size();
		}
		//cout << "n_reactants[" << r << "]=" << reactantCounts[r] << endl;
	}

	double value = cf1->evaluateOn(argMappedMolecule1, argScope1, reactantCounts, n_reactants);
	//cout << "return value=" << value << endl;

	delete [] reactantCounts;
	return value;
}


//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double DOR2RxnClass::evaluateLocalFunctions2(MappingSet *ms)
{
	//cout << "DOR2RxnClass::evaluateLocalFunctions2(" << ms << ")" << endl;
	//cout << "mapping molecule type: " << ms->get(0)->getMolecule()->getMoleculeTypeName() << endl;
	//Go through each function, and set the value of the function

	//Grab the molecules needed for the local function to evaluate
	for (int i=0; i < n_argMolecules2; i++) {
		argMappedMolecule2[i] = ms->get(argIndexIntoMappingSet2[i])->getMolecule();
	}

	// done setting molecules, so now calling the composite function evaluate method
	int * reactantCounts = new int[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		if(r==DORreactantIndex1) {
			reactantCounts[r] = reactantTree1->size();
		}
		else if(r==this->DORreactantIndex2) {
			reactantCounts[r] = reactantTree2->size();
		}
		else {
			reactantCounts[r] = reactantLists[r]->size();
		}
	}

	double value = cf2->evaluateOn(argMappedMolecule2, argScope2, reactantCounts, n_reactants);
	//cout << "return value=" << value << endl;

	delete [] reactantCounts;
	return value;
}


double DOR2RxnClass::update_a() {
	a = baseRate;
	//cout << "> DOR2RxnClass::update_a()" << endl;
	//cout << "baseRate=" << baseRate << endl;
	for (unsigned int i=0; i<n_reactants; i++) {
		if (i==(unsigned int)DORreactantIndex1) {
			a*=reactantTree1->getRateFactorSum();
			//cout << i << ":rateFactorSum1=" << reactantTree1->getRateFactorSum() << endl;
		}
		else if (i==(unsigned int)DORreactantIndex2) {
			a*=reactantTree2->getRateFactorSum();
			//cout << i << ":rateFactorSum2=" << reactantTree2->getRateFactorSum() << endl;
		}
		else {
			a*=(double)getCorrectedReactantCount(i);
			//cout << i << ":ReactantCount=" << (double)getCorrectedReactantCount(i) << endl;
		}
	}
	//cout << "update_a=" << a << endl;
	return a;
}


void DOR2RxnClass::pickMappingSets(double randNumber) const
{
	//here we cannot just select a random molecule.  This is where all of our hard
	//work pays off.  We can use the tree to correctly select the next DOR reactant
	//(as well as all the other reactants.  So here we go...
	//double rateFactorMultiplier = baseRate;
	for(unsigned int i=0; i<n_reactants; i++) {
		if( i!=(unsigned)DORreactantIndex1 && i!=(unsigned)DORreactantIndex2) {
			if ( isPopulationType[i] ) {
				reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
			} else {
				reactantLists[i]->pickRandom(mappingSet[i]);
			}
			//rateFactorMultiplier*=getReactantCount(i);
		}
	}

	double randNumber1 = NFutil::RANDOM( reactantTree1->getRateFactorSum() );
	reactantTree1->pickReactantFromValue( mappingSet[DORreactantIndex1], randNumber1, 1.0);

	double randNumber2 = NFutil::RANDOM( reactantTree2->getRateFactorSum() );
	reactantTree2->pickReactantFromValue( mappingSet[DORreactantIndex2], randNumber2, 1.0);

}


void DOR2RxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) {
	if (reactantIndex==DORreactantIndex1) {
		double newValue = evaluateLocalFunctions1(reactantTree1->getMappingSet(rxnListIndex));
		reactantTree1->updateValue(rxnListIndex,newValue);
	}
	else if (reactantIndex==DORreactantIndex2) {
		double newValue = evaluateLocalFunctions2(reactantTree2->getMappingSet(rxnListIndex));
		reactantTree2->updateValue(rxnListIndex,newValue);
	}
	else {
		cout<<"Internal Error in DORRxnClass::notifyRateFactorChange!!  : trying to change a rate\n";
		cout<<"factor of a non-DOR reactant.  That means this function was called in error!\n";
		exit(1);
	}
}


void DOR2RxnClass::printDetails() const
{
	cout<<"DOR2RxnClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		if( r==(unsigned)DORreactantIndex1) {
			cout<<"      -(DOR1) |"<< getReactantCount(r)<<" mappings|\t";
			cout<<reactantTemplates[r]->getPatternString()<<"\n";
			cout<<"             (rateFactorSum="<<reactantTree1->getRateFactorSum();
			cout<<")."<<endl;
		}
		else if( r==(unsigned)DORreactantIndex2) {
			cout<<"      -(DOR2) |"<< getReactantCount(r)<<" mappings|\t";
			cout<<reactantTemplates[r]->getPatternString()<<"\n";
			cout<<"             (rateFactorSum="<<reactantTree2->getRateFactorSum();
			cout<<")."<<endl;
		} else {
			cout<<"      -|"<< getReactantCount(r)<<" mappings|\t";
			cout<<reactantTemplates[r]->getPatternString()<<"\n";

		}
	}

	if (n_reactants==0)
		cout << "      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}





