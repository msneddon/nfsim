

#include "NFreactions.hh"



//////////////////////////////////////////
////  BINDING
///////////////////////////////////////////
//ReactionSimpleBinding::ReactionSimpleBinding(	
//	char * name, 
//	int n_reactants, 
//	TemplateMolecule ** reactantTemplates, 
//	double rate,
//	int bSiteIndex1,
//	int bSiteIndex2): ReactionClass(name, 2, reactantTemplates, rate)
//{
//	this->bSiteIndex1 = bSiteIndex1;
//	this->bSiteIndex2 = bSiteIndex2;
//}
//
//ReactionSimpleBinding::ReactionSimpleBinding(	char * name, 
//	int n_reactants, 
//	TemplateMolecule ** reactantTemplates, 
//	double rate,
//	char* bSiteName1,
//	char* bSiteName2): ReactionClass(name, 2, reactantTemplates, rate)
//{
//	this->bSiteIndex1 = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bSiteName1);
//	this->bSiteIndex2 = reactantTemplates[1]->getMoleculeType()->getBindingSiteIndex(bSiteName2);
//}
//
//
//void ReactionSimpleBinding::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]==NULL || reactants[1]==NULL) {
//		cerr<<"Binding Rxn Failed!! a reactant was null!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	Molecule::bind(reactants[0],bSiteIndex1,reactants[1],bSiteIndex2);	
//}
//
//
//
//////////////////////////////////////////
////  UNBINDING
///////////////////////////////////////////
//ReactionUnbinding::ReactionUnbinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int bSiteIndex) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = bSiteIndex;								
//}
//
//ReactionUnbinding::ReactionUnbinding(	char * name, 
//	int n_reactants, 
//	TemplateMolecule ** reactantTemplates, 
//	double rate,
//	char * bSiteName) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bSiteName);
//}
//
//void ReactionUnbinding::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]->isBindingSiteOpen(bSiteIndex)) {
//		cerr<<"UnBinding Rxn Failed!! nothing to unbind!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	Molecule::unbind(reactants[0],bSiteIndex);
//}
//
//////////////////////////////////////////
////  Symmetric UNBINDING (same as binding, but rate is cut in half
///////////////////////////////////////////
//SymmetricReactionUnbinding::SymmetricReactionUnbinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int bSiteIndex) : ReactionClass(name, 1, reactantTemplates, (rate/2.0))
//{
//	this->bSiteIndex = bSiteIndex;								
//}
//void SymmetricReactionUnbinding::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]->isBindingSiteOpen(bSiteIndex)) {
//		cerr<<"UnBinding Rxn Failed!! nothing to unbind!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	Molecule::unbind(reactants[0],bSiteIndex);
//}
//
//
//
//
//////////////////////////////////////////
////  CHANGE STATE REACTION
///////////////////////////////////////////
//ReactionChangeState::ReactionChangeState(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int stateIndex,
//								int newStateValue) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->stateIndex = stateIndex;
//	this->newStateValue = newStateValue;
//}
//
//
//ReactionChangeState::ReactionChangeState(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName,
//								int newStateValue) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->stateIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName);
//	this->newStateValue = newStateValue;
//}
//
//
//void ReactionChangeState::transformReactants(Molecule ** reactants, int nReactants)
//{
//	reactants[0]->setState(stateIndex, newStateValue);
//}
//
//
//ReactionChangeTwoStates::ReactionChangeTwoStates(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->stateIndex1 = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName1);
//	this->newStateValue1 = newStateValue1;
//	this->stateIndex2 = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName2);
//	this->newStateValue2 = newStateValue2;
//}
//void ReactionChangeTwoStates::transformReactants(Molecule ** reactants, int nReactants)
//{
//	reactants[0]->setState(stateIndex1, newStateValue1);
//	reactants[0]->setState(stateIndex2, newStateValue2);
//	
//}
//
//
//
//
//////////////////////////////////////////
////  CHANGE STATE and UNBIND REACTION
///////////////////////////////////////////
//ReactionUnbindAndChangeState::ReactionUnbindAndChangeState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName,
//								int newStateValue) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bindingSiteName);
//	this->stateIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName);
//	this->newStateValue = newStateValue;
//}
//
//void ReactionUnbindAndChangeState::transformReactants(Molecule ** reactants, int nReactants)
//{
//	//First change state	
//	reactants[0]->setState(stateIndex, newStateValue);
//	
//	//Next, unbind
//	Molecule::unbind(reactants[0],bSiteIndex);	
//}
//
//
//ReactionUnbindAndChangeBothState::ReactionUnbindAndChangeBothState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bindingSiteName);
//	this->stateIndex1 = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName1);
//	this->newStateValue1 = newStateValue1;
//	
//	int tempBsiteIndex = reactantTemplates[0]->getTemplateBsiteIndexFromMoleculeBsiteIndex(bSiteIndex);
//	TemplateMolecule * t = reactantTemplates[0]->getBondedTemplateMolecule(tempBsiteIndex);
//	this->newStateValue2 = t->getMoleculeType()->getStateIndex(stateName2);
//	
//	this->newStateValue2 = newStateValue2;
//}
//
//
//void ReactionUnbindAndChangeBothState::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]->isBindingSiteOpen(bSiteIndex)) {
//		cerr<<"UnBinding & both state change Rxn Failed!! nothing to unbind!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	
//	//First change states	
//	reactants[0]->setState(stateIndex1, newStateValue1);
//	reactants[0]->getBondedMolecule(bSiteIndex)->setState(stateIndex2,newStateValue2);
//	//Next, unbind
//	Molecule::unbind(reactants[0],bSiteIndex);	
//}
//
//
//
//////////////////////////////////////////
////  DECREMENT STATE and UNBIND REACTION
///////////////////////////////////////////
//ReactionUnbindAndDecrementState::ReactionUnbindAndDecrementState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bindingSiteName);
//	this->stateIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName);
//	
//}
//
//void ReactionUnbindAndDecrementState::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]->isBindingSiteOpen(bSiteIndex)) {
//		cerr<<"UnBinding & decrement state Rxn Failed!! nothing to unbind!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	
//	//First decrement state	
//	int currentStateValue = reactants[0]->getState(stateIndex);
//	reactants[0]->setState(stateIndex, (currentStateValue-1));
//	
//	//Next, unbind
//	Molecule::unbind(reactants[0],bSiteIndex);	
//}
//
//
//
//////////////////////////////////////////
////  INCREMENT STATE and UNBIND REACTION
///////////////////////////////////////////
//ReactionUnbindAndIncrementState::ReactionUnbindAndIncrementState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName) : ReactionClass(name, 1, reactantTemplates, rate)
//{
//	this->bSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(bindingSiteName);
//	this->stateIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName);
//	
//}
//
//void ReactionUnbindAndIncrementState::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(reactants[0]->isBindingSiteOpen(bSiteIndex)) {
//		cerr<<"UnBinding & increment state Rxn Failed!! nothing to unbind!"<<endl;
//		this->printDetails();
//		exit(1);
//	}
//	
//	//First decrement state	
//	int currentStateValue = reactants[0]->getState(stateIndex);
//	reactants[0]->setState(stateIndex, (currentStateValue+1));
//	
//	//Next, unbind
//	Molecule::unbind(reactants[0],bSiteIndex);	
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//////////////////////////////////////////
////  CHANGE STATE OF TWO MOLECULES
///////////////////////////////////////////
//ReactionChangeStateOfTwoMolecules::ReactionChangeStateOfTwoMolecules(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2) : ReactionClass(name, 2, reactantTemplates, rate)
//{
//	this->stateIndex1 = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName1);
//	this->newStateValue1 = newStateValue1;
//	this->stateIndex2 = reactantTemplates[1]->getMoleculeType()->getStateIndex(stateName2);
//	this->newStateValue2 = newStateValue2;
//}
//								
//								
//								
//void ReactionChangeStateOfTwoMolecules::transformReactants(Molecule ** reactants, int nReactants)
//{
//	reactants[0]->setState(stateIndex1, newStateValue1);
//	reactants[1]->setState(stateIndex2, newStateValue2);
//	
//	
//}
