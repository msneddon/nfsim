#ifndef NFREACTIONS_HH_
#define NFREACTIONS_HH_

//#include "../../NFcore/NFcore.hh"
//
//using namespace NFcore;
//
////Take two molecules and bind them, does not yet have a
////constructor that takes name of binding sites instead of an index
//class ReactionSimpleBinding : public ReactionClass {
//
//	public:
//		ReactionSimpleBinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int bSiteIndex1,
//								int bSiteIndex2);
//		ReactionSimpleBinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char* bSiteName1,
//								char* bSiteName2);
//		virtual ~ReactionSimpleBinding() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex1;
//		int bSiteIndex2;
//};
//
//
////Take one molecule and unbind a specific site
//class ReactionUnbinding : public ReactionClass {
//
//	public:
//		ReactionUnbinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int bSiteIndex);
//		ReactionUnbinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bSiteName);
//		virtual ~ReactionUnbinding() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//};
//
//
//////Take one molecule and unbind a specific site
//class SymmetricReactionUnbinding : public ReactionClass {
//
//	public:
//		SymmetricReactionUnbinding(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int bSiteIndex);
//		virtual ~SymmetricReactionUnbinding() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//};
//
////Take one molecule and change its state, has constructors that
////accept either the state index or state name
//class ReactionChangeState : public ReactionClass{
//
//	public:
//		ReactionChangeState(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								int stateIndex,
//								int newStateValue);
//		ReactionChangeState(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName,
//								int newStateValue);
//		~ReactionChangeState() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int stateIndex;
//		int newStateValue;	
//};
//
//class ReactionChangeTwoStates : public ReactionClass{
//
//	public:
//		ReactionChangeTwoStates(	char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2);
//		~ReactionChangeTwoStates() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int stateIndex1;
//		int newStateValue1;
//		int stateIndex2;
//		int newStateValue2;	
//};
//
//
//
////Take two molecules, bind them, and change state of one
//// !!! not implemented yet !!!
//class ReactionBindAndChangeState {
//
//	public:
//		//Input binding site 1, binding site 2, 
//		ReactionBindAndChangeState() {};
//		~ReactionBindAndChangeState() {};
//	private:
//		int bSiteIndex1;
//		int bSiteIndex2;
//		int stateChangeMoleculeReactantIndex;
//		int stateIndex;
//		int stateValue;
//};
//
//
////Take one molecule, unbind a binding site, and change its state
//class ReactionUnbindAndChangeState : public ReactionClass {
//
//	public:
//		//Input:  binding site to unbind, state to change, new state
//		ReactionUnbindAndChangeState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName,
//								int newStateValue);
//		~ReactionUnbindAndChangeState() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//		int stateIndex;
//		int newStateValue;
//};
//
//
////Take one molecule, unbind a binding site, and change its state and a state of the
////molecule that is now unbound
//class ReactionUnbindAndChangeBothState : public ReactionClass {
//
//	public:
//		//Input:  binding site to unbind, state to change, new state
//		ReactionUnbindAndChangeBothState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2);
//		~ReactionUnbindAndChangeBothState() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//		int stateIndex1;
//		int newStateValue1;
//		int stateIndex2;
//		int newStateValue2;
//};
//
////Take one molecule, unbind a binding site, and decrease the value of one of its states
//class ReactionUnbindAndDecrementState : public ReactionClass {
//
//	public:
//		//Input:  binding site to unbind, state to change, new state
//		ReactionUnbindAndDecrementState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName);
//		~ReactionUnbindAndDecrementState() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//		int stateIndex;
//};
//
////Take one molecule, unbind a binding site, and increment the value of one of its states
//class ReactionUnbindAndIncrementState : public ReactionClass {
//
//	public:
//		//Input:  binding site to unbind, state to change, new state
//		ReactionUnbindAndIncrementState(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * bindingSiteName,
//								char * stateName);
//		~ReactionUnbindAndIncrementState() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int bSiteIndex;
//		int stateIndex;
//};
//
//
//
//
//
//
//
//
////Take two molecules, change a state in each of the two molecules
//// (this is good for transfer reactions where a group is transfered from one molecule to the other)
//class ReactionChangeStateOfTwoMolecules : public ReactionClass {
//
//	public:
//		//Input:  binding site to unbind, state to change, new state
//		ReactionChangeStateOfTwoMolecules(char * name, 
//								int n_reactants, 
//								TemplateMolecule ** reactantTemplates, 
//								double rate,
//								char * stateName1,
//								int newStateValue1,
//								char * stateName2,
//								int newStateValue2);
//		~ReactionChangeStateOfTwoMolecules() {};
//		virtual void transformReactants(Molecule ** reactants, int nReactants);
//	private:
//		int stateIndex1;
//		int newStateValue1;
//		int stateIndex2;
//		int newStateValue2;
//};
//
//
//






#endif /*NFREACTIONS_HH_*/
