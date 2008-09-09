//#include "NFcore.hh"
//#include <math.h>
//
//using namespace std;
//using namespace NFcore;
//
//
//Group::Group(string groupName, System * s, int stateIndex) 
//{
//	this->groupName = groupName;
//	this->stateIndex = stateIndex;
//	this->Group_ID = s->addGroup(this);
//	this->value.push_back(0.0);
//	this->areReactionsUpToDate = true;
//}
//
//Group::~Group() 
//{
//	
//}
//		
//void Group::addToGroup(Molecule * m)
//{
//	//Add molecule to list of molecules
//	this->groupMembers.push_back(m);
//	
//	//Create a StateChangeListener to listen to the molecule
//	new StateChangeListener(m,this,stateIndex);
//	
//	//and update the value of this group based on the value of the state in m
//	value.at(0) += m->getState(stateIndex);
//}
//
//void Group::addToGroupWithoutListener(Molecule * m)
//{
//	//Add molecule to list of molecules
//	this->groupMembers.push_back(m);
//	
//	//We actually still need a listener so that the molecule can get to 
//	//the group, but it will listen to nothing (or state -1)
//	new StateChangeListener(m,this,-1);
//	
//}
//
//
//double Group::getValue(unsigned int valIndex)
//{
//	if((unsigned)valIndex>value.size())
//		return NAN;
//	return value.at(valIndex); 
//}
//
//
//
///* this gets called by a stateChangeListener when a state changes  Remember - it
// * doesn't update reaction rates!  You have to update ReactionRates with a separate
// * call to updateReactionRates! (This is because we need to finish firing the 
// * reaction before we update rates */
//void Group::notify(Molecule *changedMolecule, int oldStateValue, int newStateValue)
//{
//	//cout<<"    Group " << groupName<<"_"<<Group_ID<<" was alerted of change. ";
//	//double oldValue = value;
//	
//	//If nothing changed, catch it now and return
//	if(oldStateValue == newStateValue) return;
//	
//	value.at(0) -= oldStateValue;
//	value.at(0) += newStateValue;
//	
//	areReactionsUpToDate = false;
//	
//	//cout<<" value was "<<oldValue<<" and is "<<value<<endl;
//}
//		
//void Group::updateReactionRates()
//{
//	if(areReactionsUpToDate) { /*cout<<"     in group: rxns up to date, returning.." << endl; */ return; }
//	
//	//Loop through the reactions that
//	
//	//cout<<"       in group: updating rxns"<<endl;
//	for(molIter = groupMembers.begin(); molIter != groupMembers.end(); molIter++ )
//	{
//		(*molIter)->updateDORs();
//	}
//	areReactionsUpToDate = true;
//}
//
//void Group::updateGroupProperty(double * values, int n_values)
//{
//	cout<<"Warning, updating a group property in group class!!  - Indicates that";
//	cout<<" updateGroupProperty function was not implemented in your group class!"<<endl;
//}
//
//
//void Group::printDetails()
//{
//	cout<<"Group: "<<groupName<<"\tgID: "<<Group_ID<<"\tn: "<<groupMembers.size();
//	cout<<"\tValue: "<<value.at(0)<<endl;
//	
//}
//
//
