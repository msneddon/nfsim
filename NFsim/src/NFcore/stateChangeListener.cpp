//#include "NFcore.hh"
//
//
//using namespace std;
//using namespace NFcore;
//
//
//
//StateChangeListener::StateChangeListener(Molecule * m, Group * g, int stateIndex)
//{
//	this->molecule = m;
//	this->group = g;
//	this->stateIndex = stateIndex;
//	if(stateIndex<0)
//	{
//		this->oldStateValue = 0;
//		this->newStateValue = 0;
//	}
//	else
//	{
//		
//		this->oldStateValue = m->getComponentState(stateIndex);
//		this->newStateValue = this->oldStateValue;
//	}
//	
//	
//	
//	
//	//Tell molecule he has a listener!
//	m->addListener(this);
//}
//
//StateChangeListener::~StateChangeListener()
//{
//	molecule = 0;
//	group = 0;
//}
//	
///* called by molecule to alert that state has changed */
//void StateChangeListener::notify(Molecule *changedMolecule, int stateIndexOfChange) 
//{
//	//cout<<"  -StateChangeListener heard that! Message came from: "<< molecule->getMoleculeID()<<endl;
//	if(stateIndexOfChange!=stateIndex) {
//		//cout<<"    but state I am listening for did not change, returning..."<<endl;
//		return;
//	}
//	
//	// set old state value
//	oldStateValue = newStateValue;
//	newStateValue = molecule->getComponentState(stateIndex);
//	
//	//Tell the group that this value has changed
//	//And the group will take care of the rest
//	group->notify(changedMolecule, oldStateValue, newStateValue);
//}
//
//void StateChangeListener::updateGroupReactions()
//{
//	//cout<<"In listener: calling update reaction rate from group"<<endl;
//	group->updateReactionRates();
//}
//
