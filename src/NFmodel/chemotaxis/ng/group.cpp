#include "ng.hh"





using namespace NG;



DimerGroup::DimerGroup(char *name, System *s, int methStateIndex, unsigned int numOfMethSites) : Group(name, s, methStateIndex)
{
	isFull = false;
	this->numOfMethSites = numOfMethSites;
	//if(numOfMethSites>8)
	//{
	//	cerr<<"Error: when creating group DimerGroup, you are trying to use more than 8 methylation sites!"<<endl;
	//	cerr<<"I don't think I can properly handle that yet, so I am quitting."<<endl; exit(1);
	//}
	
	//Value vector contains our P_ON and P_OFF values, so we have to
	//add those two values to the value vector.  I clear it first just for fun.
	value.clear();
	value.push_back(0);  //pushing P_on
	value.push_back(0);  //pushing P_off
}



DimerGroup::~DimerGroup() 
{
}

//When we call this, we are always adding receptor dimers!!!
//to add CheA to this group, call the function addToGroupWithoutListener
//so that we don't listen for CheA's methylation state
void DimerGroup::addToGroup(Molecule * m)
{
	//cout<<"Trying to add: "<<endl; m->printDetails();
	if(isFull)
	{
		cerr<<"Adding more than one receptor to this DimerGroup!"<<endl;
		exit(1);	
	}
	
	//Get information about this receptor
	int ReceptorType = m->getState(m->getMoleculeType()->getStateIndex("type"));
	int methLevel = m->getState(m->getMoleculeType()->getStateIndex("m"));
	
	if(ReceptorType==TAR)
	{
	}
	else if(ReceptorType==TSR)
	{
	}
	else
	{
		cerr<<"Trying to add: "<<endl; m->printDetails();
		cerr<<"Which is not a receptor of type TAR or TSR.  I just can't do that yet."<<endl;
		exit(1);
	}
	
	//Add this molecule to the group if we have successfully updated everything
	this->groupMembers.push_back(m);
	
	//Create a StateChangeListener to listen to the molecule
	new StateChangeListener(m,this,this->stateIndex);
	
	
	//and update the value of this group
	value.at(METH_SITE_LEVEL) = methLevel * m->getDORvalueFromGroup(CLUSTER_NAME, P_ON_INDEX);
	value.at(FREE_SITE_LEVEL) = (numOfMethSites-methLevel) *m->getDORvalueFromGroup(CLUSTER_NAME, P_OFF_INDEX);
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false;

	isFull = true;
}





/* this gets called by a stateChangeListener when a state changes.  Remember - it
 * doesn't update reaction rates!  You have to update ReactionRates with a separate
 * call to updateReactionRates! (This is because we need to finish firing the 
 * reaction before we update rates)  */
void DimerGroup::notify(Molecule *changedMolecule, int oldStateValue, int newStateValue)
{
	//A simple check
	if(newStateValue>numOfMethSites || newStateValue<0) 
	{
		cerr<<"Num of meth sites: " << numOfMethSites << endl;
		cerr<<"Error in NG_dimerGroup: going from meth level "<<oldStateValue<<" to level ->"<<newStateValue<<endl; 
		exit(1);
	}


	//and update the value of this group
	
	//changedMolecule->getDORvalueFromGroup(CLUSTER_NAME, P_OFF_INDEX);
	
	//cout<<"hereFirst:"<<endl;
	value.at(METH_SITE_LEVEL) = newStateValue*changedMolecule->getDORvalueFromGroup(CLUSTER_NAME, P_ON_INDEX);
	value.at(FREE_SITE_LEVEL) = (numOfMethSites-newStateValue)* changedMolecule->getDORvalueFromGroup(CLUSTER_NAME, P_OFF_INDEX);
	//cout<<"and here:"<<endl;
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false; 
}

//This is the function that gets called if a ligand concentration
//is changed.  By convention, values will be the same length as the number
//of ligand molecules we are trying to sense.  For now, this means that
//it is only of length one, and that one value is Aspartate concentration
void DimerGroup::updateGroupProperty(double * values, int n_values)
{
	//we do need to update anything! so do nothing!
}



void DimerGroup::printDetails()
{
	cout<<"NG_DimerGroup " << this->Group_ID <<" : " <<endl;
	cout<<"   -Meth Site Count: " << value.at(METH_SITE_LEVEL)<<endl;
	cout<<"   -Free Site Count: " << value.at(FREE_SITE_LEVEL)<<endl;
	cout<<endl;
	
	
}
