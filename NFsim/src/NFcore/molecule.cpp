#include <iostream>
#include "NFcore.hh"
#include <queue>


using namespace std;
using namespace NFcore;

int Molecule::uniqueIdCount = 0;














// Molecule Constructor
//
//
Molecule::Molecule(MoleculeType * parentMoleculeType, int listId)
{
	if(DEBUG) cout<<"-creating molecule instance of type " << parentMoleculeType->getName() << endl;	
	this->parentMoleculeType = parentMoleculeType;
	
	//Set the state of this molecule to default state as defined by this MoleculeType
	this->states = new int [parentMoleculeType->getNumOfStates()];
	for(int s=0; s<parentMoleculeType->getNumOfStates(); s++)
		states[s] = parentMoleculeType->getDefaultState(s);
	
	//Create the default binding sites
	int b = 0;
	this->bonds = new Molecule * [parentMoleculeType->getNumOfBindingSites()];
	for(b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
		bonds[b] = 0;
	
	this->bSiteIndexOfBond = new int [parentMoleculeType->getNumOfBindingSites()];
	for(b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
		bSiteIndexOfBond[b] = -1;	
	
	this->hasVisitedBond = new bool [parentMoleculeType->getNumOfBindingSites()];
	for(b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
		hasVisitedBond[b] = false;
	
	hasVisitedMolecule = false;
	rxnListMappingId = 0;
	nReactions = 0;
	useComplex = parentMoleculeType->getSystem()->isUsingComplex();
	isPrepared = false;
	
	//register this molecule with moleculeType and get some ID values
	//ID_number = this->parentMoleculeType->addMolecule(this);
	ID_complex = this->parentMoleculeType->createComplex(this);
	ID_type = this->parentMoleculeType->getTypeID();
	ID_unique = Molecule::uniqueIdCount++;
	this->listId = listId;
}

// Molecule Deconstructor
//
//
Molecule::~Molecule()
{
	if(DEBUG) cout <<"   -destroying molecule instance of type " << parentMoleculeType->getName() << endl;
	delete [] bonds;
	
	parentMoleculeType = 0;

	StateChangeListener *l;
	while(listeners.size()>0)
	{
		l = listeners.back();
		listeners.pop_back();
		delete l;
	}
	
	delete [] states;
	delete [] bSiteIndexOfBond;
	delete [] hasVisitedBond;
	delete [] rxnListMappingId;
}


void Molecule::prepareForSimulation()
{
	if(isPrepared) return;
	nReactions = parentMoleculeType->getReactionCount();
	this->rxnListMappingId = new int[nReactions];
	for(int r=0; r<nReactions; r++)
		rxnListMappingId[r] = -1;
	isPrepared = true;
}


void Molecule::updateRxnMembership()
{
	parentMoleculeType->updateRxnMembership(this);
}


void Molecule::notifyGroupsThatRateMayChange()
{
	
	//if(listeners.size()>0) cout<<"    notifying groups that rate may have changed"<<endl;
	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
		(*listenerIter)->updateGroupReactions();
	
}

void Molecule::updateDORs()
{
	
	for(int r=0; r<parentMoleculeType->getDORrxnCount(); r++)
	{
		
		ReactionClass * DORrxn = parentMoleculeType->getDORrxn(r);
		int dorRxnIndex = parentMoleculeType->getDORreactantIndex(r);
		int dorRxnPos = parentMoleculeType->getDORreactantPosition(r);
		
	//	cout<<" identified DOR RXN index: "<<dorRxnIndex<<endl;
	//	cout<<" identified DOR RXN pos: "<<dorRxnPos<<endl;
		DORrxn->notifyRateFactorChange(this, dorRxnPos, rxnListMappingId[dorRxnIndex]);
	}
	
}

double Molecule::getDORvalueFromGroup(char * groupName, int valueIndex)
{
	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
	{
		if(strlen(groupName)==strlen((*listenerIter)->getGroupName()))
			if(strncmp(groupName,(*listenerIter)->getGroupName(),strlen(groupName))==0)
				return (*listenerIter)->getValue(valueIndex);
	}
	
	cerr<<"Error!! trying to get DOR value for a group, but no name match!"<<endl;
	cerr<<"    Looking for group: "<<groupName<<" from molecule ";
	cerr<<this->getMoleculeTypeName()<<"_"<<this->getMoleculeID()<<endl;
	exit(1);
}



void Molecule::removeFromObservables()
{
	parentMoleculeType->removeFromObservables(this);
}
void Molecule::addToObservables()
{
	parentMoleculeType->addToObservables(this);
}


void Molecule::setState(const char * state, int value)
{
	int stateIndex = this->parentMoleculeType->getStateIndex(state);
	this->states[stateIndex]=value;
	
	//if(listeners.size()>0) cout<<"Molecule State has changed..."<<endl;
	//Let all the listeners know that the state of a molecule has changed...
	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
		(*listenerIter)->notify(this,stateIndex);
}


void Molecule::setState(int state, int value)
{
	//cout<<"value: "<<value<<"  state: "<<state<<endl;
	this->states[state]=value;
	
	//if(listeners.size()>0) cout<<"Molecule State has changed..."<<endl;
	//Let all the listeners know that the state of a molecule has changed...
	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
		(*listenerIter)->notify(this,state);
}

void Molecule::printDetails() const
{
	int degree = 0;
	cout<<"++ Molecule instance of type: " << parentMoleculeType->getName();
	cout<< " (uId="<<ID_unique<<", mId=" << ID_number << ", tId=" << ID_type << ", cId" << ID_complex<<")"<<endl;
	cout<<"       states: ";
	for(int s=0; s<parentMoleculeType->getNumOfStates(); s++)
	{
		cout<< parentMoleculeType->getStateName(s) <<"=";
		cout<< states[s] << " ";
	}
	cout<<endl<<"       bonds: ";
	for(int b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
	{
		cout<< parentMoleculeType->getBindingSiteName(b) <<"=";
		//cout<<this->hasVisitedBond[b]<<", ";
		if(bonds[b]==0) cout<<"empty ";
		else { degree++; 
			cout<<bonds[b]->getUniqueID()<<" ";
			
			//cout<<"full ";
			} 
	}
	cout<<"  :::degree="<<degree<<endl<<endl;
}

//Get the number of molecules this molecule is bonded to
int Molecule::getDegree()
{
	int degree = 0;
	for(int b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
		if(bonds[b]!=0) degree++;
	return degree;
}

bool Molecule::isBindingSiteOpen(int bIndex) const
{
	if(this->bonds[bIndex]==0) return true;
	return false;
}

bool Molecule::isBindingSiteBonded(int bIndex) const
{
	if(this->bonds[bIndex]==0) return false;
	return true;
}

Molecule * Molecule::getBondedMolecule(int bSiteIndex) const
{
	return bonds[bSiteIndex];	
}


//Molecule * Molecule::getBindingSiteBondParent(int bIndex) const; //{ return bindingSites[bIndex]->getBondedParent(); };
void Molecule::bind(Molecule *m1, int bSiteIndex1, Molecule *m2, int bSiteIndex2)
{
	if(m1->bonds[bSiteIndex1]!=0 || m2->bonds[bSiteIndex2]!=0) { 
		cout<<"Trying to bond "<< m1->getMoleculeTypeName() << "_"<<m1->getUniqueID()<<"(";
		cout<<m1->getMoleculeType()->getBindingSiteName(bSiteIndex1)<<") & ";
		cout<< m2->getMoleculeTypeName()<<"_"<<m2->getUniqueID()<<"(";
		cout<<m2->getMoleculeType()->getBindingSiteName(bSiteIndex2)<<")\n";
		cout<<" to sites that are already occupied!  Check rxn rules!!\n";
		
		m1->printDetails();
		m2->printDetails();
		exit(1); 
	}
	
	m1->bonds[bSiteIndex1] = m2;
	m2->bonds[bSiteIndex2] = m1;
	
	m1->bSiteIndexOfBond[bSiteIndex1] = bSiteIndex2;
	m2->bSiteIndexOfBond[bSiteIndex2] = bSiteIndex1;
	
	//Handle Complexes
	if(m1->useComplex)
	{
		if(m1->getComplex()!=m2->getComplex())
			m1->getComplex()->mergeWithList(m2->getComplex());
	}
}

void Molecule::bind(Molecule *m1, const char * bSiteName1, Molecule *m2, const char * bSiteName2)
{
	int bSiteIndex1 = m1->getMoleculeType()->getBindingSiteIndex(bSiteName1);
	int bSiteIndex2 = m2->getMoleculeType()->getBindingSiteIndex(bSiteName2);
	Molecule::bind(m1, bSiteIndex1, m2, bSiteIndex2);
}


void Molecule::unbind(Molecule *m1, int bSiteIndex)
{
	//get the other molecule bound to this site
	//cout<<"I am here. "<<bSiteIndex<<endl;
	Molecule *m2 = m1->bonds[bSiteIndex];
	if(m2==NULL)
	{
		cout<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cout<<"Trying to unbind a binding site that is not bound!!  Check rxn rules! Quitting."<<endl;
		cout<<endl<<endl<<"The molecule is:"<<endl;
		m1->printDetails();
		cout<<endl<<"The site trying to be unbound was: ";
		cout<<m1->getMoleculeType()->getBindingSiteName(bSiteIndex)<<endl;
		exit(1);
	}
	
	
	//cout<<" and now here.";
	int bSiteIndex2 = m1->bSiteIndexOfBond[bSiteIndex];
	
	//break the bond
	m1->bonds[bSiteIndex] = 0;
	m2->bonds[bSiteIndex2] = 0;
	
	m1->bSiteIndexOfBond[bSiteIndex] = -1;
	m2->bSiteIndexOfBond[bSiteIndex2] = -1;
	
	//Handle Complexes
	if(m1->useComplex)
	{
		m1->getComplex()->updateComplexMembership(m1);
	}
	
	//cout<<" UnBinding!  mol1 complex: ";
	//m1->getComplex()->printDetails();
}

void Molecule::unbind(Molecule *m1, char * bSiteName)
{
	int bSiteIndex = m1->getMoleculeType()->getBindingSiteIndex(bSiteName);
	Molecule::unbind(m1,bSiteIndex);
}




void Molecule::setHasVisited(int bSiteIndex)
{
	hasVisitedBond[bSiteIndex] = true;
}

void Molecule::clear()
{
	hasVisitedMolecule = false;
	for(int i=0; i<parentMoleculeType->getNumOfBindingSites(); i++)
		hasVisitedBond[i] = false;
}




//void Molecule::traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit)
//{	
//	//cout<<"Traversing molecule: " << this->getUniqueID()<<endl;
//	
//		//If we are too deep at this point, then return
//	
//	
//	
//	//cout<<"Traversing, limit: " << traversalLimit << endl;
//	if(hasVisitedMolecule==true) { return; }
//	
//	if(traversalLimit!=ReactionClass::NO_LIMIT)
//	{
//		//cout<<"   Limit: "<<traversalLimit<<endl;
//		if(traversalLimit==0)
//		{
//			//cout<<"     -returning!!! "<<endl;
//			clear();
//			return;
//		}
//	}
//	
//	//Make sure that this is not alreadly in the member list
//	//If it was, then we exit because if it was in the member list,
//	//then we alreadly explored down all its bonds and we have the
//	//entire complex
//	list <Molecule *>::iterator molIter;
//	for( molIter = members.begin(); molIter != members.end(); molIter++ )
//  		if((*molIter)==this) { clear(); return; }
//	
//	//This molecule must then be added to the member list
//	members.push_back(this);
//	
//			
//	//Otherwise, we keep exploring...
//	hasVisitedMolecule = true;
//	
//	for(int b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
//	{
//		if(hasVisitedBond[b]==true) { continue; }
//		if(bonds[b]==0) //binding site is open, so continue
//		{
//			//cout<<"    -bond "<<b<<" is empty."<<endl;
//			continue;
//		}
//		else // binding site has a bond, so we must explore it
//		{
//			//cout<<"    -Going down bond: "<<b<<endl;
//			//get template that is bound to this binding site
//			Molecule * m2 = bonds[b];
//			
//			m2->setHasVisited(bSiteIndexOfBond[b]); //tell t2 we have visited this bond
//			setHasVisited(b); //remember that we visited through this bond.
//			if(traversalLimit==ReactionClass::NO_LIMIT)
//			{
//				//cout<<"continuing with no limit"<<endl;
//				m2->traverseBondedNeighborhood(members, ReactionClass::NO_LIMIT);
//			}
//			else
//			{
//				//cout<<"continuing with one less"<<endl;
//				m2->traverseBondedNeighborhood(members,traversalLimit-1 );
//			}
//		}		
//	}
//	clear();
//}



void Molecule::breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth)
{	
	//Create the queues
	queue <Molecule *> q;
	
	queue <int> d;
	int currentDepth = 0;
	
	//First add this molecule
	q.push(m);
	members.push_back(m);
	d.push(currentDepth+1);
	m->hasVisitedMolecule=true;
	
	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		Molecule *cM = q.front();
		currentDepth = d.front();
		q.pop(); 
		d.pop();
		
		//Make sure the depth does not exceed the limit we want to search
		if((depth!=ReactionClass::NO_LIMIT) && (currentDepth>=depth)) continue;
		
		//Loop through the bonds
		int bMax = cM->getMoleculeType()->getNumOfBindingSites();
		for(int b=0; b<bMax; b++)
		{
			if(cM->isBindingSiteBonded(b))
			{
				Molecule *neighbor = cM->getBondedMolecule(b);
				if(!neighbor->hasVisitedMolecule)
				{
					neighbor->hasVisitedMolecule=true;
					members.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
				}
			}
		}
	}
	
	
	//clear the has visitedMolecule values
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;
	
	
}


void Molecule::traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit)
{	
	
	
	Molecule::breadthFirstSearch(members, this, traversalLimit);
	return;
	//cout<<"Traversing molecule: " << this->getUniqueID()<<endl;
	
		//If we are too deep at this point, then return
	
	
	//cout<<"Traversing, limit: " << traversalLimit << endl;
	if(hasVisitedMolecule==true) { return; }
	
	if(traversalLimit!=ReactionClass::NO_LIMIT)
	{
		//cout<<"   Limit: "<<traversalLimit<<endl;
		if(traversalLimit==0)
		{
			//cout<<"     -returning!!! "<<endl;
			clear();
			return;
		}
	}
	
	//Make sure that this is not alreadly in the member list
	//If it was, then we exit because if it was in the member list,
	//then we alreadly explored down all its bonds and we have the
	//entire complex
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		if((*molIter)==this) { clear(); return; }
	
	//This molecule must then be added to the member list
	members.push_back(this);
	
			
	//Otherwise, we keep exploring...
	hasVisitedMolecule = true;
	
	bool * needToVisit = new bool[parentMoleculeType->getNumOfBindingSites()];
	
	for(int b=0; b<parentMoleculeType->getNumOfBindingSites(); b++)
	{
		if(hasVisitedBond[b]==true) { continue; }
		if(bonds[b]==0) //binding site is open, so continue
		{
			//cout<<"    -bond "<<b<<" is empty."<<endl;
			needToVisit[b]=false;
			continue;
		}
		else // binding site has a bond, so we must explore it
		{
			//cout<<"    -Going down bond: "<<b<<endl;
			//get template that is bound to this binding site
			Molecule * m2 = bonds[b];
			
			m2->setHasVisited(bSiteIndexOfBond[b]); //tell t2 we have visited this bond
			setHasVisited(b); //remember that we visited through this bond.
			if(traversalLimit==ReactionClass::NO_LIMIT)
			{
				//cout<<"continuing with no limit"<<endl;
				needToVisit[b]=true;
				m2->traverseBondedNeighborhood(members, ReactionClass::NO_LIMIT);
			}
			else
			{
				needToVisit[b]=true;
				//cout<<"continuing with one less"<<endl;
				m2->traverseBondedNeighborhood(members,traversalLimit-1 );
			}
		}		
	}
	clear();
}



void Molecule::printMoleculeList(list <Molecule *> &members)
{
	cout<<"List of molecules contains: "<<endl;
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ ) {
		cout<<"   -"<<(*molIter)->getMoleculeTypeName();
		cout<<"_"<<(*molIter)->getMoleculeID();
		cout<<"_u"<<(*molIter)->getUniqueID()<<endl;
	}
}






