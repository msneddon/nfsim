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

	//First initialize the component states and bonds
	this->numOfComponents = parentMoleculeType->getNumOfComponents();
	this->component = new int [parentMoleculeType->getNumOfComponents()];
	for(int c=0; c<numOfComponents; c++)
		component[c] = parentMoleculeType->getDefaultComponentState(c);

	this->bond = new Molecule * [numOfComponents];
	this->indexOfBond = new int [numOfComponents];
	for(int b=0; b<numOfComponents; b++) {
		bond[b]=0; indexOfBond[b]=NOBOND;
	}

	hasVisitedMolecule = false;
	hasEvaluatedMolecule = false;
	rxnListMappingId = 0;
	nReactions = 0;
	useComplex = parentMoleculeType->getSystem()->isUsingComplex();
	isPrepared = false;
	isObservable = 0;

	//register this molecule with moleculeType and get some ID values
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
	delete [] bond;

	parentMoleculeType = 0;

//	StateChangeListener *l;
//	while(listeners.size()>0)
//	{
//		l = listeners.back();
//		listeners.pop_back();
//		delete l;
//	}

	delete [] isObservable;
	delete [] component;
	delete [] indexOfBond;
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

	//We do not belong to any observable... yet.
	isObservable=new bool [parentMoleculeType->getNumOfObservables()];
	for(int o=0;o<parentMoleculeType->getNumOfObservables(); o++) {
		isObservable[o]=false;
	}

	localFunctionValues=new double[parentMoleculeType->getNumOfTypeIFunctions()];
	for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
		localFunctionValues[lf]=0;
	}

}



//Used so that this molecule can remember what its local function was
//evaluated to.  Only TypeI local functions are set up in this way
void Molecule::setLocalFunctionValue(double newValue,int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
		cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
		cout<<"index provided was out of bounds!  I shall quit now."<<endl;
		exit(1);
	}
	localFunctionValues[localFunctionIndex] = newValue;
}



void Molecule::updateRxnMembership()
{
	parentMoleculeType->updateRxnMembership(this);
}


void Molecule::notifyGroupsThatRateMayChange()
{
	//if(listeners.size()>0) cout<<"    notifying groups that rate may have changed"<<endl;
//	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
//		(*listenerIter)->updateGroupReactions();

}

//void Molecule::updateDORs()
//{
//
//	for(int r=0; r<parentMoleculeType->getDORrxnCount(); r++)
//	{
//
//		ReactionClass * DORrxn = parentMoleculeType->getDORrxn(r);
//		int dorRxnIndex = parentMoleculeType->getDORreactantIndex(r);
//		int dorRxnPos = parentMoleculeType->getDORreactantPosition(r);
//
//	//	cout<<" identified DOR RXN index: "<<dorRxnIndex<<endl;
//	//	cout<<" identified DOR RXN pos: "<<dorRxnPos<<endl;
//		DORrxn->notifyRateFactorChange(this, dorRxnPos, rxnListMappingId[dorRxnIndex]);
//	}
//
//}

double Molecule::getDORvalueFromGroup(string groupName, int valueIndex)
{
//	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
//	{
//		if(groupName==(*listenerIter)->getGroupName())
//			return (*listenerIter)->getValue(valueIndex);
//	}

	cerr<<"Error!! trying to get DOR value for a group, but no name match!"<<endl;
	cerr<<"    Looking for group: "<<groupName<<" from molecule ";
	cerr<<this->getMoleculeTypeName()<<"_"<<this->getUniqueID()<<endl;
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




void Molecule::setComponentState(int cIndex, int newValue)
{
	this->component[cIndex]=newValue;

	//if(listeners.size()>0) cout<<"Molecule State has changed..."<<endl;
	//Let all the listeners know that the state of a molecule has changed...
	//for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
	//	(*listenerIter)->notify(this,stateIndex);
}
void Molecule::setComponentState(string cName, int newValue) {
	this->component[this->parentMoleculeType->getCompIndexFromName(cName)]=newValue;
}



void Molecule::printDetails() const
{
	int degree = 0;
	cout<<"++ Molecule instance of type: " << parentMoleculeType->getName();
	cout<< " (uId="<<ID_unique << ", tId=" << ID_type << ", cId" << ID_complex<<", degree="<<degree<<")"<<endl;
	cout<<"      components: ";
	for(int c=0; c<numOfComponents; c++)
	{
		if(c!=0)cout<<"                  ";
		cout<< parentMoleculeType->getComponentName(c) <<"=";
		cout<<parentMoleculeType->getComponentStateName(c,component[c]);
		cout<<"\tbond=";
		if(bond[c]==NOBOND) cout<<"empty";
		else cout<<bond[c]->getUniqueID();
		cout<<endl;
	}

	if(parentMoleculeType->getNumOfTypeIFunctions()>0) {
		cout<<"      loc funcs:";
		for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
			if(lf!=0) cout<<"                  ";
			cout<<"  "<<parentMoleculeType->getTypeILocalFunction(lf)->getNiceName();
			cout<<"="<<localFunctionValues[lf]<<"\n";
		}
	}
}

//Get the number of molecules this molecule is bonded to
int Molecule::getDegree()
{
	int degree = 0;
	for(int c=0; c<numOfComponents; c++)
		if(bond[c]!=NOBOND) degree++;
	return degree;
}

bool Molecule::isBindingSiteOpen(int cIndex) const
{
	if(bond[cIndex]==NOBOND) return true;
	return false;
}

bool Molecule::isBindingSiteBonded(int cIndex) const
{
	if(bond[cIndex]==NOBOND) return false;
	return true;
}

Molecule * Molecule::getBondedMolecule(int cIndex) const
{
	return bond[cIndex];
}



void Molecule::bind(Molecule *m1, int cIndex1, Molecule *m2, int cIndex2)
{
	if(m1->bond[cIndex1]!=NOBOND || m2->bond[cIndex2]!=NOBOND) {
		cout<<"Trying to bond "<< m1->getMoleculeTypeName() << "_"<<m1->getUniqueID()<<"(";
		cout<<m1->getMoleculeType()->getComponentName(cIndex1)<<") & ";
		cout<< m2->getMoleculeTypeName()<<"_"<<m2->getUniqueID()<<"(";
		cout<<m2->getMoleculeType()->getComponentName(cIndex2)<<")\n";
		cout<<" to sites that are already occupied!  Check rxn rules!!\n";

		m1->printDetails();
		m2->printDetails();
		exit(1);
	}

	m1->bond[cIndex1] = m2;
	m2->bond[cIndex2] = m1;

	m1->indexOfBond[cIndex1] = cIndex2;
	m2->indexOfBond[cIndex2] = cIndex1;

	//Handle Complexes
	if(m1->useComplex)
	{
		if(m1->getComplex()!=m2->getComplex())
			m1->getComplex()->mergeWithList(m2->getComplex());
	}
}

void Molecule::bind(Molecule *m1, string compName1, Molecule *m2, string compName2)
{
	int cIndex1 = m1->getMoleculeType()->getCompIndexFromName(compName1);
	int cIndex2 = m2->getMoleculeType()->getCompIndexFromName(compName2);
	Molecule::bind(m1, cIndex1, m2, cIndex2);
}


void Molecule::unbind(Molecule *m1, int cIndex)
{
	//get the other molecule bound to this site
	//cout<<"I am here. "<<bSiteIndex<<endl;
	Molecule *m2 = m1->bond[cIndex];
	if(m2==NULL)
	{
		cout<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cout<<"Trying to unbind a binding site that is not bound!!  Check rxn rules! Quitting."<<endl;
		cout<<endl<<endl<<"The molecule is:"<<endl;
		m1->printDetails();
		cout<<endl<<"The site trying to be unbound was: ";
		cout<<m1->getMoleculeType()->getComponentName(cIndex)<<endl;
		exit(1);
	}


	int cIndex2 = m1->indexOfBond[cIndex];

	//break the bond
	m1->bond[cIndex] = NOBOND;
	m2->bond[cIndex2] = NOBOND;

	m1->indexOfBond[cIndex] = NOINDEX;
	m2->indexOfBond[cIndex2] = NOINDEX;

	//Handle Complexes
	if(m1->useComplex)
	{
		m1->getComplex()->updateComplexMembership(m1);
	}

	//cout<<" UnBinding!  mol1 complex: ";
	//m1->getComplex()->printDetails();
}

void Molecule::unbind(Molecule *m1, char * compName)
{
	int cIndex = m1->getMoleculeType()->getCompIndexFromName(compName);
	Molecule::unbind(m1,cIndex);
}








queue <Molecule *> Molecule::q;
queue <int> Molecule::d;
list <Molecule *>::iterator Molecule::molIter;

void Molecule::breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth)
{
	//Create the queues (for effeciency, now queues are a static attribute of Molecule...)
	//queue <Molecule *> q;
	//queue <int> d;
	int currentDepth = 0;

	//cout<<"traversing on:"<<endl;
	//m->printDetails();

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
		int cMax = cM->numOfComponents;
		for(int c=0; c<cMax; c++)
		{
			//cM->getComp
			if(cM->isBindingSiteBonded(c))
			{
				Molecule *neighbor = cM->getBondedMolecule(c);
				//cout<<"looking at neighbor: "<<endl;
				//neighbor->printDetails();
				if(!neighbor->hasVisitedMolecule)
				{
					neighbor->hasVisitedMolecule=true;
					members.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
					//cout<<"adding... to traversal list."<<endl;
				}
			}
		}
	}


	//clear the has visitedMolecule values
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;
}


void Molecule::traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit)
{
	Molecule::breadthFirstSearch(members, this, traversalLimit);
	return;
}



void Molecule::printMoleculeList(list <Molecule *> &members)
{
	cout<<"List of molecules contains: "<<endl;
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ ) {
		cout<<"   -"<<(*molIter)->getMoleculeTypeName();
		cout<<"_u"<<(*molIter)->getUniqueID()<<endl;
	}
}






