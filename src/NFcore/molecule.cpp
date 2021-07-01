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

	// set population type (1 for particle type, 0 for population type)
	this->population_count = ( parentMoleculeType->isPopulationType()  ?  0  :  1 );

	//First initialize the component states and bonds
	this->numOfComponents = parentMoleculeType->getNumOfComponents();
	this->component = new int [parentMoleculeType->getNumOfComponents()];
	for(int c=0; c<numOfComponents; c++)
		component[c] = parentMoleculeType->getDefaultComponentState(c);

	// initialize bond sites
	this->bond = new Molecule * [numOfComponents];
	this->indexOfBond = new int [numOfComponents];
	this->hasVisitedBond = new bool [numOfComponents];
	for(int b=0; b<numOfComponents; b++) {
		bond[b]=0; indexOfBond[b]=NOBOND;
		hasVisitedBond[b] = false;
	}


	hasVisitedMolecule = false;
	hasEvaluatedMolecule = false;
	isMatchedTo=0;
	rxnListMappingId2 = 0;
	nReactions = 0;
	useComplex = parentMoleculeType->getSystem()->isUsingComplex();
	isPrepared = false;
	isObservable = 0;
	localFunctionValues=0;
	//isDead = true;

	//register this molecule with moleculeType and get some ID values
	ID_complex = this->parentMoleculeType->createComplex(this);
	ID_type = this->parentMoleculeType->getTypeID();
	ID_unique = Molecule::uniqueIdCount++;
	this->listId = listId;
	isAliveInSim = false;
}

// Molecule Deconstructor
//
//
Molecule::~Molecule()
{
	if(DEBUG) cout <<"   -destroying molecule instance of type " << parentMoleculeType->getName() << endl;
	delete [] bond;

	parentMoleculeType = 0;


	delete [] isObservable;
	delete [] component;
	delete [] indexOfBond;
	delete [] rxnListMappingId2;
	delete [] hasVisitedBond;

	if(localFunctionValues!=0)
		delete [] localFunctionValues;
}


void Molecule::prepareForSimulation()
{
	if(isPrepared) return;
	nReactions = parentMoleculeType->getReactionCount();
	this->rxnListMappingId2 = new set<int>[nReactions];

	isPrepared = true;

	//We do not belong to any observable... yet.
	isObservable=new int [parentMoleculeType->getNumOfMolObs()];
	for(int o=0;o<parentMoleculeType->getNumOfMolObs(); o++) {
		isObservable[o]=0;
	}


}






void Molecule::setUpLocalFunctionList()
{
	if (parentMoleculeType->getNumOfTypeIFunctions() > 0)
	{
		localFunctionValues=new double[parentMoleculeType->getNumOfTypeIFunctions()];
		for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
			localFunctionValues[lf]=0;
		}
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
	//cout<<"here, mol: "<<this->getMoleculeTypeName()<<"_"<<this->getUniqueID()<<endl;
	//cout<<"localfunctionIndex given: "<< localFunctionIndex<<endl;
	//cout<<"n_localfunctionValue given: "<< localFunctionIndex<<endl;

	//cout<<"localFunctionValues="<<endl;
	//for(int k=0; k<parentMoleculeType->getNumOfTypeIFunctions(); k++) {
	//	cout<<"  "<<localFunctionValues[k]<<endl;
	//}
	//cout<<"done."<<endl;

	localFunctionValues[localFunctionIndex] = newValue;
}


double Molecule::getLocalFunctionValue(int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
		cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
		cout<<"index provided was out of bounds!  I shall quit now."<<endl;
		exit(1);
	}
	return localFunctionValues[localFunctionIndex];
}


LocalFunction * Molecule::getLocalFunction(int localFunctionIndex) {
	if(localFunctionIndex<0 || localFunctionIndex>=parentMoleculeType->getNumOfTypeIFunctions()) {
			cout<<"Error in Molecule: trying to set the value of a local function, but the\n";
			cout<<"index provided was out of bounds!  I shall quit now."<<endl;
			exit(1);
		}
	return parentMoleculeType->getTypeILocalFunction(localFunctionIndex);
}



void Molecule::updateRxnMembership(ReactionClass * r, bool useConnectivity)
{
	if (useConnectivity) {
		parentMoleculeType->updateConnectedRxnMembership(this, r);
	}
	else {
		parentMoleculeType->updateRxnMembership(this);
	}
}

void Molecule::updateTypeIIFunctions()
{
	for (int i=0; i<parentMoleculeType->getNumOfTypeIIFunctions(); i++) {
		parentMoleculeType->getTypeIILocalFunction(i)->evaluateOn(this, LocalFunction::SPECIES);
	}
}

void Molecule::updateTypeIIFunctions( vector <Complex *> & productComplexes )
{
	for (int i=0; i<parentMoleculeType->getNumOfTypeIIFunctions(); i++) {
		vector <Complex *>::iterator complexIter;
		for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter ) {
			parentMoleculeType->getTypeIILocalFunction(i)->evaluateOn( (*complexIter)->getFirstMolecule(), LocalFunction::SPECIES);
		}
	}
}

void Molecule::updateDORRxnValues()
{
	ReactionClass *rxn=0; int rxnIndex=-1, rxnPos=-1;
	//cout<<"Looping over :"<<parentMoleculeType->getNumOfDORrxns()<<" dor rxns "<<endl;
	for(int i=0; i<parentMoleculeType->getNumOfDORrxns(); i++) {
		rxn=parentMoleculeType->getDORrxn(i);
		rxnIndex=parentMoleculeType->getDORrxnIndex(i);
		rxnPos=parentMoleculeType->getDORrxnPosition(i);
		//cout<<"\t\ti="<<i<<" rxn="<<rxn->getName()<<" rxnIndex="<<rxnIndex<<" rxnPos="<<rxnPos<<endl;

		if(isPrepared) {
			//If we are in this reaction, then we have to update our value...
			if(getRxnListMappingId(rxnIndex)>=0) {
				//iterate over all mappings
				set<int> tempSet = getRxnListMappingSet(rxnIndex);
				//iterate over all agent-mappings  for the same reaction
				for(set<int>::iterator it= tempSet.begin();it!= tempSet.end(); ++it){

				//Careful here!  remember to update the propensity of this
				//reaction in the system after we notify of the rate factor change!
					double oldA = rxn->get_a();
					rxn->notifyRateFactorChange(this,rxnPos,*it);
					parentMoleculeType->getSystem()->update_A_tot(rxn,oldA,rxn->update_a());
				}
			}
		}
	}
}

///////////////
//  MOLECULE_DEPENDENT_UPDATE_ADDITION
//void Molecule::addDependentUpdateMolecule(Molecule *m) {
//	for(molIter=dependentUpdateMolecules.begin();molIter!=dependentUpdateMolecules.end();molIter++)
//		if((*molIter)->getUniqueID()==m->getUniqueID())
//			return;
//	dependentUpdateMolecules.push_back(m);
//}
//void Molecule::removeDependentUpdateMolecule(Molecule *m) {
//	for(molIter=dependentUpdateMolecules.begin();molIter!=dependentUpdateMolecules.end();molIter++)
//		if((*molIter)->getUniqueID()==m->getUniqueID()) {
//			dependentUpdateMolecules.erase(molIter);
//		}
//}
////////////////



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

//double Molecule::getDORvalueFromGroup(string groupName, int valueIndex)
//{
////	for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
////	{
////		if(groupName==(*listenerIter)->getGroupName())
////			return (*listenerIter)->getValue(valueIndex);
////	}
//
//	cerr<<"Error!! trying to get DOR value for a group, but no name match!"<<endl;
//	cerr<<"    Looking for group: "<<groupName<<" from molecule ";
//	cerr<<this->getMoleculeTypeName()<<"_"<<this->getUniqueID()<<endl;
//	exit(1);
//}



void Molecule::removeFromObservables()
{
	parentMoleculeType->removeFromObservables(this);
}
void Molecule::addToObservables()
{
	parentMoleculeType->addToObservables(this);
}


// set population
bool Molecule::setPopulation( int count )
{
	if ( isPopulationType()  &&  (count >= 0) )
	{
		population_count = count;
		return true;
	}
	else return false;
}

// get popualtion
int Molecule::getPopulation() const
{
	return population_count;
}

// increase population by one
bool Molecule::incrementPopulation()
{
	if ( isPopulationType() )
	{
		++population_count;
		return true;
	}
	else return false;
}

// decrease population by one
bool Molecule::decrementPopulation()
{
	if ( isPopulationType()  &&  (population_count > 0) )
	{
		--population_count;
		return true;
	}
	else return false;

}


void Molecule::setComponentState(int cIndex, int newValue)
{
	this->component[cIndex]=newValue;
	if (useComplex)
		// Need to manually unset canonical flag since we're not calling a Complex method
		getComplex()->unsetCanonical();

	//if(listeners.size()>0) cout<<"Molecule State has changed..."<<endl;
	//Let all the listeners know that the state of a molecule has changed...
	//for(listenerIter = listeners.begin(); listenerIter != listeners.end(); listenerIter++ )
	//	(*listenerIter)->notify(this,stateIndex);
}
void Molecule::setComponentState(string cName, int newValue) {
	this->component[this->parentMoleculeType->getCompIndexFromName(cName)]=newValue;

	if (useComplex)
		// Need to manually unset canonical flag since we're not calling a Complex method
		getComplex()->unsetCanonical();

}


void Molecule::printDetails() {
	this->printDetails(cout);
}
void Molecule::printDetails(ostream &o)
{
	int degree = 0;
	o<<"++ Molecule instance of type: " << parentMoleculeType->getName();
	o<< " (uId="<<ID_unique << ", tId=" << ID_type << ", cId=" << ID_complex<<", degree="<<degree<<")"<<endl;
	o<<"      components: ";
	for(int c=0; c<numOfComponents; c++)
	{
		// Do not print non-bonded states so that mRNA representations are compact
		// Arvind Rasi Subramaniam
		if (bond[c] == 0) continue;
		if(c!=0)o<<"                  ";
		o<< parentMoleculeType->getComponentName(c) <<"=";
		o<<parentMoleculeType->getComponentStateName(c,component[c]);
		o<<"\tbond=";
		if(bond[c]==nullptr) o<<"empty";
		else {
			o<<bond[c]->getMoleculeTypeName()<<"_"<<bond[c]->getUniqueID();
			o<<"("<<bond[c]->getMoleculeType()->getComponentName(this->indexOfBond[c])<<")";
		}

		o<<endl;
	}

	o.flush();
	if(parentMoleculeType->getNumOfTypeIFunctions()>0) {
		o<<"      loc funcs:";
		for(int lf=0; lf<parentMoleculeType->getNumOfTypeIFunctions(); lf++) {
			if(lf!=0) o<<"                  ";
			o<<"  "<<parentMoleculeType->getTypeILocalFunction(lf)->getNiceName();
			o<<"="<<localFunctionValues[lf]<<"\n";
		}
	}
}

/**
 * Prints to screen, see printDetails(ostram )
 * @author Arvind Rasi Subramaniam
 */
void Molecule::printBondDetails() {
	this->printBondDetails(cout);
}

/**
 * Print all bonded states and their details to output stream
 *
 * Iterates through all components of a molecule type and prints
 * only those components with a bond including the binding partner
 * and the site of attachment.
 *
 * Same function written separately for printing to cout and NFStream
 * @param o - Stream to write to
 * @author Arvind Rasi Subramaniam
 */
void Molecule::printBondDetails(ostream &o) {
	int degree = 0;
	o<< parentMoleculeType->getName() << "\t"<<ID_unique;
//	o<<"\t";
//	for(int c=0; c<numOfComponents; c++)
//	{
//		if(bond[c]==NOBOND) {continue;}
//		else {
//			o<<"\t";
//			o<< parentMoleculeType->getComponentName(c) <<"=";
//			o<<parentMoleculeType->getComponentStateName(c,component[c]);
//			o<<",bond=";
//			o<<bond[c]->getMoleculeTypeName()<<"_"<<bond[c]->getUniqueID();
//			o<<"("<<bond[c]->getMoleculeType()->getComponentName(this->indexOfBond[c])<<")";
//		}
//	}
	o.flush();
}

/**
 * Same as printBondDetails(ostream )
 * but writes to NFStream file instead of ostream
 * @author Arvind Rasi Subramaniam
 */
void Molecule::printBondDetails(NFstream &o)
{
	if (parentMoleculeType->getSystem()->getRxnNumberTrack()) {
		o << parentMoleculeType->getTypeID() << "\t" << ID_unique;
	} else {
		o << parentMoleculeType->getName() << "\t" << ID_unique;
	}
	if (parentMoleculeType->getSystem()->getTrackConnected()) {
		o<<"\t";
		for(int c=0; c<numOfComponents; c++)
		{
			if(bond[c] == NULL) {continue;}
			else {
				o<<"||";
				o << parentMoleculeType->getComponentName(c);
				if (parentMoleculeType->getComponentStateName(c,component[c]) != "NO_STATE") {
					o<< "-" << parentMoleculeType->getComponentStateName(c,component[c]);
				}
				o<<":";
				o<<bond[c]->getMoleculeType()->getComponentName(this->indexOfBond[c]);
				o<<"-"<<bond[c]->getMoleculeTypeName()<<"_"<<bond[c]->getUniqueID();
			}
		}
	}
	o.flush();
}

//Get the number of molecules this molecule is bonded to
int Molecule::getDegree()
{
	int degree = 0;
	for(int c=0; c<numOfComponents; c++)
		if(bond[c]!=nullptr) degree++;
	return degree;
}

// Get a label for this molecule or one of it's components (labels are not unique)
//   cIndex==-1  =>  get label for molecule, "m:typename"
//   cIndex>=0   =>  get label for component cIndex, "c:name~state"
string Molecule::getLabel ( int cIndex ) const
{
    string label("");
    if ( cIndex < 0 )
    {	// molecule label
        label += "m:" + getMoleculeTypeName();
    }
    else
    {	// component label
        label += "c:" + (   getMoleculeType()->isEquivalentComponent(cIndex)
        		          ? getMoleculeType()->getEquivalenceClassComponentNameFromComponentIndex(cIndex)
        		          : getMoleculeType()->getComponentName(cIndex) )
        		      + "~" + getMoleculeType()->getComponentStateName(cIndex, getComponentState(cIndex));
    }
    return label;
}


bool Molecule::isBindingSiteOpen(int cIndex) const
{
	if(bond[cIndex]==nullptr) return true;
	return false;
}

bool Molecule::isBindingSiteBonded(int cIndex) const
{
	if(bond[cIndex]==nullptr) return false;
	return true;
}

Molecule * Molecule::getBondedMolecule(int cIndex) const
{
	return bond[cIndex];
}

// given the component index, look up what we are bonded to.  Then
// in the molecule we are bonded to, look at what site we are bonded to
//
//  for instance
//
//    this(a!1).other(b!),   and we call this->getBondedMoleculeBindingSiteIndex(0)
//
//    where index 0 = this site a, then this function would return the
//    component index of b in molecule other.
//
int Molecule::getBondedMoleculeBindingSiteIndex(int cIndex) const
{
	return indexOfBond[cIndex];
}



void Molecule::bind(Molecule *m1, int cIndex1, Molecule *m2, int cIndex2)
{
	if(m1->bond[cIndex1]!=nullptr || m2->bond[cIndex2]!=nullptr) {
		cerr<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"Your universal traversal limit was probably set too low, so some updates were not correct!\n\n";

		cerr<<"Trying to bond "<< m1->getMoleculeTypeName() << "_"<<m1->getUniqueID()<<"(";
		cerr<<m1->getMoleculeType()->getComponentName(cIndex1)<<") & ";
		cerr<< m2->getMoleculeTypeName()<<"_"<<m2->getUniqueID()<<"(";
		cerr<<m2->getMoleculeType()->getComponentName(cIndex2)<<")\n";
		cerr<<" to sites that are already occupied!  Check rxn rules!!\n";
		cerr<<"\n";
		m1->printDetails(cerr);
		m2->printDetails(cerr);
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
		{
			// NOTE: mergeWithList will handle canonical flags
			m1->getComplex()->mergeWithList(m2->getComplex());
		}
		else
			// Need to manually unset canonical flag since we're not calling a Complex method
			m1->getComplex()->unsetCanonical();
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
		cerr<<endl<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"Your universal traversal limit was probably set too low, so some updates were not correct!"<<endl;
		cerr<<"Trying to unbind a binding site that is not bound!!  Check rxn rules, and traversal limits! Quitting."<<endl;
		cerr<<endl<<endl<<"The molecule is:"<<endl;
		m1->printDetails(cerr);
		cerr<<endl<<"The site trying to be unbound was: ";
		cerr<<m1->getMoleculeType()->getComponentName(cIndex)<<endl;
		exit(3);
	}


	int cIndex2 = m1->indexOfBond[cIndex];

	//break the bond (older compilers don't let you assign NOBOND to type molecule)
	m1->bond[cIndex] = 0; //NOBOND;
	m2->bond[cIndex2] = 0; //NOBOND;

	m1->indexOfBond[cIndex] = NOINDEX;
	m2->indexOfBond[cIndex2] = NOINDEX;

	//Handle Complexes
	if(m1->useComplex)
	{
		// NOTE: mergeWithList will handle canonical flags
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
// This might need to change to list if I'm reversing vector changes
void Molecule::breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth)
{
	if(m==0) {
		cerr<<"Error in Molecule::breadthFirstSearch, m is null.\n";
		cerr<<"Likely an internal error where a MappingSet is on a list and\n";
		cerr<<"is not actually mapped to any molecule!";
		exit(3);
	}

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
	//always call breadth first search, it is a bit faster
	//if(traversalLimit>=0)
		Molecule::breadthFirstSearch(members, this, traversalLimit);
	//else
	//	this->depthFirstSearch(members);
}

//Isn't ever called really, but is availabe.  Note that it cannot use traversal limits
//because it is depth first
void Molecule::depthFirstSearch(list <Molecule *> &members)
// void Molecule::depthFirstSearch(vector <Molecule *> &members)
{
	if(this->hasVisitedMolecule==true) {
		return;
	}

	this->hasVisitedMolecule=true;
	members.push_back(this);

	int cMax = this->numOfComponents;
	for(int c=0; c<cMax; c++)
	{
		if(hasVisitedBond[c]==true) continue;
		if(this->isBindingSiteBonded(c))
		{
			Molecule *neighbor = this->getBondedMolecule(c);
			neighbor->hasVisitedBond[indexOfBond[c]]=true;
			hasVisitedBond[c]=true;
			neighbor->depthFirstSearch(members);
		}
	}

	//clear things out
	hasVisitedMolecule = false;
	for(int c=0; c<numOfComponents; c++)
		hasVisitedBond[c] = false;
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






