#include <iostream>
#include "NFcore.hh"
#include <queue>


using namespace std;
using namespace NFcore;


int Molecule::uniqueIdCount = 0;


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
		bond[b]=NOBOND; indexOfBond[b]=0;
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

	//register this molecule with moleculeType and get some ID values
	ID_complex = this->parentMoleculeType->createComplex(this);
	ID_type = this->parentMoleculeType->getTypeID();
	ID_unique = Molecule::uniqueIdCount++;
	this->listId = listId;
	isAliveInSim = false;
}


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


void Molecule::updateRxnMembership()
{
	parentMoleculeType->updateRxnMembership(this);
}


void Molecule::updateDORRxnValues()
{
	ReactionClass *rxn=0; int rxnIndex=-1, rxnPos=-1;
	for(int i=0; i<parentMoleculeType->getNumOfDORrxns(); i++) {
		rxn=parentMoleculeType->getDORrxn(i);
		rxnIndex=parentMoleculeType->getDORrxnIndex(i);
		rxnPos=parentMoleculeType->getDORrxnPosition(i);

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


void Molecule::removeFromObservables()
{
	parentMoleculeType->removeFromObservables(this);
}


void Molecule::addToObservables()
{
	parentMoleculeType->addToObservables(this);
}


bool Molecule::setPopulation( int count )
{
	if ( isPopulationType()  &&  (count >= 0) )
	{
		population_count = count;
		return true;
	}
	else return false;
}


int Molecule::getPopulation() const
{
	return population_count;
}


bool Molecule::incrementPopulation()
{
	if ( isPopulationType() )
	{
		++population_count;
		return true;
	}
	else return false;
}


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
		if(c!=0)o<<"                  ";
		o<< parentMoleculeType->getComponentName(c) <<"=";
		o<<parentMoleculeType->getComponentStateName(c,component[c]);
		o<<"\tbond=";
		if(bond[c]==NOBOND) o<<"empty";
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

int Molecule::getDegree()
{
	int degree = 0;
	for(int c=0; c<numOfComponents; c++)
		if(bond[c]!=NOBOND) degree++;
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
	if(m1->bond[cIndex1]!=NOBOND || m2->bond[cIndex2]!=NOBOND) {
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

	m1->bond[cIndex] = 0;
	m2->bond[cIndex2] = 0;

	m1->indexOfBond[cIndex] = NOINDEX;
	m2->indexOfBond[cIndex2] = NOINDEX;

	//Handle Complexes
	if(m1->useComplex)
	{
		// NOTE: mergeWithList will handle canonical flags
		m1->getComplex()->updateComplexMembership(m1);
	}
}


queue <Molecule *> Molecule::q;
queue <int> Molecule::d;
list <Molecule *>::iterator Molecule::molIter;


void Molecule::breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth)
{
	if(m==0) {
		cerr<<"Error in Molecule::breadthFirstSearch, m is null.\n";
		cerr<<"Likely an internal error where a MappingSet is on a list and\n";
		cerr<<"is not actually mapped to any molecule!";
		exit(3);
	}

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
		int cMax = cM->numOfComponents;
		for(int c=0; c<cMax; c++)
		{
			//cM->getComp
			if(cM->isBindingSiteBonded(c))
			{
				Molecule *neighbor = cM->getBondedMolecule(c);

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
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;
}


void Molecule::traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit)
{
	Molecule::breadthFirstSearch(members, this, traversalLimit);
}
