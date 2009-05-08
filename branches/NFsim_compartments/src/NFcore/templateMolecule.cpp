#include <iostream>

#include "templateMolecule.hh"

using namespace std;
using namespace NFcore;



queue <TemplateMolecule *> TemplateMolecule::q;
queue <int> TemplateMolecule::d;
vector <TemplateMolecule *>::iterator TemplateMolecule::tmVecIter;
list <TemplateMolecule *>::iterator TemplateMolecule::tmIter;

int TemplateMolecule::TotalTemplateMoleculeCount=0;


/*! Only constructor for TemplateMolecules */
TemplateMolecule::TemplateMolecule(MoleculeType * moleculeType){
	this->moleculeType=moleculeType;
	this->uniqueTemplateID=TotalTemplateMoleculeCount++;

	this->n_mapGenerators=0;
	this->mapGenerators=new MapGenerator*[0];

	//Start everything off with no constraints, then we will reinitialize
	//every time we add some type of constraint
	this->n_emptyComps=0;
	this->emptyComps=new int[0];
	this->n_occupiedComps=0;
	this->occupiedComps=new int[0];
	this->n_compStateConstraint=0;
	this->compStateConstraint_Comp=new int[0];
	this->compStateConstraint_Constraint=new int[0];
	this->n_compStateExclusion=0;
	this->compStateExclusion_Comp=new int[0];
	this->compStateExclusion_Exclusion=new int[0];

	this->n_bonds=0;
	this->bondComp=new int[0];
	this->bondCompName=new string[0];
	this->bondPartner=new TemplateMolecule * [0];
	this->bondPartnerCompName=new string[0];
	this->bondPartnerCompIndex=new int[0];
	this->hasVisitedBond=new bool[0];

	this->n_connectedTo=0;
	this->connectedTo=new TemplateMolecule*[n_connectedTo];
	this->hasTraversedDownConnectedTo=new bool[n_connectedTo];
	this->otherTemplateConnectedToIndex=new int[n_connectedTo];
	this->connectedToHasRxnCenter=new bool[n_connectedTo];



	//Init symmetric site matchers...
	this->n_symComps=0;
	this->symCompName=new string[0];
	this->symCompUniqueId=new string[0];
	this->symCompStateConstraint=new int[0];
	this->symCompBoundState=new int[0];  //either Empty (0) or Occupied (1);
	this->symBondPartner=new TemplateMolecule * [0];
	this->symBondPartnerCompName=new string[0];
	this->symBondPartnerCompIndex=new int[0];
	this->hasTraversedDownSym=new bool[0];


	//Get information from the moleculeType as to how we might
	//match up each of its components
	this->n_totalComps=moleculeType->getNumOfComponents();
	this->isSymCompMapped=new bool[n_totalComps];
	this->compIsAlwaysMapped=new bool[n_totalComps];
	for(int i=0; i<n_totalComps; i++) {
		isSymCompMapped[i]=false;
		compIsAlwaysMapped[i]=false;
	}

	//
	this->matchMolecule=0;
	this->hasVisitedThis=false;

	//finally, we have to register this template molecule with the molecule
	//type so that we can easily destroy them at the end.
	this->moleculeType->addTemplateMolecule(this);
}


TemplateMolecule::~TemplateMolecule() {

	for(int i=0; i<n_mapGenerators; i++) {
		delete mapGenerators[i];
		mapGenerators[i]=0;
	}
	delete [] mapGenerators;

	delete [] emptyComps;
	delete [] occupiedComps;
	delete [] compStateConstraint_Comp;
	delete [] compStateConstraint_Constraint;
	delete [] compStateExclusion_Comp;
	delete [] compStateExclusion_Exclusion;
	delete [] symCompUniqueId;

	delete [] bondComp;
	delete [] bondCompName;
	delete [] bondPartner;
	delete [] bondPartnerCompName;
	delete [] bondPartnerCompIndex;
	delete [] hasVisitedBond;

	delete [] connectedTo;
	delete [] hasTraversedDownConnectedTo;
	delete [] otherTemplateConnectedToIndex;
	delete [] connectedToHasRxnCenter;

	delete [] symCompName;
	delete [] symCompBoundState;
	delete [] symCompStateConstraint;
	delete [] symBondPartner;
	delete [] symBondPartnerCompName;
	delete [] symBondPartnerCompIndex;
	delete [] hasTraversedDownSym;

	delete [] isSymCompMapped;
	delete [] compIsAlwaysMapped;

	matchMolecule=0;
	hasVisitedThis=false;
}


string TemplateMolecule::getMoleculeTypeName() const {
	return moleculeType->getName();
}


void TemplateMolecule::addMapGenerator(MapGenerator *mg) {
	MapGenerator ** newMapGenerators = new MapGenerator*[n_mapGenerators+1];
	for(int i=0; i<n_mapGenerators; i++) {
		newMapGenerators[i]=mapGenerators[i];
	}
	newMapGenerators[n_mapGenerators]=mg;
	delete [] mapGenerators;
	mapGenerators = newMapGenerators;
	n_mapGenerators++;
}












void TemplateMolecule::addEmptyComponent(string cName) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add empty binding site of a symmetric component with addEmptyComponent() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);
	int *newEmptyCompArray=new int[n_emptyComps+1];
	for(int i=0; i<n_emptyComps; i++)
		newEmptyCompArray[i]=this->emptyComps[i];
	newEmptyCompArray[n_emptyComps]=compIndex;
	delete [] emptyComps;
	emptyComps=newEmptyCompArray;
	n_emptyComps++;
	compIsAlwaysMapped[compIndex]=true;
}
void TemplateMolecule::addBoundComponent(string cName) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add bound binding site of a symmetric component with addBoundComponent() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);
	int *newOccupiedCompArray=new int[n_occupiedComps+1];
	for(int i=0; i<n_occupiedComps; i++)
		newOccupiedCompArray[i]=this->occupiedComps[i];
	newOccupiedCompArray[n_occupiedComps]=compIndex;
	delete [] occupiedComps;
	occupiedComps=newOccupiedCompArray;
	n_occupiedComps++;
	compIsAlwaysMapped[compIndex]=true;
}
void TemplateMolecule::addComponentConstraint(string cName, string stateName) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add component constraint of a symmetric component with addComponentConstraint() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);
	int stateValue=moleculeType->getStateValueFromName(compIndex,stateName);
	addComponentConstraint(cName,stateValue);

}
void TemplateMolecule::addComponentConstraint(string cName, int stateValue) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add component constraint of a symmetric component with addComponentConstraint() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);

	int *newConstraint_Comp=new int[n_compStateConstraint+1];
	int *newConstraint_Constraint=new int[n_compStateConstraint+1];
	for(int i=0; i<n_compStateConstraint; i++) {
		newConstraint_Comp[i]=compStateConstraint_Comp[i];
		newConstraint_Constraint[i]=compStateConstraint_Constraint[i];
	}
	newConstraint_Comp[n_compStateConstraint]=compIndex;
	newConstraint_Constraint[n_compStateConstraint]=stateValue;
	delete [] compStateConstraint_Comp;
	delete [] compStateConstraint_Constraint;
	compStateConstraint_Comp=newConstraint_Comp;
	compStateConstraint_Constraint=newConstraint_Constraint;
	n_compStateConstraint++;
	compIsAlwaysMapped[compIndex]=true;
}
void TemplateMolecule::addComponentExclusion(string cName, string stateName) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add component exclusion of a symmetric component with addComponentExclusion() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);
	int stateValue=moleculeType->getStateValueFromName(compIndex,stateName);
	addComponentExclusion(cName,stateValue);
}
void TemplateMolecule::addComponentExclusion(string cName, int stateValue) {
	if(moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add component exclusion of a symmetric component with addComponentExclusion() function.");
	}
	int compIndex=moleculeType->getCompIndexFromName(cName);

	int *newExclusion_Comp=new int[n_compStateExclusion+1];
	int *newExclusion_Exclusion=new int[n_compStateExclusion+1];
	for(int i=0; i<n_compStateExclusion; i++) {
		newExclusion_Comp[i]=compStateExclusion_Comp[i];
		newExclusion_Exclusion[i]=compStateExclusion_Exclusion[i];
	}
	newExclusion_Comp[n_compStateExclusion]=compIndex;
	newExclusion_Exclusion[n_compStateExclusion]=stateValue;
	delete [] compStateExclusion_Comp;
	delete [] compStateExclusion_Exclusion;
	compStateExclusion_Comp=newExclusion_Comp;
	compStateExclusion_Exclusion=newExclusion_Exclusion;
	n_compStateExclusion++;
	compIsAlwaysMapped[compIndex]=true;
}

void TemplateMolecule::clearConnectedTo()
{
	delete [] connectedTo;
	delete [] hasTraversedDownConnectedTo;
	delete [] otherTemplateConnectedToIndex;
	delete [] connectedToHasRxnCenter;
	this->n_connectedTo=0;
	this->connectedTo=new TemplateMolecule*[n_connectedTo];
	this->hasTraversedDownConnectedTo=new bool[n_connectedTo];
	this->otherTemplateConnectedToIndex=new int[n_connectedTo];
	this->connectedToHasRxnCenter=new bool[n_connectedTo];
}

void TemplateMolecule::addConnectedTo(TemplateMolecule *t2, int otherConToIndex) {
	addConnectedTo(t2,otherConToIndex,false);
}
void TemplateMolecule::addConnectedTo(TemplateMolecule *t2, int otherConToIndex,bool otherHasRxnCenter) {

	TemplateMolecule **newConnectedTo = new TemplateMolecule * [n_connectedTo+1];
	bool * newHasTraversedDownConnectedTo = new bool[n_connectedTo+1];
	int * newOtherTemplateConnectedToIndex = new int[n_connectedTo+1];
	bool *newConnectedToHasRxnCenter = new bool[n_connectedTo+1];
	for(int i=0; i<n_connectedTo; i++) {
		newConnectedTo[i]=connectedTo[i];
		newHasTraversedDownConnectedTo[i] = false;
		newOtherTemplateConnectedToIndex[i]=otherTemplateConnectedToIndex[i];
		newConnectedToHasRxnCenter[i]=connectedToHasRxnCenter[i];
	}
	newConnectedTo[n_connectedTo]=t2;
	newHasTraversedDownConnectedTo[n_connectedTo]=false;
	newOtherTemplateConnectedToIndex[n_connectedTo] = otherConToIndex;
	newConnectedToHasRxnCenter[n_connectedTo] = otherHasRxnCenter;
	delete [] connectedTo;
	delete [] hasTraversedDownConnectedTo;
	delete [] otherTemplateConnectedToIndex;
	delete [] connectedToHasRxnCenter;
	connectedTo=newConnectedTo;
	hasTraversedDownConnectedTo=newHasTraversedDownConnectedTo;
	otherTemplateConnectedToIndex=newOtherTemplateConnectedToIndex;
	connectedToHasRxnCenter = newConnectedToHasRxnCenter;
	n_connectedTo++;
}



void TemplateMolecule::printErrorAndExit(string message) {
	cerr<<"Error in TemplateMolecule!"<<endl;
	cerr<<message<<endl;
	cerr<<"\n In TemplateMolecule of type: "<< moleculeType->getName()<<endl;
	exit(4);
}


void TemplateMolecule::printDetails() {
	printDetails(cout);
}

void TemplateMolecule::printDetails(ostream &o) {
	o<<"-----------------------------\n";
	o<<"TemplateMolecule of type:   "<< moleculeType->getName();
	o<<", with id: "<<this->uniqueTemplateID;

	o<<"\n  Connected-to:                       ";
	if(n_connectedTo==0)o<<"none";
	for(int i=0;i<n_connectedTo; i++) {
		o<<"["<<connectedTo[i]->getMoleculeTypeName()<<"("<<connectedTo[i]->uniqueTemplateID<<")";
	    o<<",other index:"<<otherTemplateConnectedToIndex[i];
	    o<<","<<this->connectedToHasRxnCenter[i]<<"]   ";
	}

	o<<"\n  Empty Binding Site Constraints:      ";
	if(n_emptyComps==0)o<<"none";
	for(int i=0;i<n_emptyComps; i++)
		o<<moleculeType->getComponentName(emptyComps[i])<<"(index="<<emptyComps[i]<<")  ";

	o<<"\n  Occupied Binding Site Constraints:   ";
	if(n_occupiedComps==0)o<<"none";
	for(int i=0;i<n_occupiedComps; i++)
		o<<moleculeType->getComponentName(occupiedComps[i])<<"(index="<<occupiedComps[i]<<")  ";

	o<<"\n  State Value Constraints:             ";
	if(n_compStateConstraint==0)o<<"none";
	for(int i=0;i<n_compStateConstraint; i++) {
		o<<moleculeType->getComponentName(compStateConstraint_Comp[i])<<"(index="<<compStateConstraint_Comp[i]<<")=";
		o<<moleculeType->getComponentStateName(compStateConstraint_Comp[i],compStateConstraint_Constraint[i])<<" ";
	}

	o<<"\n  State Value Exclusions:              ";
	if(n_compStateExclusion==0)o<<"none";
	for(int i=0;i<n_compStateExclusion; i++) {
		o<<moleculeType->getComponentName(compStateExclusion_Comp[i])<<"(index="<<compStateExclusion_Exclusion[i]<<")=";
		o<<moleculeType->getComponentStateName(compStateExclusion_Comp[i],compStateExclusion_Exclusion[i])<<" ";
	}

	o<<"\n  Bond Constraints:                    ";
	if(n_bonds==0)o<<"none";
	for(int i=0;i<n_bonds; i++) {
		o<<moleculeType->getComponentName(bondComp[i])<<"(index="<<bondComp[i]<<"):";
		o<<bondPartner[i]->getMoleculeTypeName()<<"("<<bondPartnerCompName[i]<<")";
	}

	o<<"\n  Symmetric Constraints:               ";
	if(n_symComps==0)o<<"none";
	for(int i=0;i<n_symComps; i++) {

		o<<symCompName[i];
		if(symCompStateConstraint[i]>=0) {
			o<<"~"<<symCompStateConstraint[i];
		}
		if(symCompBoundState[i]==EMPTY) {
			o<<",notBonded";
		} else if (this->symCompBoundState[i]==OCCUPIED) {
			if(symBondPartner[i]==0) {
				o<<"!+";
			} else {
				o<<"!"<<symBondPartner[i]->getMoleculeTypeName();
				o<<"("<<symBondPartnerCompName[i]<<")";

			}
		} else if (symCompBoundState[i]==NO_CONSTRAINT) {
			o<<"!?";
		}
		o<<"  ";
	}

	o<<"\n  Map Generators:                      ";
	o<<n_mapGenerators<<" generators.";

	o<<endl;
	o<<endl;

}


void TemplateMolecule::addSymCompConstraint(string cName, string uniqueId,
		int bondState,int stateConstraint)
{
	if(!moleculeType->isEquivalentComponent(cName)) {
		printErrorAndExit("Cannot add a symmetric component constraint on a component that is not symmetric!");
	}

	//Error checking on bond state!!!

	//Same idea as before, we just have to do some more work for
	//symmetric constraints...  so first create some temporary arrays
	string *newSymCompName = new string[n_symComps+1];
	string *newSymCompUniqueId = new string[n_symComps+1];
	int *newSymCompStateConstraint = new int[n_symComps+1];
	int *newSymBoundState = new int[n_symComps+1];
	TemplateMolecule **newSymBondPartner = new TemplateMolecule *[n_symComps+1];
	string *newSymBondPartnerCompName = new string[n_symComps+1];
	int *newSymBondPartnerCompIndex = new int[n_symComps+1];
	bool *newHasTraversedDownSym = new bool[n_symComps+1];

	//copy over the old information
	for(int k=0; k<n_symComps; k++) {
		newSymCompName[k] = symCompName[k];
		newSymCompUniqueId[k] = symCompUniqueId[k];
		newSymCompStateConstraint[k] = symCompStateConstraint[k];
		newSymBoundState[k] = symCompBoundState[k];
		newSymBondPartner[k] = symBondPartner[k];
		newSymBondPartnerCompName[k] = symBondPartnerCompName[k];
		newSymBondPartnerCompIndex[k] = symBondPartnerCompIndex[k];
		newHasTraversedDownSym[k]=hasTraversedDownSym[k];
	}

	//store the new information
	newSymCompName[n_symComps] = cName;
	newSymCompUniqueId[n_symComps] = uniqueId;
	newSymCompStateConstraint[n_symComps] = stateConstraint;
	newSymBoundState[n_symComps] = bondState;
	newSymBondPartner[n_symComps] = 0;
	newSymBondPartnerCompName[n_symComps] = "";
	newSymBondPartnerCompIndex[n_symComps] = -1;
	newHasTraversedDownSym[n_symComps]=false;


	//delete the duplicated information
	delete [] symCompName;
	delete [] symCompUniqueId;
	delete [] symCompStateConstraint;
	delete [] symCompBoundState;
	delete [] symBondPartner;
	delete [] symBondPartnerCompName;
	delete [] symBondPartnerCompIndex;
	delete [] hasTraversedDownSym;

	//transfer over the information
	symCompName = newSymCompName;
	symCompUniqueId = newSymCompUniqueId;
	symCompStateConstraint = newSymCompStateConstraint;
	symCompBoundState = newSymBoundState;
	symBondPartner = newSymBondPartner;
	symBondPartnerCompName = newSymBondPartnerCompName;
	symBondPartnerCompIndex = newSymBondPartnerCompIndex;
	hasTraversedDownSym = newHasTraversedDownSym;

	n_symComps++;
	vector <int> canBeMappedToVector;
	canBeMappedTo.push_back(canBeMappedToVector);
}

void TemplateMolecule::addSymBond(string thisBsiteName, string thisCompId,
		TemplateMolecule *t2, string bSiteName2)
{
	//find thisCompId in our set of symmetric sites, and make changes
	//to it.  We have to register this symmetric site before we can make
	//a bond to it...
	for(int i=0; i<n_symComps; i++) {
		if(symCompUniqueId[i].compare(thisCompId)==0) {
			if(symBondPartner[i]!=0) {
				printErrorAndExit("Trying to bind a symmetric site that is already bound!");
			}
			symBondPartner[i]=t2;
			symBondPartnerCompName[i]=bSiteName2;
			if(t2->moleculeType->isEquivalentComponent(bSiteName2)) {
				symBondPartnerCompIndex[i]=-1;
			} else {
				symBondPartnerCompIndex[i]=t2->moleculeType->getCompIndexFromName(bSiteName2);
			}
			symCompBoundState[i]=OCCUPIED;
			return;
		}
	}

	//If we got here, then we couldn't find the compId in our symmetric
	//component constraints, so we must abort.
	string error="Cannot add a bond to a symmetric component if that component has not\n";
	error+="been registered yet! (meaning there is no matching uniqueComponentId.";
	printErrorAndExit(error);

}


void TemplateMolecule::addBond(string thisBsiteName,
		TemplateMolecule *t2, string bSiteName2)
{
	//If we called this, then we are adding a bond to a nonsymmetric site

	//First, initialize the new arrays
	int *newBondComp = new int[n_bonds+1];
	string *newBondCompName = new string[n_bonds+1];
	TemplateMolecule **newBondPartner = new TemplateMolecule *[n_bonds+1];
	string *newBondPartnerCompName = new string[n_bonds+1];
	int *newBondPartnerCompIndex = new int[n_bonds+1];
	bool *newHasVisitedBond = new bool[n_bonds+1];

	//Copy over original the information to the new arrays
	for(int k=0; k<n_bonds; k++) {
		newBondComp[k] = bondComp[k];
		newBondCompName[k] = bondCompName[k];
		newBondPartner[k] = bondPartner[k];
		newBondPartnerCompName[k] = bondPartnerCompName[k];
		newBondPartnerCompIndex[k] = bondPartnerCompIndex[k];
		newHasVisitedBond[k] = false; //this should always be false here!
	}

	//Insert the new information
	int compIndex=moleculeType->getCompIndexFromName(thisBsiteName);
	newBondComp[n_bonds] = compIndex;
	newBondCompName[n_bonds] = thisBsiteName;
	newBondPartner[n_bonds] = t2;
	newBondPartnerCompName[n_bonds] = bSiteName2;
	if(t2->moleculeType->isEquivalentComponent(bSiteName2)) {
		newBondPartnerCompIndex[n_bonds] = -1;
	} else {
		newBondPartnerCompIndex[n_bonds]=t2->moleculeType->getCompIndexFromName(bSiteName2);
	}
	newHasVisitedBond[n_bonds] = false;


	//Delete the duplicated information
	delete [] bondComp;
	delete [] bondCompName;
	delete [] bondPartner;
	delete [] bondPartnerCompName;
	delete [] bondPartnerCompIndex;
	delete [] hasVisitedBond;

	//Reassign the new information
	bondComp = newBondComp;
	bondCompName = newBondCompName;
	bondPartner=newBondPartner;
	bondPartnerCompName=newBondPartnerCompName;
	bondPartnerCompIndex=newBondPartnerCompIndex;
	hasVisitedBond=newHasVisitedBond;
	n_bonds++;
	compIsAlwaysMapped[compIndex]=true;
}




void TemplateMolecule::bind(TemplateMolecule *t1, string bSiteName1, string compId1,
				TemplateMolecule *t2, string bSiteName2, string compId2)
{
	if(t1->moleculeType->isEquivalentComponent(bSiteName1)) {
		t1->addSymBond(bSiteName1, compId1, t2, bSiteName2);
		//cout<<"cannot handle symmetric binding yet."<<endl;
	} else {
		t1->addBond(bSiteName1, t2, bSiteName2);
	}

	if(t2->moleculeType->isEquivalentComponent(bSiteName2)) {
		t2->addSymBond(bSiteName2, compId2, t1, bSiteName1);
		//cout<<"cannot handle symmetric binding yet."<<endl;
	} else {
		t2->addBond(bSiteName2, t1, bSiteName1);
	}
}


bool TemplateMolecule::contains(TemplateMolecule *tempMol)
{
	bool found = false;
	//the queues and lists should be static for efficiency
	//queue Q, depth queue D, and list T
	//queue <TemplateMolecule *> q;
	list <TemplateMolecule *> t;
	//queue <int> d;

	int currentDepth = 0;

	//First add this molecule
	q.push(this);
	t.push_back(this);
	d.push(currentDepth+1);
	this->hasVisitedThis=true;

	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		TemplateMolecule *cTM = q.front();
		currentDepth = d.front();
		q.pop();
		d.pop();

		//Check if we found the molecule
		if(cTM==tempMol) {
			found=true;
			break;
		}

		//Loop through the bonds for the non-symmetric case
		for(int b=0; b<cTM->n_bonds; b++) {

			//If we are connected through this bond, retrieve the connection
			if(cTM->bondPartner[b]!=0) {
				TemplateMolecule *neighbor = cTM->bondPartner[b];
				if(!neighbor->hasVisitedThis) {
					neighbor->hasVisitedThis=true;
					t.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
		}	}	}

		//Loop through the bonds for the symmetric sites, because
		//we should add those as well
		for(int b=0; b<cTM->n_symComps;b++) {
			if(cTM->symBondPartner[b]!=0) {
				TemplateMolecule *neighbor = cTM->symBondPartner[b];
				if(!neighbor->hasVisitedThis) {
					neighbor->hasVisitedThis=true;
					t.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
		}	}	}

		//Finally, also loop through the "connected-to" molecules
		for(int b=0; b<cTM->n_connectedTo;b++) {
			if(cTM->connectedTo[b]!=0) {
				TemplateMolecule *neighbor = cTM->connectedTo[b];
				if(!neighbor->hasVisitedThis) {
					neighbor->hasVisitedThis=true;
					t.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
		}	}	}
	}

	//clear the hasVisitedMolecule values
	for( tmIter = t.begin(); tmIter != t.end(); tmIter++ )
		(*tmIter)->hasVisitedThis=false;

	//the queues should be empty, but we should clear the list of
	//Template molecules that we have aggregated.
	t.clear();
	while(!q.empty()) q.pop();
	while(!d.empty()) d.pop();

	return found;
}


void TemplateMolecule::traverse(TemplateMolecule *tempMol, vector <TemplateMolecule *> &tmList, bool skipConnectedTo)
{
	//cout<<"traversing"<<endl;
	//the queues and lists should be static for efficiency
	//queue Q, depth queue D
	//queue <TemplateMolecule *> q;
	//queue <int> d;
	int currentDepth = 0;

	q.push(tempMol);
	tmList.push_back(tempMol);
	d.push(currentDepth+1);
	tempMol->hasVisitedThis=true;

	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		TemplateMolecule *cTM = q.front();
		currentDepth = d.front();
		q.pop();
		d.pop();

		//Loop through the bonds for the non-symmetric case
		for(int b=0; b<cTM->n_bonds; b++) {

		//If we are connected through this bond, retrieve the connection
		if(cTM->bondPartner[b]!=0) {
			TemplateMolecule *neighbor = cTM->bondPartner[b];
			if(!neighbor->hasVisitedThis) {
				neighbor->hasVisitedThis=true;
				tmList.push_back(neighbor);
				q.push(neighbor);
				d.push(currentDepth+1);
		}	}	}

		//Loop through the bonds for the symmetric sites, because
		//we should add those as well
		for(int b=0; b<cTM->n_symComps;b++) {
			if(cTM->symBondPartner[b]!=0) {
				TemplateMolecule *neighbor = cTM->symBondPartner[b];
				if(!neighbor->hasVisitedThis) {
					neighbor->hasVisitedThis=true;
					tmList.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
		}	}	}

		if(!skipConnectedTo) {
			//Finally, also loop through the "connected-to" molecules
			for(int b=0; b<cTM->n_connectedTo;b++) {
				if(cTM->connectedTo[b]!=0) {
					TemplateMolecule *neighbor = cTM->connectedTo[b];
					if(!neighbor->hasVisitedThis) {
						neighbor->hasVisitedThis=true;
						tmList.push_back(neighbor);
						q.push(neighbor);
						d.push(currentDepth+1);
			}	}	}
		}
	}

	//clear the has visitedMolecule values
	for( tmVecIter = tmList.begin(); tmVecIter != tmList.end(); tmVecIter++ ) {
		(*tmVecIter)->hasVisitedThis=false;
	}

	while(!q.empty()) q.pop();
	while(!d.empty()) d.pop();
}



void TemplateMolecule::clear() {
	if(matchMolecule!=0) {
		matchMolecule->isMatchedTo=0;
		matchMolecule=0;
	}
	for(int b=0; b<n_bonds; b++) hasVisitedBond[b]=false;

	if(n_symComps>0) {
		for(int c=0; c<n_symComps; c++) {
			canBeMappedTo.at(c).clear();
		}
	}

	for(int t=0; t<n_connectedTo; t++) {
		hasTraversedDownConnectedTo[t]=false;
	}

	//cout<<"clear"<<endl;
}


bool TemplateMolecule::compare(Molecule *m)
{
	return compare(m,0,0);
}


bool TemplateMolecule::tryToMap(Molecule *toMap, string toMapComponent,
		Molecule *mappedFrom, string mappedFromComponent)
{
	//cout<<"trying to map: "<<toMapComponent<<endl;
	//cout<<"to me: "<<endl;
	//this->printDetails(cout);
	bool canMapToSomething=false;
	for(int c=0; c<n_symComps; c++)
	{
		//cout<<"analyzing sym constraint: "<< this->symCompName[c]<<endl;
		if(symCompName[c].compare(toMapComponent)!=0) continue;

		int *molEqComp; int n_molEqComp=0;
		moleculeType->getEquivalencyClass(molEqComp,n_molEqComp,this->symCompName[c]);
		for(int sc=0; sc<n_molEqComp; sc++)
		{
			//cout<<"  comparing to site: "<<this->moleculeType->getComponentName(molEqComp[sc])<<endl;
			//first make sure that we can map to this component
			if(compIsAlwaysMapped[molEqComp[sc]]) continue;

			//Now check each constraint to see if things are bound correctly
			if(this->symCompBoundState[c]==TemplateMolecule::EMPTY) {
				if(!toMap->isBindingSiteOpen(molEqComp[sc])) continue;
			} else if (symCompBoundState[c]==TemplateMolecule::OCCUPIED) {
				if(!toMap->isBindingSiteBonded(molEqComp[sc])) continue;
			}

			//cout<<"sites are properly set up."<<endl;

			//make sure the states match up
			if(symCompStateConstraint[c]!=TemplateMolecule::NO_CONSTRAINT) {
				if(toMap->getComponentState(molEqComp[sc])!=this->symCompStateConstraint[c]) {
					continue;
				}
			}

			//cout<<"state constraints are fine."<<endl;

			if(symBondPartner[c]!=0) {
				Molecule *m2=toMap->getBondedMolecule(molEqComp[sc]);
				if(m2!=0) {
					if(m2!=mappedFrom) {
						continue;
					}
				} else { continue; }
			} else {
				continue;
			}

			canMapToSomething=true;
			canBeMappedTo.at(c).push_back(molEqComp[sc]);
			//cout<<"Looks like we can map."<<endl;
		}


	}
	//If we couldn't map this symmetric component, then we must quit
	if(!canMapToSomething) {
		clear(); return false;
	}


	return true;
}

int TemplateMolecule::getNumDisjointSets(vector < TemplateMolecule * > &tMolecules,
				vector <vector <TemplateMolecule *> > &sets,
				vector <int> &uniqueSetId)
{
	int setCount=0;
	for(unsigned int i=0; i<tMolecules.size(); i++)
	{
		//First see if this template was already found in a previous set.
		//if it was, then we don't have to traverse
		bool alreadyFound = false;
		for(unsigned int j=0; j<i; j++) {
			//search set J for this template
			for(unsigned int kj=0; kj<sets.at(j).size(); kj++) {
				if(sets.at(j).at(kj)==tMolecules.at(i)) {
					alreadyFound = true;
					break;
				}
			}
			//If we found it, remember the uniqueSetId of set J
			if(alreadyFound) {
				uniqueSetId.push_back(uniqueSetId.at(j));
				vector <TemplateMolecule *> thisSet;
				sets.push_back(thisSet);
				break;
			}
		}
		if(alreadyFound) { continue; }
		else {
			//If we have not found this molecule before, then
			//it must be in a new set, so we traverse and remember that set

			uniqueSetId.push_back(setCount);
			setCount++;
			vector <TemplateMolecule *> thisSet;
			TemplateMolecule::traverse(tMolecules.at(i),thisSet,TemplateMolecule::SKIP_CONNECTED_TO);
			sets.push_back(thisSet);
		}
	}
	return setCount;
}

bool TemplateMolecule::isSymMapValid()
{
	//we have to find if a non overlapping set of mappings exists
	//for instance, if symComp1 and symComp2 can only be mapped to
	//the same site, then we do not have a match, because we have to
	//map both.

	//KEY ASSUMPTION!!  This appears to be true:
	// if c1 can be mapped to pos 0,1,2
	// and c2 can be mapped to pos 0, then c1 and c2
	// must be identical and therefore c2 must also map
	// to c1.
	// In other words, if we assume that c1 and c2 can only map to ALL the
	// same sites, or completely disjoint sets, then the method below
	// works (which I think is reasonable).  If not, then the method below
	// will miss certain permutations and we will get into trouble
	// certain

	//Print the entire canBeMappedTo arrays for debugging
	//cout<<"checking if sym map is valid. "<<endl;
	//for(int i=0; i<n_symComps; i++) {
	//	cout<<"comp: "<<i<<" ("<<this->symCompName[i]<<"): ";
	//	for(unsigned int k=0; k<canBeMappedTo.at(i).size();k++) {
	//		cout<<"  "<<canBeMappedTo.at(i).at(k);
	//	} cout<<endl;
	//}

	//Start at the current position equal to zero
	int * curPos = new int[this->n_symComps];
	curPos[0]=0;
	//for(int i=0; i<n_symComps; i++) { curPos[i]=0; }


	int startComp=0; //Arbitrarily map the first component to its first posibility
	int counter=0;
	for(int curComp=startComp+1; curComp<n_symComps; curComp++)
	{
		//try to find a match for this current component
		bool foundValidComp=false;
		for(unsigned int i=0; i<canBeMappedTo.at(curComp).size(); i++) {
			bool isOkPos=true;
			for(int k=0; k<curComp;k++) {
				if(canBeMappedTo.at(k).at(curPos[k]) == canBeMappedTo.at(curComp).at(i)) {
					isOkPos=false; break;
				}
			}
			if(isOkPos) {
				curPos[curComp]=i;
				foundValidComp=true;
				break;
			}
		}

		//Debugging: this prints out the current mapping that is possible
		//for(int k=0; k<curComp; k++) {
		//		if(k==0) cout<<"\t"<<counter<<": [";
		//		else cout<<"[";
		//		cout<<canBeMappedTo.at(k).at(curPos[k])<<"]";
		//} cout<<endl;

		if(!foundValidComp) {
			//cout<<"nope.  sym map is not valid."<<endl;
			return false;
			//by the above assumption, if we cannot find a permutation that matches
			//this site, then no valid permutation can exist ever! So quit immediately
		}

		//cout<<"this works so far: ";
		//for(int k=0; k<curComp; k++) {
		//		if(k==0) cout<<"\t"<<counter<<": [";
		//		else cout<<"[";
		//		cout<<canBeMappedTo.at(k).at(curPos[k])<<"]";
		//} cout<<endl;

		counter++;
	}

	//If we got here, then we matched every site succesfully, so we
	//can generate a unique mapping.  This symMap is valid.
	//cout<<"yes, it is valid."<<endl;
	return true;



//  Here is code for a full traversal, which is usually fast, but can be
//  exceptionally long in certain circumstances, (for instance, if there are 6 identical
//  sites, then this will consider 6^6 permutations), which is why we make the above
//  assumption given in the beginning....
//	//Loop until we've found a valid permutation, or have exhausted all possibilities
//	int counter = 0;
//	bool foundValid=false;
//	bool lastCycle=false;
//	while(!foundValid)
//	{
//		//Print current permutation (for debugging)
//		//for(int k=0; k<n_symComps; k++) {
//		//	if(k==0) cout<<"\t"<<counter<<": [";
//		//	else cout<<"[";
//		//	cout<<canBeMappedTo.at(k).at(curPos[k])<<"]";
//		//}
//
//		//Check if it is valid (note, might be more efficient to sort, then check
//		//neighbors - but if there is only two or three unique sites, this will be
//		//just as fast, because there is no overhead).
//		bool isValid=true;
//		for(int j=1;j<n_symComps; j++) {
//			for(int k=0; k<j; k++) {
//				if(canBeMappedTo.at(k).at(curPos[k]) == canBeMappedTo.at(j).at(curPos[j])) {
//					isValid=false; break;
//				}
//			} if(!isValid) break;
//		}
//
//		//If it is valid, then we can just return true, because we are ok
//		cout<<" is valid? ";
//		if(isValid) { cout<<"yes"; return true; } else { cout<<"no"; }
//		cout<<endl;
//		//
//
//		counter++;
//
//		//Cycle to the next permutation
//		if(lastCycle) break;
//
//		int curComponent = n_symComps-1;
//		do {
//			curPos[curComponent]++;
//			if(curPos[curComponent]>=(int)canBeMappedTo.at(curComponent).size()) {
//				curPos[curComponent] = 0;
//				curComponent--;
//			} else { break; }
//			if(curComponent<0) {
//				lastCycle = true; break;
//			}
//		} while(true);
//	}
//  return true;
}

bool TemplateMolecule::compare(Molecule *m, ReactantContainer *rc, MappingSet *ms)
{
	//this->printDetails();
	//cout<<"comparing to: "<<endl;
	//m->printDetails();
	//if(m->isMatchedTo!=0) {
	//	if(m->isMatchedTo==this) {return true; }
	//	else { clear(); return false; }
	//}

	//First check if we've been here before, and return accordingly
	if(this->matchMolecule!=0) {
		if(matchMolecule==m) { return true; }
		else { clear(); return false; }
	}

	//Make sure we are of the same type
	if(m->getMoleculeType()!=this->moleculeType) {
		clear(); return false;
	}

	//Check all the basic components first to get them out of the way
	//First check that all of our states match
	for(int c=0; c<n_compStateConstraint; c++) {
		if(m->getComponentState(compStateConstraint_Comp[c]) != compStateConstraint_Constraint[c]) {
			clear(); return false;
		}
	}
	//Check that all of our exclusions are indeed not present (for state!=value checks)
	for(int c=0; c<n_compStateExclusion; c++) {
		if(m->getComponentState(compStateExclusion_Comp[c]) == compStateExclusion_Exclusion[c]) {
			clear(); return false;
		}
	}
	//Make sure binding sites that are open / occupied are
	for(int c=0; c<n_emptyComps; c++) {
		if(!m->isBindingSiteOpen(emptyComps[c])) {
			clear(); return false;
		}
	}
	for(int c=0; c<n_occupiedComps; c++) {
		if(!m->isBindingSiteBonded(occupiedComps[c])) {
			clear(); return false;
		}
	}


	//cout<<"all the basic things match."<<endl;
	//Good, good - everything matches so let's set our match molecule
	matchMolecule = m;
	m->isMatchedTo=this;


	//Now for the tricky and fun part.  The actual traversal....
	//Cycle through the bonds
	for(int b=0; b<n_bonds; b++)
	{
		//continue on if we've been down this bond before
		if(hasVisitedBond[b]) { continue; }

		//The binding site must be occupied!
		if(m->isBindingSiteOpen(bondComp[b])) {
			clear(); return false;
		}

		//Grab the template molecule and the actual molecule that we have to compare
		TemplateMolecule *t2=bondPartner[b];
		Molecule *m2=m->getBondedMolecule(bondComp[b]);

		//If this template has been matched already, then it should be matched
		//to m2.  IF not, we have problems.  If it has already been matched, then continue.
		if(t2->matchMolecule!=0) {
			if(t2->matchMolecule!=m2) {
				//cout<<"match molecules do not match"<<endl;
				clear(); return false;
			} else { continue; }
		}

		//Now we have to make sure that the correct bond is present on the other side
		//of the interaction. This can be tricky because of symmetric sites!
		if(bondPartnerCompIndex[b]<0) { //Argh!  this means we are symmetric!!

			//See if we can map this bond to one of the symmetric components on m2
			bool canMap = t2->tryToMap(m2,bondPartnerCompName[b],m,bondCompName[b]);
			if(canMap) {
				//cout<<"from non sym, can map, going down sym site"<<endl;
				this->hasVisitedBond[b]=true;
				bool match=t2->compare(m2,rc,ms);
				if(!match) {
					clear(); return false;
				}
			} else {
			//	cout<<"cannot map sym site.."<<endl;
				clear(); return false;
			}

		} else { //Phew!  we can check this guy normally.

			//First, make sure the opposite molecule is connected to this molecule correctly
			//and has not been matched before.
			Molecule *potentialMatch=m2->getBondedMolecule(bondPartnerCompIndex[b]);
			if(potentialMatch==Molecule::NOBOND) {
				clear(); return false;
			} else if(potentialMatch!=matchMolecule) {
				clear(); return false;
			}

			//Remember that we've visited this bond before, should we ever come back to it.
			this->hasVisitedBond[b]=true;

			//Now traverse onto this molecule, and make sure we match down the list
			bool match=t2->compare(m2,rc,ms);
			if(!match) {
				clear(); return false;
			}
		}

	}


	//cout<<"non-symmetric bonds match"<<endl;


	//////////////////////////////////////////////////////////////////////////
	//Now go through each of the symmetric sites and try to map them
	for(int c=0; c<n_symComps; c++)
	{
		//cout<<"comparing symComp: "<<symCompName[c]<<endl;
		//Loop through each of the equivalent components to see if we can match them
		int *molEqComp; int n_molEqComp=0;
		moleculeType->getEquivalencyClass(molEqComp,n_molEqComp,this->symCompName[c]);
		for(int sc=0; sc<n_molEqComp; sc++)
		{
			//cout<<"  -to molecule comp: "<<m->getMoleculeType()->getComponentName(molEqComp[sc])<<endl;

			//first make sure that we can map to this component
			if(compIsAlwaysMapped[molEqComp[sc]]) continue;

			//Now check each constraint to see if things are bound correctly
			if(this->symCompBoundState[c]==TemplateMolecule::EMPTY) {
				if(!m->isBindingSiteOpen(molEqComp[sc])) continue;
			} else if (symCompBoundState[c]==TemplateMolecule::OCCUPIED) {
				if(!m->isBindingSiteBonded(molEqComp[sc])) continue;
			}

			//make sure the states match up
			if(symCompStateConstraint[c]!=TemplateMolecule::NO_CONSTRAINT) {
				if(m->getComponentState(molEqComp[sc])!=this->symCompStateConstraint[c]) {
					continue;
				}
			}

			//cout<<"  -basic bound states match"<<endl;

			//Now make sure binding sites match up.  This can be tricky!
			if(symBondPartner[c]!=0) {
				//cout<<"checking if bond partner matches."<<endl;


				//Grab the template molecule and the actual molecule that we have to compare
				TemplateMolecule *t2=symBondPartner[c];
				Molecule *m2=m->getBondedMolecule(molEqComp[sc]);

				//Check for a matched molecule, if we are already connected to something
				//that has been matched, then we have already compared this guy from the
				//other end (or we will...)
				if(t2->matchMolecule!=0) {
					if(t2->matchMolecule!=m2) {
						continue;
					} else {
						//this->canBeMappedTo.at(c).push_back(molEqComp[sc]);
						// don't need to do this here, because the other end already did.
						continue;
					}
				}


				if(symBondPartnerCompIndex[c]<0) { //Argh!  this means we are symmetric!!
					//See if we can map this bond to one of the symmetric components on m2
					bool canMap = t2->tryToMap(m2,symBondPartnerCompName[c],m,this->symCompName[c]);
					if(canMap) {
						//cout<<"comparing down symmetric site.."<<endl;
						bool match=t2->compare(m2,rc,ms);
						if(!match) { continue; } //keep going if we can't match
					}
				} else { //Phew!  we can check this guy normally.

					//First, make sure the opposite molecule is connected to this molecule correctly
					Molecule *potentialMatch=m2->getBondedMolecule(symBondPartnerCompIndex[c]);
					if(potentialMatch==Molecule::NOBOND) {
						continue;
					} else if(potentialMatch!=matchMolecule) {
						continue;
					}

					//Now traverse onto this molecule, and make sure we match down the list
					bool match=t2->compare(m2,rc,ms);
					if(!match) { continue; }
				}
			}

			//if we got here, then by golly, I think we got a match!  So remember it!
			canBeMappedTo.at(c).push_back(molEqComp[sc]);
		}

		//If we couldn't map this symmetric component, then we must quit
		if(canBeMappedTo.at(c).size()==0) {
			//cout<<"could not find a mapping. (canMapThisComponent=false)"<<endl;
			clear(); return false;
		}
		//else {
			//cout<<"can be mapped"<<endl;
		//}
	}

	//Great, if we got here, everything matched up, all components can be mapped, and
	//we just have to double check if our mappings are valid...
	if(this->n_symComps>1) {
		if(!isSymMapValid()) {
			//oh no!  we were so close, but in the end, we couldn't get a unique
			//mapping onto all of the identical components
			clear(); return false;
		}
	}



	//If we were given a mappingSet, then map this molecule
	//with all the generators we've got (NOTE that we have to do this BEFORE we
	//look at connected molecules, so that when we clone, we also clone maps onto
	//this molecule
	if(ms!=0) {
		for(int i=0;i<n_mapGenerators; i++) {
			mapGenerators[i]->map(ms,m);
		}
	}


	//Check connected-to molecules
	if(n_connectedTo>0) {

		vector <MappingSet *> lastMappingSets;
		lastMappingSets.push_back(ms);

		for(int cTo=0; cTo<this->n_connectedTo; cTo++) {

			if(hasTraversedDownConnectedTo[cTo]) continue;

			list <Molecule *> molList;
			list <Molecule *>::iterator molIter;
			m->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);

			bool canMatch=false;
			for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
				if((*molIter)->isMatchedTo!=0) continue;

				//remember that we went down this route before, so we don't just go back and forth
				//between connectedTo bonds...
				bool canMatchThis = false;
				connectedTo[cTo]->hasTraversedDownConnectedTo[otherTemplateConnectedToIndex[cTo]]=true;

				//If the other set has a reaction center, then the other set
				//is merely context, so we only have to find a single instance of it
				//and return as soon as we have matched.
				if(!connectedToHasRxnCenter[cTo]) {
					canMatchThis=connectedTo[cTo]->compare((*molIter));
					if(canMatchThis) { canMatch=true; continue; }
				}

				// otherwise, we have to do more work...
				else {
					canMatchThis=connectedTo[cTo]->compare((*molIter),rc,lastMappingSets.at(lastMappingSets.size()-1));
					if(canMatchThis) {
						rc->notifyPresenceOfClonedMappings();
						canMatch = true;
						MappingSet * newMS = rc->pushNextAvailableMappingSet();
						MappingSet::clone(lastMappingSets.at(lastMappingSets.size()-1),newMS);
						lastMappingSets.push_back(newMS);
					}
				}
			}
			if(!canMatch) {
				//cout<<"could not match this!"<<endl;
				//Clear out all mapping sets that we tried to assign.  We can do this
				//easily by just clearing the clone pointer from the original mappingset
				//
				//ms->clearClonedMapping();
				//if(lastMappingSets.size()>1)
				//	rc->removeMappingSet(lastMappingSets.at(1)->getId());
				clear(); return false;
			}

		}
		//if(lastMappingSets.size()>2) {
		//	rc->notifyPresenceOfClonedMappings();
		//}

		if(lastMappingSets.size()>1) {
			lastMappingSets.at(lastMappingSets.size()-2)->clearClonedMapping();
			rc->removeMappingSet(lastMappingSets.at(lastMappingSets.size()-1)->getId());
		}
	}
	///  End handle connected-to

	clear();
	return true;
}


