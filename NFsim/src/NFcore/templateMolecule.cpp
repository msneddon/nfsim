#include <iostream>

#include "NFcore.hh"

using namespace std;
using namespace NFcore;


TemplateMolecule::TemplateMolecule(MoleculeType * parentMoleculeType)
{
	if(DEBUG) cout<<"-creating template molecule of type " << parentMoleculeType->getName() << endl;
	this->parentMoleculeType = parentMoleculeType;
	this->parentMoleculeType->addTemplateMolecule(this);
	clear();
}

TemplateMolecule::~TemplateMolecule()
{
	if(DEBUG) cout <<"   -destroying template molecule of type " << parentMoleculeType->getName() << endl;
	
	
	MapGenerator *mg;
	while(mapGenerators.size()>0) {
		mg = mapGenerators.back();
		mapGenerators.pop_back();
		delete mg;
	}
}

void TemplateMolecule::setHasVisited(int bSiteCompareIndex)
{
	hasVisitedBond.at(bSiteCompareIndex) = true;
}

//void TemplateMolecule::addBindingSite(int bSiteIndex)
//{
//	// get the binding site from the MoleculeType
//	this->bSiteIndex.push_back(bSiteIndex); // saves index into bSite array in MoleculeType
////	TemplateSite * b = new TemplateSite(parentMoleculeType->getBindingSiteName(bSiteIndex),this, bonds.size());
//	bonds.push_back(b);
//	hasVisitedBond.push_back(false);
//}


void TemplateMolecule::addStateValue(string cName, int stateValue)
{
	this->addStateValue(parentMoleculeType->getCompIndexFromName(cName), stateValue);
}

void TemplateMolecule::addStateValue(string cName, string stateValue)
{
	int cIndex = parentMoleculeType->getCompIndexFromName(cName);
	this->addStateValue(cIndex, parentMoleculeType->getStateValueFromName(cIndex,stateValue));
}

void TemplateMolecule::addStateValue(int cIndex, int stateValue)
{
	this->stateIndex.push_back(cIndex);
	this->stateValue.push_back(stateValue);
}

void TemplateMolecule::addNotStateValue(char * stateName, int notStateValue)
{
	this->notStateIndex.push_back(parentMoleculeType->getCompIndexFromName(stateName));
	this->notStateValue.push_back(notStateValue);
}



bool TemplateMolecule::isBindingSiteOpen(int bIndex) const
{
	if(bonds.at(bIndex)==0) return true;
	return false;
}

bool TemplateMolecule::isBindingSiteBonded(int bIndex) const
{
	if(bonds.at(bIndex)==0) return false;
	return true;	
}

TemplateMolecule * TemplateMolecule::getBondedTemplateMolecule(int bIndex) const
{
	return bonds.at(bIndex);	
}

int TemplateMolecule::getBindingSiteIndex(int bIndex) const
{
	return bSiteIndex.at(bIndex);
}
int TemplateMolecule::getBindingSiteIndexOfBondedTemplate(int bIndex) const
{
	return bSiteIndexOfBond.at(bIndex);
}

unsigned int TemplateMolecule::addEmptyBindingSite(const char * bSiteName)
{
	return this->addEmptyBindingSite(parentMoleculeType->getCompIndexFromName(bSiteName));	
}

unsigned int TemplateMolecule::addEmptyBindingSite(int bSiteIndex)
{
	this->bSiteIndex.push_back(bSiteIndex);
	bonds.push_back(0);
	hasVisitedBond.push_back(false);
	bSiteIndexOfBond.push_back(-1);
	return (this->bSiteIndex.size()-1);  //return the index that we added this bSiteIndex to
}


void TemplateMolecule::addOccupiedBindingSite(const char * bSiteName)
{
	this->sitesThatMustBeOccupied.push_back(parentMoleculeType->getCompIndexFromName(bSiteName));
}


bool TemplateMolecule::isStateValue(const char * stateName, int stateVal)
{
	int sIndex = parentMoleculeType->getCompIndexFromName(stateName);
	for(unsigned int s=0; s<stateIndex.size(); s++)
	{
			if(stateIndex.at(s)==sIndex)
			{ 
				if(stateValue.at(s)==stateVal)
					return true;
				return false;
			}
	}
	return false;
}

bool TemplateMolecule::isBonded(const char * bSiteName)
{
	int bIndex = parentMoleculeType->getCompIndexFromName(bSiteName);
	for(unsigned int b=0; b<sitesThatMustBeOccupied.size(); b++)
	{
		if(sitesThatMustBeOccupied.at(b)==bIndex)
			return true;
	}
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		if(bSiteIndex.at(b)==bIndex)
		{
			if(bonds.at(b)!=NULL)
				return true;
			else
				return false;
		}
	}
	return false;
	
}

void TemplateMolecule::bind(TemplateMolecule *t1, int bSiteIndex1, TemplateMolecule *t2, int bSiteIndex2)
{
	//Find the site of bSiteIndex1 and bSiteIndex2
	int bIndex1 = -1;
	for(unsigned int b1=0; b1<t1->bSiteIndex.size(); b1++)
	{
		if(t1->bSiteIndex.at(b1)==bSiteIndex1) bIndex1 = b1;	
	}
	int bIndex2 = -1;
	for(unsigned int b2=0; b2<t2->bSiteIndex.size(); b2++)
	{
		if(t2->bSiteIndex.at(b2)==bSiteIndex2) bIndex2 = b2;	
	}
	
	//IF we didn't have the bonds created already, create them now
	if(bIndex1==-1) bIndex1 = t1->addEmptyBindingSite(bSiteIndex1);
	if(bIndex2==-1) bIndex2 = t2->addEmptyBindingSite(bSiteIndex2);
	
	//Finally, we can set the bonds
	t1->bonds.at(bIndex1) = t2;
	t1->bSiteIndexOfBond.at(bIndex1) = bIndex2;
	t2->bonds.at(bIndex2) = t1;
	t2->bSiteIndexOfBond.at(bIndex2) = bIndex1;
}

void TemplateMolecule::bind(TemplateMolecule *t1, const char * bSiteName1, TemplateMolecule *t2, const char * bSiteName2)
{
	int bSiteIndex1 = t1->parentMoleculeType->getCompIndexFromName(bSiteName1);
	int bSiteIndex2 = t2->parentMoleculeType->getCompIndexFromName(bSiteName2);
	
	bind(t1,bSiteIndex1, t2, bSiteIndex2); 
}


queue <TemplateMolecule *> TemplateMolecule::tmq;
queue <Molecule *> TemplateMolecule::mq;
list <TemplateMolecule *> TemplateMolecule::tml;
queue <int> TemplateMolecule::d;
list <TemplateMolecule *>::iterator TemplateMolecule::tmIter;

bool TemplateMolecule::compareBreadthFirst(TemplateMolecule *tm, Molecule *m)
{
	bool match=true;
	
	//Create the queues and lists
	//queue <TemplateMolecule *> tmq;
	//queue <Molecule *> mq;
	//list <TemplateMolecule *> tml;
	//queue <int> d;
	int currentDepth = 0;
	
	tmq.push(tm);
	mq.push(m);
	tml.push_back(tm);
	d.push(currentDepth+1);
	tm->hasVisited=true;

		
	//Look at children until the queue is empty
	while(!tmq.empty())
	{
		//Get the next parent to look at (currentMolecule)
		TemplateMolecule *cTM = tmq.front();
		Molecule *cM = mq.front();
		currentDepth = d.front();
		tmq.pop(); 
		mq.pop();
		d.pop();
			
		//First check if we've been here before
		if(cTM->matchMolecule!=0) {
			if(cTM->matchMolecule == m) { continue; }
			else { match=false; break; }
		}
			
		//Next, check if we have the same type
		if(cM->getMoleculeType()!=cTM->parentMoleculeType){
			match=false; break;
		}
		
		//Check that states are correct
		for(unsigned int s=0; s<cTM->stateIndex.size(); s++) {
			if(cM->getComponentState(cTM->stateIndex.at(s)) != cTM->stateValue.at(s)) { 
					match=false; break;
			}
		}
			
		for(unsigned int s=0; s<cTM->notStateIndex.size(); s++) {
			if(cM->getComponentState(cTM->notStateIndex.at(s)) == cTM->notStateValue.at(s)) { 
					match=false; break;
			}
		}
			
			
		//Check if the sites that must be occupied by anything are occupied
		for(unsigned int b=0; b<cTM->sitesThatMustBeOccupied.size(); b++)
		{
			//cout<<"Checking site that must be occupied: "<< sitesThatMustBeOccupied.at(b) <<endl;
			//If the site is not occupied in the molecule, no match
			if(!cM->isBindingSiteBonded(cTM->sitesThatMustBeOccupied.at(b)))
			{
				match=false; break;
			}
		}
		
		//We have matched so far, so set our match
		cTM->matchMolecule = cM;
		
		

		
		//Loop through the bonds to get the next set of Templates to search
		int bMax = cTM->getNumBindingSites();
		for(int b=0; b<bMax; b++)
		{
			if(cTM->bonds.at(b)==0) //template binding site is open, so molecule binding site must be too.
			{
				if( cM->isBindingSiteBonded(cTM->bSiteIndex.at(b))) {
					match=false; break;
				}
			}
			
			
			if(cTM->isBindingSiteBonded(b)) //the template site is bonded, 
			{
				if(cM->isBindingSiteOpen(cTM->bSiteIndex.at(b))) { //so the molecule site must be bonded too
					match=false; break;
				}
				
				//Grab the bonded neighbors of both the template and the molecule
				TemplateMolecule *tempNeighbor = cTM->getBondedTemplateMolecule(b);
				Molecule *molNeighbor = cM->getBondedMolecule(cTM->bSiteIndex.at(b));
				
				//Make sure the back binding (of molNeighbor to cM is through the correct site on molNeighbor)
				if(cM->getComponentIndexOfBond(cTM->bSiteIndex.at(b)) != tempNeighbor->bSiteIndex.at(cTM->bSiteIndexOfBond.at(b))) {
					match=false; break;
				}
				
				//Looks like everything is good to go, so add this guy to the list
				if(!tempNeighbor->hasVisited)
				{
					tempNeighbor->hasVisited=true;
					tml.push_back(tempNeighbor);
					tmq.push(tempNeighbor);
					mq.push(molNeighbor);
					d.push(currentDepth+1);
				}
			}
		}
		if(match=false) break;
	}
	
	
	
	
	//clear the has visitedMolecule values
	//list <TemplateMolecule *>::iterator tmIter;
	for( tmIter = tml.begin(); tmIter != tml.end(); tmIter++ ) {
  		(*tmIter)->hasVisited=false;
  		(*tmIter)->matchMolecule=0;
	}
	
	
	while(!tmq.empty()) tmq.pop();
	while(!mq.empty()) mq.pop();
	while(!d.empty()) d.pop();
	tml.clear();
	return match;
}


bool TemplateMolecule::compare(Molecule * m)
{
	//return compareBreadthFirst(this, m);
	//if(DEBUG) {
	//	cout<<"Calling Template Molecule compare: "<<endl;
	//	m->printDetails();
	//	this->printDetails();
	//}
	
	
	//First check if we've been here before
	if(matchMolecule!=0)
	{
		if(matchMolecule == m) { return true;	}
		else { clear(); return false; }
	}
	
	//Next, check if we have the same type
	if(m->getMoleculeType()!=this->parentMoleculeType)
	{
		clear();
		return false;
	}
	//Check that states are correct
	for(unsigned int s=0; s<stateIndex.size(); s++)
	{
		if(m->getComponentState(stateIndex.at(s)) != stateValue.at(s))
		{ 
			clear();
			return false; 
		}
	}
	
	for(unsigned int s=0; s<notStateIndex.size(); s++)
	{
		if(m->getComponentState(notStateIndex.at(s)) == notStateValue.at(s))
		{ 
			clear();
			return false; 
		}
	}
	
	
	//Check if the sites that must be occupied by anything are occupied
	for(unsigned int b=0; b<sitesThatMustBeOccupied.size(); b++)
	{
		//cout<<"Checking site that must be occupied: "<< sitesThatMustBeOccupied.at(b) <<endl;
		//If the site is not occupied in the molecule, no match
		if(!m->isBindingSiteBonded(sitesThatMustBeOccupied.at(b)))
		{
			clear();
			return false;
		}
	}
	
	
	
	//We have matched so far, so set our match
	matchMolecule = m;
	
	//Cycle through the bonds
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		if(hasVisitedBond.at(b)==true) continue;
		if(bonds.at(b)==0) //template binding site is open, so molecule binding site must be too.
		{
			if( m->isBindingSiteBonded(bSiteIndex.at(b)))
			{
				clear();
				return false;
			}
		}
		else //template binding site has a bond
		{
			if(m->isBindingSiteOpen(bSiteIndex.at(b)))
			{
				clear();
				return false;
			}
			
			//get template that is bound to this binding site
			TemplateMolecule * t2 = this->getBondedTemplateMolecule(b);
			Molecule * m2 = m->getBondedMolecule(bSiteIndex.at(b));
			
			//We have to remember to check that it was the right bond site on the other molecule
			//too.  This means that the bond site in m2 that m is connected to must be the same
			//bond site in t2 that t1 is connected to.  If not, then we must say false
			if(m->getComponentIndexOfBond(bSiteIndex.at(b)) != t2->bSiteIndex.at(bSiteIndexOfBond.at(b)))
			{
				clear();
				return false;
			}
			
			//tell t2 we have visited this bond before, now that we are sure it is correct
			t2->setHasVisited(bSiteIndexOfBond.at(b));
			
			//remember that we visited through this bond in this template, in case we get to
			//this molecule again
			this->setHasVisited(b); 
			
			bool goodMatch = t2->compare(m2);
			if(!goodMatch)
			{
				clear();
				return false;
			}
		}		
	}
	
	// we have a match!  so clear and return
	//cout<<"Match!"<<endl;
	clear();
	return true;
}


//void TemplateMolecule::addTemplateMapping(TemplateMapping *tm)
//{
//	tMappings.push_back(tm);
//}
/*
bool TemplateMolecule::compare(Molecule * m, MappingSet *mappingSet)
{	//First check if we've been here before
	if(matchMolecule!=0)
	{
		if(matchMolecule == m) { return true;	}
		else { clear(); return false; }
	}
	
	//Next, check if we have the same type
	if(m->getMoleculeType()!=this->parentMoleculeType)
	{
		clear();
		return false;
	}
	//Check that states are correct
	for(unsigned int s=0; s<stateIndex.size(); s++)
	{
		if(m->getState(stateIndex.at(s)) != stateValue.at(s))
		{ 
			clear();
			return false; 
		}
	}
	
	for(unsigned int s=0; s<notStateIndex.size(); s++)
	{
		if(m->getState(notStateIndex.at(s)) == notStateValue.at(s))
		{ 
			clear();
			return false; 
		}
	}
	
	
	//Check if the sites that must be occupied by anything are occupied
	for(unsigned int b=0; b<sitesThatMustBeOccupied.size(); b++)
	{
		//cout<<"Checking site that must be occupied: "<< sitesThatMustBeOccupied.at(b) <<endl;
		//If the site is not occupied in the molecule, no match
		if(!m->isBindingSiteBonded(sitesThatMustBeOccupied.at(b)))
		{
			clear();
			return false;
		}
	}
	
	
	
	//We have matched so far, so set our match
	matchMolecule = m;
	
	//Cycle through the bonds
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		if(hasVisitedBond.at(b)==true) continue;
		if(bonds.at(b)==0) //template binding site is open, so molecule binding site must be too.
		{
			if( m->isBindingSiteBonded(bSiteIndex.at(b)))
			{
				clear();
				return false;
			}
		}
		else //template binding site has a bond
		{
			if(m->isBindingSiteOpen(bSiteIndex.at(b)))
			{
				clear();
				return false;
			}
			
			//get template that is bound to this binding site
			TemplateMolecule * t2 = this->getBondedTemplateMolecule(b);
			Molecule * m2 = m->getBondedMolecule(bSiteIndex.at(b));
			
			//We have to remember to check that it was the right bond site on the other molecule
			//too.  This means that the bond site in m2 that m is connected to must be the same
			//bond site in t2 that t1 is connected to.  If not, then we must say false
			if(m->getBsiteIndexOfBond(bSiteIndex.at(b)) != t2->bSiteIndex.at(bSiteIndexOfBond.at(b)))
			{
				clear();
				return false;
			}
			
			//tell t2 we have visited this bond before, now that we are sure it is correct
			t2->setHasVisited(bSiteIndexOfBond.at(b));
			
			//remember that we visited through this bond in this template, in case we get to
			//this molecule again
			this->setHasVisited(b); 
			
			bool goodMatch = t2->compare(m2, mappingSet);
			if(!goodMatch)
			{
				clear();
				return false;
			}
		}		
	}
	clear();
	
	
	//If we got here, we are going to return true, so create the new mappings we will need
	for(tMapIter = tMappings.begin(); tMapIter != tMappings.end(); tMapIter++ )
		mappingSet->add( (*tMapIter)->createNewMapping(m) );
	
	return true;
}

*/









bool TemplateMolecule::contains(TemplateMolecule *tempMol)
{
	bool found = false;
	
	//Create the queues and lists
	queue <TemplateMolecule *> q;
	list <TemplateMolecule *> t;
	queue <int> d;
	int currentDepth = 0;
	
	//First add this molecule
	q.push(this);
	t.push_back(this);
	d.push(currentDepth+1);
	this->hasVisited=true;
	
	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		TemplateMolecule *cM = q.front();
		currentDepth = d.front();
		q.pop(); 
		d.pop();
		
		//Check if we found the molecule
		if(cM==tempMol) {
			found=true;
			break;
		}
		
		//Loop through the bonds to get the next set of Templates to search
		int bMax = cM->getNumBindingSites();
		for(int b=0; b<bMax; b++)
		{
			if(cM->isBindingSiteBonded(b))
			{
				TemplateMolecule *neighbor = cM->getBondedTemplateMolecule(b);
				if(!neighbor->hasVisited)
				{
					neighbor->hasVisited=true;
					t.push_back(neighbor);
					q.push(neighbor);
					d.push(currentDepth+1);
				}
				

			}
		}
	}
	
	
	//clear the has visitedMolecule values
	list <TemplateMolecule *>::iterator tmIter;
	for( tmIter = t.begin(); tmIter != t.end(); tmIter++ )
  		(*tmIter)->hasVisited=false;
	
	return found;
}













void TemplateMolecule::printDetails() const
{
	cout<<"~~ Template Molecule of type: " << parentMoleculeType->getName() << endl;
	cout<<"      constraints on states: ";
	if(stateIndex.size()==0) cout<<"none.";
	for(unsigned int s=0; s<stateIndex.size(); s++)
	{
		cout<< parentMoleculeType->getComponentName(stateIndex.at(s)) <<"=";
		cout<< stateValue.at(s) << " ";
	}
	cout<<endl<<"      constraints on bonds: ";
	if(bSiteIndex.size()==0) cout<<"none.";
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		cout<< parentMoleculeType->getComponentName(bSiteIndex.at(b)) <<"=";
		if(bonds.at(b)==0) cout<<"empty ";
		else cout<<"full "; 
	}
	cout<<endl<<"      has match molecule: ";
	if(matchMolecule==0) cout<< "no";
	else cout<<"yes";
	
	cout<<endl<<"      have traversed bond: ";
	if(bSiteIndex.size()==0) cout<<"n/a";
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		cout<< parentMoleculeType->getComponentName(bSiteIndex.at(b)) <<"=";
		if(hasVisitedBond.at(b)==0) cout<<"no ";
		else cout<<"yes "; 
	}
	
	
	
	cout<<endl<<endl;
}



unsigned int TemplateMolecule::getTemplateBsiteIndexFromMoleculeBsiteIndex(int molBsiteIndex)
{
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
		if(bSiteIndex.at(b)==molBsiteIndex)
			return b;
	cerr<<"Could not find the binding site index in a TemplateMolecule: so I could"<<endl;
	cerr<<"Not get you your Template Binding site index.  I am thus quitting."<<endl;
	exit(1);
	return 0;
}


















void TemplateMolecule::addMapGenerator(MapGenerator *mg)
{
	mapGenerators.push_back(mg);
}


bool TemplateMolecule::compare(Molecule * m, MappingSet *mappingSet)
{	
	//First check if we've been here before
	if(matchMolecule!=0) {
		if(matchMolecule == m) { return true;	}
		else { clear(); return false; }
	}
	
	//Next, check if we have the same type
	if(m->getMoleculeType()!=this->parentMoleculeType) {
		clear();
		return false;
	}
	//Check that states are correct
	int s;
	for(s=0, intVecIter=stateIndex.begin(); intVecIter!=stateIndex.end(); intVecIter++, s++) {
		if(m->getComponentState((*intVecIter)) != stateValue.at(s)) { 
			clear();
			return false; 
		}
	}
	
	for(s=0, intVecIter=notStateIndex.begin(); intVecIter!=notStateIndex.end(); intVecIter++, s++) {
		if(m->getComponentState((*intVecIter)) == notStateValue.at(s)) { 
			clear();
			return false; 
		}
	}
	
	
	//Check if the sites that must be occupied by anything are occupied
	for(intVecIter=sitesThatMustBeOccupied.begin(); intVecIter!=sitesThatMustBeOccupied.end(); intVecIter++)
	{
		//cout<<"Checking site that must be occupied: "<< sitesThatMustBeOccupied.at(b) <<endl;
		//If the site is not occupied in the molecule, no match
		if(!m->isBindingSiteBonded((*intVecIter))) {
			clear();
			return false;
		}
	}
	
	
	
	
	//We have matched so far, so set our match
	matchMolecule = m;
	
	//Cycle through the bonds
	for(unsigned int b=0; b<bSiteIndex.size(); b++)
	{
		if(hasVisitedBond.at(b)==true) continue;
		if(bonds.at(b)==0) //template binding site is open, so molecule binding site must be too.
		{
			if( m->isBindingSiteBonded(bSiteIndex.at(b))) {
				clear();
				return false;
			}
		}
		else //template binding site has a bond
		{
			if(m->isBindingSiteOpen(bSiteIndex.at(b))) {
				clear();
				return false;
			}
			
			//get template that is bound to this binding site
			TemplateMolecule * t2 = this->getBondedTemplateMolecule(b);
			Molecule * m2 = m->getBondedMolecule(bSiteIndex.at(b));
			
			//We have to remember to check that it was the right bond site on the other molecule
			//too.  This means that the bond site in m2 that m is connected to must be the same
			//bond site in t2 that t1 is connected to.  If not, then we must say false
			if(m->getComponentIndexOfBond(bSiteIndex.at(b)) != t2->bSiteIndex.at(bSiteIndexOfBond.at(b))) {
				clear();
				return false;
			}
			
			//tell t2 we have visited this bond before, now that we are sure it is correct
			t2->setHasVisited(bSiteIndexOfBond.at(b));
			
			//remember that we visited through this bond in this template, in case we get to
			//this molecule again
			this->setHasVisited(b); 
			
			bool goodMatch = t2->compare(m2, mappingSet);
			if(!goodMatch)
			{
				clear();
				return false;
			}
		}		
	}
	clear();
	

	//If we got here, we are going to return true, so make the mappings we need
	for(mgIter = mapGenerators.begin(); mgIter != mapGenerators.end(); mgIter++ )
		(*mgIter)->map(mappingSet,m);
	
	return true;
}





