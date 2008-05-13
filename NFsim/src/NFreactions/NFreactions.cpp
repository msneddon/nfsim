

#include "NFreactions.hh"

using namespace NFcore;
using namespace std;


void NFcore::test()
{
	cout<<"Testing rxns..."<<endl;
	
	
	System *s = new System("boo");
	// create MoleculeType X  with one binding site and one state 
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "y";
	
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = "p";
	
	//This is the default state value that new molecules are created with
	int * stateValues = new int [numOfStates];
	stateValues[0] = 0;
	
	//When we create a molecule, it automatically adds itself to the system, 
	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	
	
	TemplateMolecule *xTemp1 = new TemplateMolecule(molX);
	xTemp1->addStateValue("p",1);
	TemplateMolecule *xTemp2 = new TemplateMolecule(molX);
	xTemp2->addStateValue("p",1);
	TemplateMolecule::bind(xTemp1,"y",xTemp2,"y");
	xTemp2->printDetails();
	
	TemplateMolecule *xTemp3 = new TemplateMolecule(molX);
	xTemp3->addStateValue("p",1);
	TemplateMolecule *xTemp4 = new TemplateMolecule(molX);
	xTemp4->addStateValue("p",1);
	
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp1 );
	templates.push_back( xTemp2 );
	templates.push_back( xTemp3 );
	templates.push_back( xTemp4 );
	
	
	

	
	
	TransformationSet *t = new TransformationSet(templates);
	t->addStateChangeTransform(xTemp4,"p",1);
	t->addStateChangeTransform(xTemp4,"p",1);
	delete t;
	
	
}

















TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates)
{
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);
	
	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [n_reactants];
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<n_reactants; r++)  {
		cout<<"Count in ["<<r<<"]: "<<this->transformations[r].size()<<endl;
	}
	
	
	this->n_reactants = 0;
	
	
	for(unsigned int r=0; r<n_reactants; r++)  {
		this->transformations = new vector <Transformation *>;
	}
}

bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string stateName, int finalStateValue)
{
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		exit(1);
	}
	
	// 2) Create a Transformation object to remember the information
	Transformation *transformation = Transformation::genStateChangeTransform(finalStateValue);
	
	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);
	
	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size());
	t->addMapGenerator(mg);
	return true;
}
int TransformationSet::find(TemplateMolecule *t)
{
	int findIndex = -1;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(this->reactants[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r;
				cout<<"found at index: "<<r<<endl;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	return findIndex;
}
bool TransformationSet::transform(MappingSet **mappingSets)
{
	for(unsigned int r=0; r<n_reactants; r++)  {
		MappingSet *ms = mappingSets[r];
		for(unsigned int t=0; t<transformations[r].size(); t++)
		{
			unsigned int type = transformations[r].at(t)->getType();
			if(type == Transformation::SKIP) {
				continue;
			}
			else if(type == Transformation::STATE_CHANGE) {
				Mapping *m = ms->get(t);
				m->getMolecule()->setState(m->getIndex(),transformations[r].at(t)->getNewStateValue());
			}
			else if(type == Transformation::UNBINDING) {
				Mapping *m = ms->get(t);
				Molecule::unbind(m->getMolecule(),m->getIndex());
			}
			else if(type == Transformation::BINDING) {
				Mapping *m1 = ms->get(t);
				Mapping *m2 = mappingSets[transformations[r].at(t)->getBiPartnerReactantIndex()]->get(transformations[r].at(t)->getBiPartnerMappingIndex());
				Molecule::bind(m1->getMolecule(),m1->getIndex(), m2->getMolecule(), m2->getIndex());
			}
		}
	}
}
MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(reactantIndex>=n_reactants) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(transformations[reactantIndex]);
}








NFcore::Transformation::Transformation()
{
	type=Transformation::SKIP;
	newStateValue = -1;
	reactant=0;
	mappingIndex=0;
}
NFcore::Transformation * NFcore::Transformation::genStateChangeTransform(int newStateValue)
{
	Transformation *t = new Transformation();
	t->type = Transformation::STATE_CHANGE;
	t->newStateValue = newStateValue;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genBindingTransform1(unsigned int reactant, unsigned int mappingIndex)
{	
	Transformation *t = new Transformation();
	t->type = Transformation::BINDING;
	t->reactant = reactant;
	t->mappingIndex = mappingIndex;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genBindingTransform2()
{
	return new Transformation();
}
NFcore::Transformation * NFcore::Transformation::genUnbindingTransform()
{
	return 0;
}
NFcore::Transformation * NFcore::Transformation::genAddMoleculeTransform()
{
	return 0;
}
NFcore::Transformation * NFcore::Transformation::genRemoveMoleculeTransform()
{
	return 0;
}









NFcore::MapGenerator::MapGenerator(unsigned int mappingIndex)
{
	this->mappingIndex = mappingIndex;
}
NFcore::MapGenerator::~MapGenerator()
{
}
bool NFcore::MapGenerator::map(MappingSet *mappingSet, Molecule *molecule)
{
	mappingSet->set(mappingIndex,molecule);
}








NFcore::MappingSet::MappingSet(vector <Transformation *> &transformations)
{
	
}
NFcore::MappingSet::~MappingSet()
{
	
}
			
bool NFcore::MappingSet::set(unsigned int mappingIndex, Molecule *m)
{
	
}
Mapping *NFcore::MappingSet::get(unsigned int mappingIndex)
{
	
}
bool NFcore::MappingSet::clear()
{
	
}
	













NFcore::Mapping::Mapping(unsigned int type, unsigned int index)
{
	
	
	
}
NFcore::Mapping::~Mapping()
{
	
}
		
unsigned int NFcore::Mapping::getType() const
{
	
	return 0;
}
unsigned int NFcore::Mapping::getIndex() const
{
	
	return 0;
}
Molecule * NFcore::Mapping::getMolecule() const
{
	return 0;
}
		
void NFcore::Mapping::clear()
{
	
}
bool NFcore::Mapping::setMolecule(Molecule *m)
{
	return false;
}
		
		
		
















//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


ReactantList::ReactantList(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity=50)
{
	this->n_mappingSets = 0;
	this->capacity = init_capacity;
	this->reactantIndex = reactantIndex;
	this->ts=ts;
	this->mappingSets = new MappingSet *[init_capacity];
	this->msPositionMap = new unsigned int [init_capacity];
	for(unsigned int i=0; i<this->capacity; i++)
	{
		mappingSets[i]=ts->generateBlankMappingSet(reactantIndex,i);
		msPositionMap[i]=i;
	}
}

ReactantList::~ReactantList()
{
	for(unsigned int i=0; i<capacity; i++)
	{
		delete mappingSets[i];
		msPositionMap[i]=0;
	}
	delete [] mappingSets;
	delete [] msPositionMap;
	this->n_mappingSets = 0;
	this->capacity = 0;
}

/*
void ReactantList::pickRandom(MappingSet *&ms, unsigned int &moleculeID)
{
	unsigned int rand = NFutil::RANDOM_INT(0,n_reactants);
	ms = reactants[rand];
	moleculeID = reactantMoleculeIDs[rand];
}
*/


unsigned int ReactantList::size() const
{
	return n_mappingSets;
}


MappingSet * ReactantList::pushNextAvailableMappingSet()
{
	//Check if we are going to exceed capacity
	if(n_mappingSets>=capacity)
	{
		//Copy everything over to new arrays that are double the size
		MappingSet ** new_mappingSets = new MappingSet *[capacity*2];
		unsigned int * new_msPositionMap = new unsigned int [capacity*2];
		for(unsigned int i=0; i<capacity; i++)  {
			new_mappingSets[i] = mappingSets[i];
			new_msPositionMap[i] = msPositionMap[i];
		}
		for(unsigned int i=capacity; i<capacity*2; i++)  {
			new_mappingSets[i] = ts->generateBlankMappingSet(reactantIndex, i);
			new_msPositionMap[i] = i;
		}
		
		//Swap the copied data with the real data and double the capacity
		delete [] mappingSets;
		delete [] msPositionMap;
		mappingSets = new_mappingSets;
		msPositionMap = new_msPositionMap;
		capacity*=2;
	}
	
	//Increase the number of reactants, and return the activated mappingSet
	n_mappingSets++;
	return mappingSets[n_mappingSets-1];
}


void ReactantList::popLastMappingSet()
{
	//Clear out the mappingSet (just in case) and decrease the count
	mappingSets[n_mappingSets-1]->clear();
	n_mappingSets--;
}

void ReactantList::removeMappingSet(unsigned int mappingSetId)
{
	//First, get the position of the mappingSet we need to remove
	unsigned int pos = msPositionMap[mappingSetId];
	
	//If the array has only one element, or we just happened to select the last element,
	//then just remove the last element without a swap
	if( pos == (n_mappingSets-1) )
	{
		popLastMappingSet();
		return;
	}
	
	//Otherwise, we have to swap with the last element in the list
	MappingSet *tempMappingSet = mappingSets[pos];
	mappingSets[pos] = mappingSets[n_mappingSets-1];
	mappingSets[n_mappingSets-1] = tempMappingSet;
	msPositionMap[pos] = n_mappingSets-1;
	msPositionMap[n_mappingSets-1] = pos;
	
	
	//Make sure we clear what we don't need
	mappingSets[n_mappingSets-1]->clear();

	//Remember to remove 
	n_mappingSets--;
}


void ReactantList::printDetails()
{
	cout<<"ReactantList that contains: "<<size()<<" MappingSets and has a capacity for "<<capacity<<" total sets."<<endl;
	/*cout<<"Members: [molecule unique Id, private mappingSet index]"<<endl;
	
	
	for(unsigned int i=0; i<capacity; i++)
		cout<<i<<" ";
	cout<<endl;
	for(unsigned i=0; i<capacity; i++)
		if(reactants[i]!=NULL)
			cout<<1<<" ";
		else
			cout<<0<<" ";
	cout<<endl;
	for(unsigned i=0; i<capacity; i++)
			cout<<reactantMoleculeIDs[i]<<" ";
		cout<<endl;
	
	multimap <unsigned int, unsigned int>::iterator it;
	for ( it=moleculeLookupTable.begin(); it != moleculeLookupTable.end(); it++ )
	   cout << "  [" << (*it).first << ", " << (*it).second << "]" << endl;
	
	cout<<endl;*/
}
		
		
