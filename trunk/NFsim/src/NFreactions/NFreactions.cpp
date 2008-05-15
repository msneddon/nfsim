

#include "NFreactions.hh"

using namespace NFcore;
using namespace std;


void NFcore::test()
{
	cout<<"Testing rxns..."<<endl;
	
	
	System *s = new System("boo");
	
	
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "y";
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = "p";
	int * stateValues = new int [numOfStates];
	stateValues[0] = 1;
	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	
	numOfBsites = 1;
	bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "x";
	numOfStates = 1;
	stateNames = new string [numOfStates];
	stateNames[0] = "m";
	stateValues = new int [numOfStates];
	stateValues[0] = 1;
	MoleculeType *molY = new MoleculeType("MolY",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	
	molX->populateWithDefaultMolecules(500);
	molY->populateWithDefaultMolecules(500);
	
	
	molX->getMolecule(0)->printDetails();
	molY->getMolecule(0)->printDetails();
	
	
	
	TemplateMolecule *xTemp1 = new TemplateMolecule(molX);
	xTemp1->addStateValue("p",1);
	TemplateMolecule *yTemp2 = new TemplateMolecule(molY);
	yTemp2->addStateValue("m",1);
	//TemplateMolecule::bind(xTemp1,"y",yTemp2,"x");
	//xTemp2->printDetails();
	
	TemplateMolecule *xTemp3 = new TemplateMolecule(molX);
	xTemp3->addStateValue("p",1);
	TemplateMolecule *xTemp4 = new TemplateMolecule(molX);
	xTemp4->addStateValue("p",1);
	
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp1 );
	templates.push_back( yTemp2 );
	//templates.push_back( xTemp3 );
	//templates.push_back( xTemp4 );
	
	
	
	
	
	
	TransformationSet *t = new TransformationSet(templates);
	//t->addStateChangeTransform(xTemp1,"p",0);
	//t->addStateChangeTransform(yTemp2,"m",3);
	t->addBindingTransform(yTemp2,"x",xTemp1,"y");
	cout<<endl;
	
	
	
	
	
	t->finalize();
	
	ReactantList * rl_0 = new ReactantList(0, t, 15);
	ReactantList * rl_1 = new ReactantList(1, t, 15);
	//rl->printDetails();
	
	
	
	
	for(int i=0; i<12; i++)
	{
		MappingSet *ms = rl_0->pushNextAvailableMappingSet();
		if(!xTemp1->compare(molX->getMolecule(20+i),ms)) {
			rl_0->popLastMappingSet();
		}
	}
	for(int i=0; i<12; i++)
	{
		MappingSet *ms = rl_1->pushNextAvailableMappingSet();
		if(!yTemp2->compare(molY->getMolecule(20+i),ms)) {
			rl_1->popLastMappingSet();
		}
	}
	
	
	rl_0->printDetails();
	rl_1->printDetails();
	
	//cout<<"ms->getId() = "<<ms->getId()<<endl;
	//rl->removeMappingSet(5);
	

	//rl->printDetails();
	
	
	
	MappingSet **ms = new MappingSet *[2];
	rl_0->pickRandom(ms[0]);
	rl_1->pickRandom(ms[1]);
	cout<<"Retrieving randomly: " <<ms[0]->get(0)->getMolecule()->getUniqueID()<<endl;
	cout<<"Retrieving randomly: " <<ms[1]->get(0)->getMolecule()->getUniqueID()<<endl;
	
	ms[0]->get(0)->getMolecule()->printDetails();
	ms[1]->get(0)->getMolecule()->printDetails();
	t->transform(ms);
	cout<<endl<<"--------"<<endl;
	ms[0]->get(0)->getMolecule()->printDetails();
	ms[1]->get(0)->getMolecule()->printDetails();
	
	cout<<endl<<endl<<endl;
	delete [] ms;
	delete t; delete rl_0; 
	delete rl_1;
	
	delete s;
}








void NFcore::test_simple()
{
	cout<<"Testing rxns..."<<endl;
	
	
	System *s = new System("boo");
	
	
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "y";
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = "p";
	int * stateValues = new int [numOfStates];
	stateValues[0] = 1;
	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	molX->populateWithDefaultMolecules(500);
	molX->getMolecule(0)->printDetails();
	
	
	
	TemplateMolecule *xTemp1 = new TemplateMolecule(molX);
	xTemp1->addStateValue("p",1);
	TemplateMolecule *xTemp2 = new TemplateMolecule(molX);
	xTemp2->addStateValue("p",1);
	//TemplateMolecule::bind(xTemp1,"y",yTemp2,"x");
	//xTemp2->printDetails();
	
	
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp1 );
	
	
	
	
	
	
	TransformationSet *t = new TransformationSet(templates);
	t->addStateChangeTransform(xTemp1,"p",0);
	//t->addBindingTransform(xTemp1,"y",yTemp2,"x");
	cout<<endl;
	
	
	
	
	
	t->finalize();
	
	
	ReactantList * rl = new ReactantList(0, t, 15);
	//rl->printDetails();
	
	
	
	
	for(int i=0; i<12; i++)
	{
		MappingSet *ms = rl->pushNextAvailableMappingSet();
		if(!xTemp1->compare(molX->getMolecule(20+i),ms)) {
			rl->popLastMappingSet();
		}
	}
	
	
	rl->printDetails();
	//cout<<"ms->getId() = "<<ms->getId()<<endl;
	//rl->removeMappingSet(5);
	

	//rl->printDetails();
	
	
	
	MappingSet **ms = new MappingSet *[2];
	rl->pickRandom(ms[0]);
	cout<<"Retrieving randomly: " <<ms[0]->get(0)->getMolecule()->getUniqueID()<<endl;
	
	ms[0]->get(0)->getMolecule()->printDetails();
	t->transform(ms);
	cout<<endl<<"--------"<<endl;
	ms[0]->get(0)->getMolecule()->printDetails();
	
	cout<<endl<<endl<<endl;
	delete t; delete rl; 
	delete s;
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
	finalized = false;
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<n_reactants; r++)  {
		cout<<"Count in ["<<r<<"]: "<<this->transformations[r].size()<<endl;
	}
	
	for(unsigned int r=0; r<n_reactants; r++)  {
		Transformation *t;
		while(transformations[r].size()>0)
		{
			t = transformations[r].back();
			transformations[r].pop_back();
			delete t;
		}
	}
	delete [] transformations;
	delete [] reactants;
	this->n_reactants = 0;
}

bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string stateName, int finalStateValue)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		exit(1);
	}
	
	// 2) Create a Transformation object to remember the information
	unsigned int stateIndex = t->getMoleculeType()->getStateIndex(stateName);
	Transformation *transformation = Transformation::genStateChangeTransform(stateIndex, finalStateValue);
	
	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);
	
	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}
bool TransformationSet::addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex2==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!"<<endl;
		exit(1);
	}
	
	//Find the index of the respective binding sites
	unsigned int bSiteIndex1 = t1->getMoleculeType()->getBindingSiteIndex(bSiteName1);
	unsigned int bSiteIndex2 = t2->getMoleculeType()->getBindingSiteIndex(bSiteName2);
	
	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	if(reactantIndex1==reactantIndex2)
		transformation1 = Transformation::genBindingTransform1(bSiteIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
	else
		transformation1 = Transformation::genBindingTransform1(bSiteIndex1, reactantIndex2, transformations[reactantIndex2].size());

	Transformation *transformation2 = Transformation::genBindingTransform2(bSiteIndex2);
	
	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);
	
	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t2->addMapGenerator(mg2);
	
	return true;
}
bool TransformationSet::addUnbindingTransform(TemplateMolecule *t, string bSiteName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		exit(1);
	}
	
	// 2) Create a Transformation object to remember the information
	unsigned int bSiteIndex = t->getMoleculeType()->getBindingSiteIndex(bSiteName);
	Transformation *transformation = Transformation::genUnbindingTransform(bSiteIndex);
	
	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);
	
	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	
	return true;
}






int TransformationSet::find(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot search for a templateMolecule once it has been finalized!"<<endl; exit(1); }
	int findIndex = -1;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(this->reactants[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r;
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
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }
	
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
				Mapping *m2 = mappingSets[transformations[r].at(t)->getPartnerReactantIndex()]->get(transformations[r].at(t)->getPartnerMappingIndex());
				Molecule::bind(m1->getMolecule(),m1->getIndex(), m2->getMolecule(), m2->getIndex());
			}
		}
	}
	return true;
}
bool TransformationSet::getListOfProducts(MappingSet **mappingSets, list<Molecule *> &products, int traversalLimit)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }
	
	for(unsigned int r=0; r<n_reactants; r++)  {
		MappingSet *ms = mappingSets[r];
		ms->get(0)->getMolecule()->traverseBondedNeighborhood(products, traversalLimit);
		
	}
	return true;
}


MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(!finalized) { cerr<<"TransformationSet cannot generate blank mapping if it is not finalized!"<<endl; exit(1); }
	if(reactantIndex>=n_reactants) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(mappingSetId, transformations[reactantIndex]);
}

void TransformationSet::finalize()
{
	//Be sure to add at least a blank transformation to every reactant if there is no transformation
	//specified so that we count the reactants even if we don't do anything to it.
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(transformations[r].size()==0) {
			transformations[r].push_back(Transformation::genEmptyTransform());
			MapGenerator *mg = new MapGenerator(transformations[r].size()-1);
			getTemplateMolecule(r)->addMapGenerator(mg);
		}
	}
	finalized = true;
}






NFcore::Transformation::Transformation()
{
	type=Transformation::SKIP;
	newStateValue = -1;
	stateORsiteIndex=0;
	otherReactantIndex=0;
	otherMappingIndex=0;
}
NFcore::Transformation::~Transformation()
{
	type=Transformation::SKIP;
	newStateValue = -1;
	stateORsiteIndex=0;
	otherReactantIndex=0;
	otherMappingIndex=0;
}
NFcore::Transformation * Transformation::genEmptyTransform()
{
	return new Transformation();
}
NFcore::Transformation * Transformation::genStateChangeTransform(unsigned int stateIndex, int newStateValue)
{
	Transformation *t = new Transformation();
	t->type = Transformation::STATE_CHANGE;
	t->stateORsiteIndex=stateIndex;
	t->newStateValue = newStateValue;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex)
{	
	Transformation *t = new Transformation();
	t->type = Transformation::BINDING;
	t->stateORsiteIndex=bSiteIndex;
	t->otherReactantIndex = otherReactantIndex;
	t->otherMappingIndex = otherMappingIndex;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genBindingTransform2(unsigned int bSiteIndex)
{
	Transformation *t = new Transformation();
	t->type = Transformation::SKIP;
	t->stateORsiteIndex=bSiteIndex;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genUnbindingTransform(unsigned int bSiteIndex)
{
	Transformation *t = new Transformation();
	t->type = Transformation::UNBINDING;
	t->stateORsiteIndex=bSiteIndex;
	return t;
}
NFcore::Transformation * NFcore::Transformation::genAddMoleculeTransform()
{
	cerr<<"Warning genAddMoleculeTransform() not yet working!!"<<endl;
	return 0;
}
NFcore::Transformation * NFcore::Transformation::genRemoveMoleculeTransform()
{
	cerr<<"Warning genRemoveMoleculeTransform() not yet working!!"<<endl;
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
	return true;
}








NFcore::MappingSet::MappingSet(unsigned int id, vector <Transformation *> &transformations)
{
	this->id = id;
	this->isSet = false;
	this->n_mappings = transformations.size();
	this->mappings = new Mapping *[n_mappings];
	
	for(unsigned int t=0; t<n_mappings; t++) {
		mappings[t] = new Mapping(transformations.at(t)->getType(), transformations.at(t)->getStateOrSiteIndex() );
	}
}
NFcore::MappingSet::~MappingSet()
{
	for(unsigned int t=0; t<n_mappings; t++) {
		delete mappings[t];
	}
	delete [] mappings;
}
			
bool NFcore::MappingSet::set(unsigned int mappingIndex, Molecule *m)
{
	if(mappingIndex>=n_mappings) {
		cerr<<"Out of bounds access to a mapping in set function of mapping set!"<< mappingIndex<<" but max is " << n_mappings<<endl;
		return false;
	}
	mappings[mappingIndex]->setMolecule(m);
	isSet = true;
	return true;
}
Mapping *NFcore::MappingSet::get(unsigned int mappingIndex)
{
	if(mappingIndex>=n_mappings) {
		cerr<<"Out of bounds access to a mapping in get function of a mapping set!"<<endl;
		return false;
	}
	return mappings[mappingIndex];
}
bool NFcore::MappingSet::clear()
{
	isSet = false;
	
	//This is not entirely necessary, but for now can be used for debugging
	for(unsigned int t=0; t<n_mappings; t++) {
		mappings[t]->clear();
	}
	
	return true;
}
	













NFcore::Mapping::Mapping(unsigned int type, unsigned int index)
{
	this->type = type;
	this->index = index;
}
NFcore::Mapping::~Mapping()
{
	type=Transformation::SKIP;
	index=0;
	clear();
}
		
unsigned int NFcore::Mapping::getType() const
{
	return this->type;
}
unsigned int NFcore::Mapping::getIndex() const
{
	return this->index;
}
Molecule * NFcore::Mapping::getMolecule() const
{
	if(m==NULL) {
		cout<<"Trying to get a molecule from a null mapping!!"<<endl;
		return 0;
	}
	return m;
}
		
void NFcore::Mapping::clear()
{
	this->m=NULL;
}
bool NFcore::Mapping::setMolecule(Molecule *m)
{
	this->m = m;
	return true;
}
		
		
		










		
		
