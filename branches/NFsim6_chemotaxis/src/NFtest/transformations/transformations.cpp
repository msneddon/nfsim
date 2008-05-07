



#include "transformations.hh"

using namespace NFtest_transformations;
using namespace NFcore;




void NFtest_transformations::run()
{
	cout<<"Testing Transformation Related Classes and Functions..."<<endl;
	
	System *s = new System("transform test");
	MoleculeType *X = tc_createX(s, 1000);
	MoleculeType *Y = tc_createY(s, 1000);
	
	
	
	//Create 1 template molecule
//	TemplateMolecule *xTemp1 = new TemplateMolecule(X);
//	xTemp1->addStateValue("m",0);
//	
//	//Create the reaction
//	vector <TemplateMolecule *> templates;
//	templates.push_back( xTemp1 );
//	ReactionClass *r1 = new ReactionClass("R1",templates, 2.0);
//	s->addReaction(r1);
//	
//	//Give the transformations for this reaction
//	Transformation *t1 = Transformation::genStateChangeTransform(xTemp1,"m",5, r1);
//	t1->printDetails();
//	
//	r1->printFullDetails();
//	for(int i=0; i<20; i++)
//		r1->tryToAdd(X->getMolecule(i),0);
//	r1->printFullDetails();
//	
//	vector <MappingSet *> mappingSets;
//	//for(int i=0; i<5; i++)
//		r1->pickMappingSets(0,mappingSets);
//	
//	MappingSet::transform(mappingSets);
	
	
	
	
	
	//Test binding Reaction ....
	TemplateMolecule *xTemp2 = new TemplateMolecule(X);
	xTemp2->addEmptyBindingSite("a2");
	TemplateMolecule *yTemp2 = new TemplateMolecule(Y);
	yTemp2->addEmptyBindingSite("b2");
	
	vector <TemplateMolecule *> templates2;
	templates2.push_back( xTemp2 );
	templates2.push_back( yTemp2 );
	ReactionClass *r2 = new ReactionClass("R2",templates2, 2.0);
	s->addReaction(r2);
	Transformation *t2 = Transformation::genBindingTransform(xTemp2,yTemp2,"a2","b2",r2);
	
	for(int i=50; i<100; i++)
	{
		r2->tryToAdd(X->getMolecule(i),0);
		r2->tryToAdd(Y->getMolecule(i),1);
	}
	
	
	
	//vector <MappingSet *> mappingSets2;
		//for(int i=0; i<5; i++)
	//r2->pickMappingSets(0,mappingSets2);
		
	//MappingSet::transform(mappingSets2);
	
	
	
	r2->fire2(0);
	
	
	
	
	
	
	
	
	/* test ReactantList....
	int molecule = 0;
	MappingSet *mappingSet = new MappingSet();
	xTemp1->compare(X->getMolecule(molecule),mappingSet);
	vector <MappingSet *> mSets;
	mSets.push_back(mappingSet);
	MappingSet::transform(mSets);
	
	ReactantList * rl = new ReactantList(3);
	rl->push(X->getMolecule(molecule), mappingSet);
	
	///////////////////////////////////
	for(int k=1; k<10; k++)
	{
	mappingSet = new MappingSet();
	xTemp1->compare(X->getMolecule(k),mappingSet);
	mSets.clear();
	mSets.push_back(mappingSet);
	MappingSet::transform(mSets);
		
	rl->push(X->getMolecule(k), mappingSet);
	}
	
	
	rl->printDetails();
	rl->pop(X->getMolecule(5));
	
	rl->printDetails();
	delete rl;
	*/
	
	
	
	
	
	
	
	//Transformation *t2 = Transformation::genBindingTransform(xTemp1,xTemp2,"a0","a1", r1);
	//t2->printDetails();
	
	
	//Transformation *t3 = Transformation::genUnbindingTransform(xTemp3,"a0");
		
	//t3->printDetails();
	
	
	
	

	
	//MappingSet *mappingSet = new MappingSet();
	
//	TemplateMapping * tempMapping = new TemplateMapping(xTemp,Mapping::STATE,"m");
//	mappingSet.add(tempMapping->createNewMapping(X->getMolecule(0)));
	
//	Mapping * mapping = mappingSet.getMapping(0);
	
	//mappingSet->add();
	
	//xTemp->printDetails();
	
	
	//xTemp1->compare(X->getMolecule(0),mappingSet);
	//mappingSet->printDetails();
	
	//mappingSet.clear();

	//xTemp2->compare(X->getMolecule(3),mappingSet);
	//mappingSet->printDetails();
	//mappingSet.clear();
	
	
	
	//cout<<"-----------------"<<endl<<endl;
//	vector <MappingSet *> mSets;
	//mSets.push_back(mappingSet);
	
	
	
	//MappingSet::transform(mSets);
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	//xTemp->compare(
	//xTemp->compare(X->getMolecule(0),mappingSet);
	
//	
//	
//	
//	TemplateMapping *tempMapping = new TemplateMapping(xTemp, Mapping::STATE,"m");
//	
//	tempMapping->printDetails();
//	
//	
//	
//	Transformation * t1 = new StateChangeTransform(tempMapping, 5);
//	Transformation * t2 = new BindingTransform(tempMapping, tempMapping);
//	
//	tempMapping->printDetails();
//	
//	Mapping *m = tempMapping->createNewMapping(X->getMolecule(0));
//	
//	
//	
//	
//	Mapping ** mappings = new Mapping *[1];
//	mappings[0] = m;
//	
//	
//	X->getMolecule(0)->printDetails();
//	
//	m->getTransformation()->transform(mappings);
//	
//	X->getMolecule(0)->printDetails();
	
	
	
	
	delete s;
	
}



MoleculeType * NFtest_transformations::tc_createX(System * s, int count)
{
	int numOfBsites = 4;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "a0";
	bSiteNames[1] = "a1";
	bSiteNames[2] = "a2";
	bSiteNames[3] = "a3";
	char  numOfStates = 2;
	const char ** stateNames = new const char * [numOfStates];
	stateNames[0] = "p";
	stateNames[1] = "m";
	int * stateValues = new int [numOfStates];
	stateValues[0] = 0;
	stateValues[1] = 0;
	MoleculeType *L = new MoleculeType("X",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}

MoleculeType * NFtest_transformations::tc_createY(System * s, int count)
{
	int numOfBsites = 4;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "b0";
	bSiteNames[1] = "b1";
	bSiteNames[2] = "b2";
	bSiteNames[3] = "b3";
	char  numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *L = new MoleculeType("Y",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}