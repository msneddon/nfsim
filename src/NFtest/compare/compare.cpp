



#include "compare.hh"

using namespace NFtest_compare;
using namespace NFcore;

bool testCompare01();
bool testCompare02();
bool testCompare03();
bool testCompare04();







void NFtest_compare::run()
{
	cout<<"Testing Template Compare Function..."<<endl;
	
	int testNum = 1;
	
	if(testCompare01())
		cout<<"Test "<< testNum++ <<": " << "Pass"<<endl;
	else
		cout<<"Test "<< testNum++ <<": " << "FAIL!!!!"<<endl;
	if(testCompare03())
		cout<<"Test "<< testNum++ <<": " << "Pass"<<endl;
	else
		cout<<"Test "<< testNum++ <<": " << "FAIL!!!!"<<endl;
	
}

//Test that see if compare accurately checks "Type" of molecule
bool testCompare01()
{
	System *s=new System("test system only");
	MoleculeType * X = tc_createX(s,1);
	MoleculeType * Y = tc_createY(s,1);
	
	int n_templates = 2;
	TemplateMolecule ** templates = new TemplateMolecule *[n_templates];
	templates[0] = new TemplateMolecule(X);
	templates[1] = new TemplateMolecule(Y);
	
	if(templates[0]->compare(X->getMolecule(0))!=true) { delete s; return false; }
	if(templates[0]->compare(Y->getMolecule(0))!=false) { delete s; return false; }
	if(templates[1]->compare(X->getMolecule(0))!=false) { delete s; return false; }
	if(templates[1]->compare(Y->getMolecule(0))!=true) { delete s; return false; }
	
	delete s;
	return true;	
}


//Test that see if simple binding / not binding works
bool testCompare02()
{
	System *s=new System("test system only");
	MoleculeType * X = tc_createX(s,5);
	MoleculeType * Y = tc_createY(s,5);
	
	int n_templates = 6;
	TemplateMolecule ** templates = new TemplateMolecule *[n_templates];
	templates[0] = new TemplateMolecule(X);
	TemplateMolecule * boundY = new TemplateMolecule(Y);
	TemplateMolecule::bind(templates[0], "a2", boundY, "b2");
	
	templates[1] = new TemplateMolecule(Y);
	TemplateMolecule * boundX = new TemplateMolecule(X);
	TemplateMolecule::bind(templates[0], "a2", boundX, "b2");
	
	templates[2] = new TemplateMolecule(X);
	TemplateMolecule * boundY2 = new TemplateMolecule(X);
	TemplateMolecule::bind(templates[0], "a2", boundY, "b2");
	
	templates[3] = new TemplateMolecule(Y);
	TemplateMolecule * boundX2 = new TemplateMolecule(X);
	TemplateMolecule::bind(templates[0], "a2", boundY, "b2");
	
	
	
	
	
	if(templates[0]->compare(X->getMolecule(0))!=true) { delete s; return false; }
	if(templates[0]->compare(Y->getMolecule(0))!=false) { delete s; return false; }
	if(templates[1]->compare(X->getMolecule(0))!=false) { delete s; return false; }
	if(templates[1]->compare(Y->getMolecule(0))!=true) { delete s; return false; }
	
	delete s;
	return true;	
}

// Test the occupied binding site checks
bool testCompare03()
{
	System *s=new System("test system only");
	MoleculeType * X = tc_createX(s,5);
	MoleculeType * Y = tc_createY(s,5);
	
	int n_templates = 1;
	TemplateMolecule ** templates = new TemplateMolecule *[n_templates];
	templates[0] = new TemplateMolecule(X);
	templates[0]->addOccupiedBindingSite("a0");
	
	Molecule::bind(X->getMolecule(0),"a1", Y->getMolecule(0), "b2");
	if(templates[0]->compare(X->getMolecule(0))!=false) { delete s; return false; }
	Molecule::bind(X->getMolecule(0),"a0", Y->getMolecule(1), "b2");
	if(templates[0]->compare(X->getMolecule(0))!=true) { delete s; return false; }
	templates[0]->addOccupiedBindingSite("a3");
	if(templates[0]->compare(X->getMolecule(0))!=false) { delete s; return false; }
	Molecule::bind(X->getMolecule(0),"a3", Y->getMolecule(2), "b2");
	if(templates[0]->compare(X->getMolecule(0))!=true) { delete s; return false; }
	Molecule::unbind(X->getMolecule(0),"a1");
	if(templates[0]->compare(X->getMolecule(0))!=true) { delete s; return false; }
	Molecule::unbind(X->getMolecule(0),"a0");
	if(templates[0]->compare(X->getMolecule(0))!=false) { delete s; return false; }
	
	delete s;
	return true;	
}



MoleculeType * NFtest_compare::tc_createX(System * s, int count)
{
	int numOfBsites = 4;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "a0";
	bSiteNames[1] = "a1";
	bSiteNames[2] = "a2";
	bSiteNames[3] = "a3";
	int numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *L = new MoleculeType("X",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}

MoleculeType * NFtest_compare::tc_createY(System * s, int count)
{
	int numOfBsites = 4;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "b0";
	bSiteNames[1] = "b1";
	bSiteNames[2] = "b2";
	bSiteNames[3] = "b3";
	int numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *L = new MoleculeType("Y",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}

