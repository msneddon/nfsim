

#include "NFreactions.hh"

using namespace NFcore;
using namespace std;


void NFcore::test()
{
//	cout<<"Testing rxns..."<<endl;
//
//
//	System *s = new System("boo");
//
//
//	int numOfBsites = 1;
//	string * bSiteNames = new string [numOfBsites];
//	bSiteNames[0] = "y";
//	int numOfStates = 1;
//	string * stateNames = new string [numOfStates];
//	stateNames[0] = "p";
//	int * stateValues = new int [numOfStates];
//	stateValues[0] = 1;
//	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
//
//	numOfBsites = 1;
//	bSiteNames = new string [numOfBsites];
//	bSiteNames[0] = "x";
//	numOfStates = 1;
//	stateNames = new string [numOfStates];
//	stateNames[0] = "m";
//	stateValues = new int [numOfStates];
//	stateValues[0] = 1;
//	MoleculeType *molY = new MoleculeType("MolY",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
//
//	molX->populateWithDefaultMolecules(500);
//	molY->populateWithDefaultMolecules(500);
//
//
//	molX->getMolecule(0)->printDetails();
//	molY->getMolecule(0)->printDetails();
//
//
//
//	TemplateMolecule *xTemp1 = new TemplateMolecule(molX);
//	xTemp1->addStateValue("p",1);
//	TemplateMolecule *yTemp2 = new TemplateMolecule(molY);
//	yTemp2->addStateValue("m",1);
//	//TemplateMolecule::bind(xTemp1,"y",yTemp2,"x");
//	//xTemp2->printDetails();
//
//	TemplateMolecule *xTemp3 = new TemplateMolecule(molX);
//	xTemp3->addStateValue("p",1);
//	TemplateMolecule *xTemp4 = new TemplateMolecule(molX);
//	xTemp4->addStateValue("p",1);
//
//	vector <TemplateMolecule *> templates;
//	templates.push_back( xTemp1 );
//	templates.push_back( yTemp2 );
//	//templates.push_back( xTemp3 );
//	//templates.push_back( xTemp4 );
//
//
//
//
//
//
//	TransformationSet *t = new TransformationSet(templates);
//	//t->addStateChangeTransform(xTemp1,"p",0);
//	//t->addStateChangeTransform(yTemp2,"m",3);
//	t->addBindingTransform(yTemp2,"x",xTemp1,"y");
//	cout<<endl;
//
//
//
//
//
//	t->finalize();
//
//	ReactantList * rl_0 = new ReactantList(0, t, 15);
//	ReactantList * rl_1 = new ReactantList(1, t, 15);
//	//rl->printDetails();
//
//
//
//
//	for(int i=0; i<12; i++)
//	{
//		MappingSet *ms = rl_0->pushNextAvailableMappingSet();
//		if(!xTemp1->compare(molX->getMolecule(20+i),ms)) {
//			rl_0->popLastMappingSet();
//		}
//	}
//	for(int i=0; i<12; i++)
//	{
//		MappingSet *ms = rl_1->pushNextAvailableMappingSet();
//		if(!yTemp2->compare(molY->getMolecule(20+i),ms)) {
//			rl_1->popLastMappingSet();
//		}
//	}
//
//
//	rl_0->printDetails();
//	rl_1->printDetails();
//
//	//cout<<"ms->getId() = "<<ms->getId()<<endl;
//	//rl->removeMappingSet(5);
//
//
//	//rl->printDetails();
//
//
//
//	MappingSet **ms = new MappingSet *[2];
//	rl_0->pickRandom(ms[0]);
//	rl_1->pickRandom(ms[1]);
//	cout<<"Retrieving randomly: " <<ms[0]->get(0)->getMolecule()->getUniqueID()<<endl;
//	cout<<"Retrieving randomly: " <<ms[1]->get(0)->getMolecule()->getUniqueID()<<endl;
//
//	ms[0]->get(0)->getMolecule()->printDetails();
//	ms[1]->get(0)->getMolecule()->printDetails();
//	t->transform(ms);
//	cout<<endl<<"--------"<<endl;
//	ms[0]->get(0)->getMolecule()->printDetails();
//	ms[1]->get(0)->getMolecule()->printDetails();
//
//	cout<<endl<<endl<<endl;
//	delete [] ms;
//	delete t; delete rl_0;
//	delete rl_1;
//
//	delete s;
}





void NFcore::test_tree()
{
	System *s = new System("boo");



}


void NFcore::test_simple()
{
//	cout<<"Testing rxns..."<<endl;
//
//
//	System *s = new System("boo");
//
//
//	int numOfBsites = 1;
//	string * bSiteNames = new string [numOfBsites];
//	bSiteNames[0] = "y";
//	int numOfStates = 1;
//	string * stateNames = new string [numOfStates];
//	stateNames[0] = "p";
//	int * stateValues = new int [numOfStates];
//	stateValues[0] = 1;
//	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
//	molX->populateWithDefaultMolecules(500);
//	molX->getMolecule(0)->printDetails();
//
//
//
//	TemplateMolecule *xTemp1 = new TemplateMolecule(molX);
//	xTemp1->addStateValue("p",1);
//	TemplateMolecule *xTemp2 = new TemplateMolecule(molX);
//	xTemp2->addStateValue("p",1);
//	//TemplateMolecule::bind(xTemp1,"y",yTemp2,"x");
//	//xTemp2->printDetails();
//
//
//	vector <TemplateMolecule *> templates;
//	templates.push_back( xTemp1 );
//
//
//
//
//
//
//	TransformationSet *t = new TransformationSet(templates);
//	t->addStateChangeTransform(xTemp1,"p",0);
//	//t->addBindingTransform(xTemp1,"y",yTemp2,"x");
//	cout<<endl;
//
//
//
//
//
//	t->finalize();
//
//
//	ReactantList * rl = new ReactantList(0, t, 15);
//	//rl->printDetails();
//
//
//
//
//	for(int i=0; i<12; i++)
//	{
//		MappingSet *ms = rl->pushNextAvailableMappingSet();
//		if(!xTemp1->compare(molX->getMolecule(20+i),ms)) {
//			rl->popLastMappingSet();
//		}
//	}
//
//
//	rl->printDetails();
//	//cout<<"ms->getId() = "<<ms->getId()<<endl;
//	//rl->removeMappingSet(5);
//
//
//	//rl->printDetails();
//
//
//
//	MappingSet **ms = new MappingSet *[2];
//	rl->pickRandom(ms[0]);
//	cout<<"Retrieving randomly: " <<ms[0]->get(0)->getMolecule()->getUniqueID()<<endl;
//
//	ms[0]->get(0)->getMolecule()->printDetails();
//	t->transform(ms);
//	cout<<endl<<"--------"<<endl;
//	ms[0]->get(0)->getMolecule()->printDetails();
//
//	cout<<endl<<endl<<endl;
//	delete t; delete rl;
//	delete s;
}


