#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;


Observable::Observable(const char* aliasName, TemplateMolecule * templateMolecule)
{
	this->aliasName = aliasName;
	this->templateMolecule = templateMolecule;
	this->count = 0;
}

Observable::~Observable()
{	
	templateMolecule = 0;
}

bool Observable::isObservable(Molecule * m) const { 
	
	//cout<<"comparing obs: "<<this->aliasName<<" to: "<<endl;
//	m->printDetails();
	bool answer = templateMolecule->compare(m);
	//if(answer) cout<<"Match."<<endl;
	//else cout<<"nope"<<endl;
	
	return answer; 
	

}