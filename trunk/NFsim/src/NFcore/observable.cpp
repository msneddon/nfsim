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

