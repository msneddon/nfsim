


#include "mapping.hh"

using namespace NFcore;


NFcore::Mapping::Mapping(unsigned int type, unsigned int index)
{
	this->type = type;
	this->index = index;
}
NFcore::Mapping::~Mapping()
{
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
		

void NFcore::Mapping::clear()
{
	//Clearing the Mapping only requires us to set the molecule to null
	this->m=NULL;
}
bool NFcore::Mapping::setMolecule(Molecule *m)
{
	this->m = m;
	return true;
}
		