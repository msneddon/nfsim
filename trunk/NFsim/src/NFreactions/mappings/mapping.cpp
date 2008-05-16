


#include "mapping.hh"


using namespace NFcore;





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
		