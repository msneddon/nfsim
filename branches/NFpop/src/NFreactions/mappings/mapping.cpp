


#include "mapping.hh"

using namespace NFcore;


NFcore::Mapping::Mapping(unsigned int type, int index)
{
	this->type = type;
	this->index = index;
	this->m=NULL;
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
int NFcore::Mapping::getIndex() const
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


void NFcore::Mapping::printDetails() const { printDetails(cout); }
void NFcore::Mapping::printDetails(ostream &o) const
{
	o<<"M("<<index<<","<<type<<"): mapped to: ";
	//cout<<"here"<<endl;
	if(m!=NULL) {
		o<<m->getMoleculeTypeName()<<"_"<<m->getUniqueID()<<"  ";
		m->printDetails();
	}
	else o<<"nothing.";

}


void NFcore::Mapping::clone(Mapping *original, Mapping *newClone)
{
	//Before we clone, we have to make sure the newClone Mapping
	//is the right type and can accept the clone
	if(original->index!=newClone->index || original->type!=newClone->type) {
		cerr<<"Error in Mapping!! : When cloning an existing Mapping into a new Mapping,\n";
		cerr<<"the new Mapping has a different index and/or type!  That means you cannot\n";
		cerr<<"clone the original onto this Mapping!"<<endl;
		exit(1);
	}

	//Assign properly the molecule
	newClone->m=original->m;
}


