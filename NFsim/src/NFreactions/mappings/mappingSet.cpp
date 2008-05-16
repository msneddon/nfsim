


#include "mappingSet.hh"


using namespace NFcore;




MappingSet::MappingSet(unsigned int id, vector <Transformation *> &transformations)
{
	this->id = id;
	this->isSet = false;
	this->n_mappings = transformations.size();
	this->mappings = new Mapping *[n_mappings];
	
	for(unsigned int t=0; t<n_mappings; t++) {
		mappings[t] = new Mapping(transformations.at(t)->getType(), transformations.at(t)->getStateOrSiteIndex() );
	}
}
MappingSet::~MappingSet()
{
	for(unsigned int t=0; t<n_mappings; t++) {
		delete mappings[t];
	}
	delete [] mappings;
}

bool MappingSet::set(unsigned int mappingIndex, Molecule *m){
				mappings[mappingIndex]->setMolecule(m);
				return true;
			}

// Defined inline			
//bool NFcore::MappingSet::set(unsigned int mappingIndex, Molecule *m)
//{
//	mappings[mappingIndex]->setMolecule(m);
//	return true;
//}
//Mapping *NFcore::MappingSet::get(unsigned int mappingIndex)
//{
//	return mappings[mappingIndex];
//}
//bool NFcore::MappingSet::clear()
//{
//	return true;
//}



//  Slower version, but good for debuggin
//
//NFcore::MappingSet::MappingSet(unsigned int id, vector <Transformation *> &transformations)
//{
//	this->id = id;
//	this->isSet = false;
//	this->n_mappings = transformations.size();
//	this->mappings = new Mapping *[n_mappings];
//	
//	for(unsigned int t=0; t<n_mappings; t++) {
//		mappings[t] = new Mapping(transformations.at(t)->getType(), transformations.at(t)->getStateOrSiteIndex() );
//	}
//}
//NFcore::MappingSet::~MappingSet()
//{
//	for(unsigned int t=0; t<n_mappings; t++) {
//		delete mappings[t];
//	}
//	delete [] mappings;
//}
//			
//bool NFcore::MappingSet::set(unsigned int mappingIndex, Molecule *m)
//{
//	if(mappingIndex>=n_mappings) {
//		cerr<<"Out of bounds access to a mapping in set function of mapping set!"<< mappingIndex<<" but max is " << n_mappings<<endl;
//		return false;
//	}
//	mappings[mappingIndex]->setMolecule(m);
//	isSet = true;
//	return true;
//}
//Mapping *NFcore::MappingSet::get(unsigned int mappingIndex)
//{
//	if(mappingIndex>=n_mappings) {
//		cerr<<"Out of bounds access to a mapping in get function of a mapping set!"<<endl;
//		return false;
//	}
//	return mappings[mappingIndex];
//}
//bool NFcore::MappingSet::clear()
//{
//	isSet = false;
//	
//	//This is not entirely necessary, but for now can be used for debugging
//	for(unsigned int t=0; t<n_mappings; t++) {
//		mappings[t]->clear();
//	}
//	
//	return true;
//}
	




