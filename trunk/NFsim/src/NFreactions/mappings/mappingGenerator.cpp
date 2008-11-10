

#include "mappingGenerator.hh"


using namespace NFcore;


MapGenerator::MapGenerator(unsigned int mappingIndex)
{
	this->mappingIndex = mappingIndex;
}

bool MapGenerator::map(MappingSet *mappingSet, Molecule *molecule)
{
	mappingSet->set(mappingIndex,molecule);
	return true;
}
