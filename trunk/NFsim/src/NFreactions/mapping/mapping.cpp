#include "mapping.hh"

#include <map>

using namespace NFcore;
using namespace std;




TemplateMapping::TemplateMapping(TemplateMolecule *tm, unsigned int type, const char *name)
{
	if(type==Mapping::BOND)
		this->index = tm->getMoleculeType()->getBindingSiteIndex(name);
	else if(type==Mapping::STATE)
		this->index = tm->getMoleculeType()->getStateIndex(name);
	else
	{
		cout<<"Creating TemplateMapping, but Type parameter is invalid!!"<<endl;
		exit(1);
	}
	this->moleculeType = tm->getMoleculeType();
	this->type = type;
	
	//Add myself to the templateMolecule so that this mapping will be used
	tm->addTemplateMapping(this);
}



TemplateMapping::TemplateMapping(TemplateMolecule *tm, unsigned int type, unsigned int index)
{
	this->moleculeType = tm->getMoleculeType();
	this->type=type;
	this->index = index;
}



Mapping *TemplateMapping::createNewMapping(Molecule *m)
{
	if(transformation==NULL)
	{
		cout<<"Trying to use TemplateMapping to create a Mapping, but this TemplateMapping has"<<endl;
		cout<<"not been associated with a transformation yet!!"<<endl;
		exit(1);
	}
	
	return new Mapping(m, transformation,type,index);
}


void TemplateMapping::setTransformation(Transformation *t)
{
	this->transformation = t;
}

void TemplateMapping::printDetails()
{
	cout<<"TemplateMapping:  creates mappings onto MoleculeType: "<<moleculeType->getName()<<endl;
	if(this->type==Mapping::BOND)
		cout<<"     Maps to binding site: "<<moleculeType->getBindingSiteName(index)<<endl;
	else if(this->type==Mapping::STATE)
			cout<<"\tMaps to state: "<<moleculeType->getStateName(index)<<endl;
	if(this->transformation!=NULL)
	{
		cout<<"\tIs registered with Transformation number: "<<transformation->getIndex();
		cout<<" of type: "<<transformation->getTypeName()<<endl;
	}
	else
		cout<<"\tIs not registered with a Transformation."<<endl;
	cout<<endl;
}



Mapping::Mapping(Molecule *m, Transformation *t, unsigned int type, unsigned int index)
{
	this->molecule = m;
	this->transformation = t;
	this->type = type;
	this->index = index;
}


void Mapping::printDetails()
{
	Molecule *m = getMolecule();
	cout<<"\tMapping: ";
	if(getType()==Mapping::STATE)
		cout<<"State '"<<m->getMoleculeType()->getStateName(getIndex());
	else if (getType()==Mapping::BOND)
		cout<<"Binding Site '"<<m->getMoleculeType()->getBindingSiteName(getIndex());
	else
		cout<<"!!! Type of Mapping cannot be identified !!!  ";
	cout<<"' of Molecule '"<< m->getUniqueID()<<"' of type '"<<m->getMoleculeTypeName();
	cout<<"' linked to Transformation '"<<getTransformation()->getIndex()<<"'"<<endl;
}







MappingSet::MappingSet()
{
	
}


MappingSet::~MappingSet()
{
	Mapping *m;
	while(mappings.size()>0)
	{
		m = mappings.back();
		mappings.pop_back();
		delete m;
	}
}
				
void MappingSet::add(Mapping * m)
{
	mappings.push_back(m);
}

void MappingSet::start()
{
	//this
	
}

Mapping * MappingSet::getNextMapping()
{
	
	return NULL;
}

int MappingSet::getNumOfMappings()
{
	return mappings.size();
}

Mapping * MappingSet::getMapping(unsigned int index)
{
	try {
		return 0;//mappings.at(index);
	} catch(exception& e) {
		cout<<"In MappingSet::getMapping(int index) index given is out of bounds.  Quitting."<<endl;
		exit(1);
	}
}
				
				
void MappingSet::transform(vector <MappingSet *> &mappingSets)
{
	map <int,Mapping *> mapLookupTable;
	map<int,Mapping *>::iterator mapLookupTableIter;
	Mapping ** mArray;

	//cout<<"Transforming "<<mappingSets.size()<<" MappingSets."<<endl;
	
	vector <MappingSet *>::iterator mapSetIter;
	vector <Mapping *>::iterator mapIter;
	for(mapSetIter = mappingSets.begin(); mapSetIter != mappingSets.end(); mapSetIter++ )
	{
		//cout<<"In this mappingSet there are "<<(*mapSetIter)->mappings.size()<<" Mappings:"<<endl;
		
		for(mapIter = (*mapSetIter)->mappings.begin(); mapIter != (*mapSetIter)->mappings.end(); mapIter++ )
		{
			//(*mapIter)->printDetails();
			//(*mapIter)->getMolecule()->printDetails();
			
			Transformation * t = (*mapIter)->getTransformation();
			unsigned int n = t->getNumOfMappings();
			if(n==1)
			{
				//Apply Transformation because this is a unimolecular reaction
				//cout<<"Unimolecular"<<endl;
				mArray = new Mapping * [1];
				mArray[0]=(*mapIter);
				t->transform(mArray);
				delete [] mArray;
			}
			else if(n==2)
			{
				//This is a bimolecular reaction
				//cout<<"Bimolecular"<<endl;
				
				mapLookupTableIter = mapLookupTable.find(t->getIndex());
				if(mapLookupTableIter==mapLookupTable.end())
					mapLookupTable[t->getIndex()] = (*mapIter);
				else
				{
					mArray = new Mapping *[2];
					mArray[0] = (*mapIter);
					mArray[1] = mapLookupTableIter->second;
					t->transform(mArray);
					delete [] mArray;
					mapLookupTable.erase(mapLookupTableIter);
				}
			}
			else
			{
				cout<<"!!  Transformation is not well formed!  Transformations can now only have 1 or 2"<<endl;
				cout<<"reactants, but this transformation has "<<n<<"!!  Quitting."<<endl;
				(*mapIter)->getTransformation()->printDetails();
				exit(1);
				
			}
			
			//(*mapIter)->getMolecule()->printDetails();
		}
	}
	mapLookupTable.clear();
}



void MappingSet::getPossibleProducts(vector <MappingSet *> &mappingSets, list <Molecule *> &products, unsigned int traversalLimit)
{
	vector <MappingSet *>::iterator mapSetIter;
	vector <Mapping *>::iterator mapIter;
	
	//Loop through the mapping sets....
	for(mapSetIter = mappingSets.begin(); mapSetIter != mappingSets.end(); mapSetIter++ )
	{
	//	cout<<"In this mappingSet there are "<<(*mapSetIter)->mappings.size()<<" Mappings:"<<endl;
			
		//Loop through each of the mappings, and extract out the molecule
		for(mapIter = (*mapSetIter)->mappings.begin(); mapIter != (*mapSetIter)->mappings.end(); mapIter++ )
		{
			//For each of the molecules that we possibly affect, traverse the neighborhood
			Molecule * molecule = (*mapIter)->getMolecule();
			
			bool isPresent=false;
			list <Molecule *>::iterator molIter;
			for( molIter = products.begin(); molIter != products.end(); molIter++ )
			  	if((*molIter)==molecule) { isPresent = true; }
			
			if(!isPresent)
				molecule->traverseBondedNeighborhood(products,traversalLimit);
		}
	}
}


	
void MappingSet::clear()
{
	mappings.clear();
	
}

void MappingSet::printDetails()
{
	cout<<"MappingSet contains: "<<endl;
	
	//Molecule * getMolecule() const {return molecule;};
	//			Transformation * getTransformation() const {return transformation;};
	//			unsigned int getType() const {return type;};
	//			unsigned int getIndex() const {return index;};
	//			void printDetails();
	
	
	vector <Mapping *>::iterator mapIter;
	for(mapIter = mappings.begin(); mapIter != mappings.end(); mapIter++ )
	{
		Molecule *m = (*mapIter)->getMolecule();
		cout<<"\t-";
		if((*mapIter)->getType()==Mapping::STATE)
			cout<<"State '"<<m->getMoleculeType()->getStateName((*mapIter)->getIndex());
		else if ((*mapIter)->getType()==Mapping::BOND)
			cout<<"Binding Site '"<<m->getMoleculeType()->getBindingSiteName((*mapIter)->getIndex());
		else
			cout<<"!!! Type of Mapping cannot be identified !!!  ";
		cout<<"' of Molecule '"<< m->getUniqueID()<<"' of type '"<<m->getMoleculeTypeName();
		cout<<"' linked to Transformation '"<<(*mapIter)->getTransformation()->getIndex()<<"'"<<endl;
		
	}
	
	
}





