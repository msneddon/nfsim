

#include "transformation.hh"

using namespace NFcore;



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
StateChangeTransform::StateChangeTransform(int stateIndex, int newStateValue) :
	Transformation(TransformationFactory::STATE_CHANGE)
{
	this->stateIndex = stateIndex;
	this->newStateValue = newStateValue;
}
void StateChangeTransform::apply(Mapping *m, MappingSet **ms)
{
	m->getMolecule()->setState(stateIndex,newStateValue);
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
BindingTransform::BindingTransform(int siteIndex, int otherReactantIndex, int otherMappingIndex) :
	Transformation(TransformationFactory::BINDING)
{
	this->siteIndex=siteIndex;
	this->otherReactantIndex=otherReactantIndex;
	this->otherMappingIndex=otherMappingIndex;
}

void BindingTransform::apply(Mapping *m, MappingSet **ms) 
{
	Mapping *m2 = ms[this->otherReactantIndex]->get(this->otherMappingIndex);
	Molecule::bind(m->getMolecule(),m->getIndex(), m2->getMolecule(), m2->getIndex());
	
}
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


void BindingSeparateComplexTransform::apply(Mapping *m, MappingSet **ms) 
{
	Mapping *m2 = ms[this->otherReactantIndex]->get(this->otherMappingIndex);
	//cout<<"complex ID: "<<m->getMolecule()->getComplexID()<<" "<<m2->getMolecule()->getComplexID()<<endl;
	if(m->getMolecule()->getComplexID()!=m2->getMolecule()->getComplexID()) {
		Molecule::bind(m->getMolecule(),m->getIndex(), m2->getMolecule(), m2->getIndex());
	}
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
UnbindingTransform::UnbindingTransform(int siteIndex) :
	Transformation(TransformationFactory::UNBINDING)
{
	this->siteIndex=siteIndex;
}
void UnbindingTransform::apply(Mapping *m, MappingSet **ms)
{
	Molecule::unbind(m->getMolecule(),m->getIndex());
}
	





///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
AddMoleculeTransform::AddMoleculeTransform(SpeciesCreator *sc) :
	Transformation(TransformationFactory::ADD)
{
	this->sc=sc;
}
AddMoleculeTransform::~AddMoleculeTransform()
{
	delete sc;
}
void AddMoleculeTransform::apply(Mapping *m, MappingSet **ms)
{
	this->sc->create();
}
	


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void RemoveMoleculeTransform::apply(Mapping *m, MappingSet **ms)
{
	cout<<"!! Warning: calling apply on a RemoveMoleculeTransform!  This cannot be handled here!"<<endl;
	cout<<"!! This function should not be called.  The TransformationSet object should handle this!"<<endl;
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
NFcore::Transformation * TransformationFactory::genEmptyTransform()
{
	return new EmptyTransform();
}
NFcore::Transformation * TransformationFactory::genStateChangeTransform(unsigned int stateIndex, int newStateValue)
{
	return new StateChangeTransform(stateIndex, newStateValue);
}
NFcore::Transformation * TransformationFactory::genBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex)
{	
	return new BindingTransform(bSiteIndex, otherReactantIndex, otherMappingIndex);
}
NFcore::Transformation * TransformationFactory::genBindingSeparateComplexTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex)
{	
	return new BindingSeparateComplexTransform(bSiteIndex, otherReactantIndex, otherMappingIndex);
}

NFcore::Transformation * TransformationFactory::genBindingTransform2(unsigned int bSiteIndex)
{
	return new EmptyTransform(bSiteIndex);
}
NFcore::Transformation * TransformationFactory::genUnbindingTransform(unsigned int bSiteIndex)
{
	return new UnbindingTransform(bSiteIndex);
}
NFcore::Transformation * TransformationFactory::genAddMoleculeTransform(SpeciesCreator *sc)
{
	return new AddMoleculeTransform(sc);
}
NFcore::Transformation * TransformationFactory::genRemoveMoleculeTransform()
{
	return new RemoveMoleculeTransform();
}
