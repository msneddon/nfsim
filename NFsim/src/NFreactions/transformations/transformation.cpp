

#include "transformation.hh"

using namespace NFcore;



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
StateChangeTransform::StateChangeTransform(int cIndex, int newValue) :
	Transformation(TransformationFactory::STATE_CHANGE)
{
	this->cIndex = cIndex;
	this->newValue = newValue;
}
void StateChangeTransform::apply(Mapping *m, MappingSet **ms)
{
	m->getMolecule()->setComponentState(cIndex,newValue);
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
BindingTransform::BindingTransform(int cIndex, int otherReactantIndex, int otherMappingIndex) :
	Transformation(TransformationFactory::BINDING)
{
	this->cIndex=cIndex;
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
	} else {
		System::NULL_EVENT_COUNTER++;
	}
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
UnbindingTransform::UnbindingTransform(int cIndex) :
	Transformation(TransformationFactory::UNBINDING)
{
	this->cIndex=cIndex;
}
void UnbindingTransform::apply(Mapping *m, MappingSet **ms)
{   //cout<<"unbinding.."<<endl;
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

