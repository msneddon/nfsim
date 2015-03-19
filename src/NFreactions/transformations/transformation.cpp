

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
IncrementStateTransform::IncrementStateTransform(unsigned int cIndex) :
	Transformation(TransformationFactory::INCREMENT_STATE)
{
	this->cIndex = cIndex;
}
void IncrementStateTransform::apply(Mapping *m, MappingSet **ms)
{
	int oldValue = m->getMolecule()->getComponentState(cIndex);
	m->getMolecule()->setComponentState(cIndex,oldValue+1);
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
DecrementStateTransform::DecrementStateTransform(unsigned int cIndex) :
	Transformation(TransformationFactory::DECREMENT_STATE)
{
	this->cIndex = cIndex;
}
void DecrementStateTransform::apply(Mapping *m, MappingSet **ms)
{
	int oldValue = m->getMolecule()->getComponentState(cIndex);
	m->getMolecule()->setComponentState(cIndex,oldValue-1);
}




///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

NewMoleculeBindingTransform::NewMoleculeBindingTransform(int cIndex, int otherReactantIndex, int otherMappingIndex) :
		BindingTransform(cIndex, otherReactantIndex, otherMappingIndex)
{ }

bool NewMoleculeBindingTransform::checkForNullCondition(Mapping *m, MappingSet **ms)
{
	// binding to a new molecule cannot break a null condition (because null conditions are thus far
	// dependent on molecularity, which must be correct if the molecule did not exist in the reactant pattern.
	return false;
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
	//cout<<" cIndex: "<<cIndex;
	//cout<<" otherReactantIndex: "<<otherReactantIndex;
	//cout<<" otherMappingIndex: "<<otherMappingIndex;

	Mapping *m2 = ms[this->otherReactantIndex]->get(this->otherMappingIndex);

	//Currently, this is set to block all binding events that happen internally to a single
	//molecule.  I think this is reasonable to do...
	// (Intra-molecular binding is probably ok. BNGL supports it. --Justin.)
	// (this is commented out, and so nfsim now supports internal binding events --michael)
	//
	//if(m->getMolecule()->getUniqueID()==m2->getMolecule()->getUniqueID()) { // && m->getIndex() == m2->getIndex()) {
	//	System::NULL_EVENT_COUNTER++;
	//} else {
	Molecule::bind(m->getMolecule(),m->getIndex(), m2->getMolecule(), m2->getIndex());
}


bool BindingTransform::checkForNullCondition(Mapping *m, MappingSet **ms)
{
	// the null condition in this case is if a molecule is trying to bind a site to the same site on itself!
	Mapping *m2 = ms[this->otherReactantIndex]->get(this->otherMappingIndex);
	if(m->getMolecule()->getUniqueID()==m2->getMolecule()->getUniqueID() && m->getIndex() == m2->getIndex())
	{
		System::NULL_EVENT_COUNTER++;
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// deprecated
//
//void BindingSeparateComplexTransform::apply(Mapping *m, MappingSet **ms)
//{
//	cerr<<"Using BindingSeparateComplexTransform!!  This transformation is deprecated in v1.09+!"<<endl;
//	exit(1);
//
//	Mapping *m2 = ms[this->otherReactantIndex]->get(this->otherMappingIndex);
//	//cout<<"complex ID: "<<m->getMolecule()->getComplexID()<<" "<<m2->getMolecule()->getComplexID()<<endl;
//
//	if(m->getMolecule()->getComplexID()!=m2->getMolecule()->getComplexID()) {
//		Molecule::bind(m->getMolecule(),m->getIndex(), m2->getMolecule(), m2->getIndex());
//	} else {
//		System::NULL_EVENT_COUNTER++;
//	}
//
//}


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
AddSpeciesTransform::AddSpeciesTransform(SpeciesCreator *sc) :
	Transformation(TransformationFactory::ADD)
{
	this->sc=sc;
}


AddSpeciesTransform::~AddSpeciesTransform()
{
	delete sc;
}


void AddSpeciesTransform::apply(Mapping *m, MappingSet **ms)
{
	this->sc->create();
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
AddMoleculeTransform::AddMoleculeTransform( MoleculeCreator * _mc )
 : Transformation( TransformationFactory::ADD )
{
	this->mc = _mc;
	new_molecule = NULL;
}


AddMoleculeTransform::~AddMoleculeTransform()
{
	delete mc;
}


bool
AddMoleculeTransform::isPopulationType() const
{
	return mc->isPopulationType();
};


// get pointer to population molecule
Molecule *
AddMoleculeTransform::get_population_pointer() const
{
	return mc->get_population_pointer();
};


void AddMoleculeTransform::apply_and_map(MappingSet *ms)
{
	// create molecule and get pointer
	new_molecule = this->mc->create_molecule();

	// point mappings to the new molecule
	unsigned int n_mappings = ms->getNumOfMappings();
	for ( unsigned int im = 0;  im < n_mappings;  ++im )
	{
		ms->set( im, new_molecule );
	}
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
RemoveMoleculeTransform::RemoveMoleculeTransform(int removalType) :
	Transformation(TransformationFactory::REMOVE) {
	this->removalType=removalType;
}


void RemoveMoleculeTransform::apply(Mapping *m, MappingSet **ms)
{
	cout << "!! Warning: calling apply from a RemoveMoleculeTransform!"
	     << "!! This cannot be handled here! The TransformationSet object should handle this!" << endl;
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
DecrementPopulationTransform::DecrementPopulationTransform() :
	Transformation(TransformationFactory::DECREMENT_POPULATION)
{
	this->cIndex = -1;
}
void DecrementPopulationTransform::apply(Mapping *m, MappingSet **ms)
{
	m->getMolecule()->decrementPopulation();
}


/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
LocalFunctionReference::LocalFunctionReference(string PointerName, int scope, TemplateMolecule *tm)
	: Transformation(TransformationFactory::LOCAL_FUNCTION_REFERENCE) {

	this->PointerName=PointerName;
	this->scope=scope;
	this->tm=tm;

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
NFcore::Transformation * TransformationFactory::genNewMoleculeBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex)
{
	return new NewMoleculeBindingTransform(bSiteIndex, otherReactantIndex, otherMappingIndex);
}

//NFcore::Transformation * TransformationFactory::genBindingSeparateComplexTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex)
//{
//	return new BindingSeparateComplexTransform(bSiteIndex, otherReactantIndex, otherMappingIndex);
//}

NFcore::Transformation * TransformationFactory::genBindingTransform2(unsigned int bSiteIndex)
{
	return new EmptyTransform(bSiteIndex);
}
NFcore::Transformation * TransformationFactory::genUnbindingTransform(unsigned int bSiteIndex)
{
	return new UnbindingTransform(bSiteIndex);
}
NFcore::AddSpeciesTransform * TransformationFactory::genAddSpeciesTransform(SpeciesCreator *sc)
{
	return new AddSpeciesTransform(sc);
}
NFcore::AddMoleculeTransform * TransformationFactory::genAddMoleculeTransform(MoleculeCreator *mc)
{
	return new AddMoleculeTransform(mc);
}
NFcore::Transformation * TransformationFactory::genRemoveMoleculeTransform(int removalType)
{
	return new RemoveMoleculeTransform(removalType);
}


NFcore::Transformation * TransformationFactory::genIncrementStateTransform(unsigned int cIndex)
{
	return new IncrementStateTransform(cIndex);
}
NFcore::Transformation * TransformationFactory::genDecrementStateTransform(unsigned int cIndex)
{
	return new DecrementStateTransform(cIndex);
}

NFcore::Transformation * TransformationFactory::genDecrementPopulationTransform()
{
	return new DecrementPopulationTransform();
}


Transformation * TransformationFactory::genLocalFunctionReference(string PointerName, int type, TemplateMolecule *tm)
{
	return new LocalFunctionReference(PointerName, type, tm);
}












