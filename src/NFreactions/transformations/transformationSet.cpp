#include "transformationSet.hh"
#include "transformation.hh"

using namespace NFcore;


vector <Molecule *> TransformationSet::deleteList;
vector <Molecule *> TransformationSet::updateAfterDeleteList;
vector <Molecule *>::iterator TransformationSet::it;

TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates) // @suppress("Class members should be properly initialized")
{
	this->hasSymUnbinding=false;
	this->hasSymBinding = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->n_addmol  = 0;

	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	this->addmol = new TemplateMolecule *[n_addmol];

	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [n_reactants];
	finalized = false;
}


TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates, // @suppress("Class members should be properly initialized")
		                             vector <TemplateMolecule *> addMoleculeTemplates )
{
	this->hasSymUnbinding = false;
	this->hasSymBinding   = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	//Remember our add molecules
	this->n_addmol = addMoleculeTemplates.size();
	this->addmol = new TemplateMolecule *[n_addmol];
	for(unsigned int r=0; r<n_addmol; r++)
		this->addmol[r] = addMoleculeTemplates.at(r);

	// complex bookkeeping is off by default
	this->complex_bookkeeping = false;

	// for now, symmetry factor is off by default
	this->useSymmetryFactor = false;
	this->symmetryFactor = 1.0;

	// check collisions is off by default
	this->check_collisions = false;

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [ this->getNmappingSets() ];
	finalized = false;
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<getNmappingSets(); r++)  {
		Transformation *t;
		while(transformations[r].size()>0)
		{
			t = transformations[r].back();
			transformations[r].pop_back();
			delete t;
		}
	}

	Transformation *t;
	while(addMoleculeTransformations.size()>0)
	{
		t = addMoleculeTransformations.back();
		addMoleculeTransformations.pop_back();
		delete t;
	}

	while(addSpeciesTransformations.size()>0)
	{
		t = addSpeciesTransformations.back();
		addSpeciesTransformations.pop_back();
		delete t;
	}

	delete [] transformations;
	delete [] reactants;
	delete [] addmol;

	this->n_reactants = 0;
	this->n_addmol = 0;
}


TemplateMolecule * TransformationSet::getTemplateMolecule( unsigned int reactantIndex ) const
{
	if ( reactantIndex < n_reactants )
	{
		return reactants[reactantIndex];
	}
	else if ( reactantIndex < getNmappingSets() )
	{
		return addmol[reactantIndex-n_reactants];
	}
}



bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, int finalStateValue)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	//cout<<"Adding state change transform to value: "<<finalStateValue<<endl;
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genStateChangeTransform(cIndex, finalStateValue, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

bool TransformationSet::addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addStateChangeTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	Transformation *transformation = TransformationFactory::genLocalFunctionReference(PointerName,scope,t);
	transformations[reactantIndex].push_back(transformation);
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}


bool TransformationSet::addIncrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addIncrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genIncrementStateTransform(cIndex, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}
bool TransformationSet::addDecrementStateTransform(TemplateMolecule *t, string cName)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementStateTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	Transformation *transformation = TransformationFactory::genDecrementStateTransform(cIndex, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}



bool TransformationSet::addStateChangeTransform(TemplateMolecule *t, string cName, string finalStateValue)
{
	int cIndex = t->getMoleculeType()->getCompIndexFromName(cName);
	int fStateValue = t->getMoleculeType()->getStateValueFromName(cIndex,finalStateValue);
	return TransformationSet::addStateChangeTransform(t,cName, fStateValue);
}


bool TransformationSet::addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex2==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	//Find the index of the respective binding sites
	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);


	//Check for symmetric binding
	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
	if( isSymmetric )
		hasSymBinding = true;


	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	if(reactantIndex1==reactantIndex2)
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1, t1);
	else
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size(), t1);

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2, t2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}



bool TransformationSet::addNewMoleculeBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	//Again, first find the reactants that the binding pertains to
	int reactantIndex1 = find(t1);
	int reactantIndex2 = find(t2);
	if(reactantIndex2==-1 || reactantIndex2==-1) {
		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}

	//Find the index of the respective binding sites
	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);


	//Check for symmetric binding
	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
	if( isSymmetric )
		hasSymBinding = true;


	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
	//equal to the size.
	Transformation *transformation1;
	if(reactantIndex1==reactantIndex2)
		transformation1 = TransformationFactory::genNewMoleculeBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1, t1);
	else
		transformation1 = TransformationFactory::genNewMoleculeBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size(), t1);

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2, t2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}




// deprecated
//
//bool TransformationSet::addBindingSeparateComplexTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
//{
//	cout<<"adding separate complex binding"<<endl;
//	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
//	//Again, first find the reactants that the binding pertains to
//	int reactantIndex1 = find(t1);
//	int reactantIndex2 = find(t2);
//	if(reactantIndex2==-1 || reactantIndex2==-1) {
//		cerr<<"Couldn't find one of the templates you gave me!  In transformation set - addBindingTransform!\n";
//		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
//		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
//		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
//		return false;
//	}
//
//	//Find the index of the respective binding sites
//	unsigned int cIndex1 = t1->getMoleculeType()->getCompIndexFromName(bSiteName1);
//	unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);
//
//
//	//Check for symmetric binding
//	bool isSymmetric = TemplateMolecule::checkSymmetry(t1,t2,bSiteName1,bSiteName2);
//	if( isSymmetric )
//		hasSymBinding = true;
//
//
//
//	//Add transformation 1: Note that if both molecules involved with this bond are in the same reactant list, then
//	//the mappingIndex will be size()+1.  But if they are on different reactant lists, then the mappingIndex will be exactly
//	//equal to the size.
//	Transformation *transformation1;
//	if(reactantIndex1==reactantIndex2)
//		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
//	else
//		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());
//
//	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);
//
//	transformations[reactantIndex1].push_back(transformation1);
//	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
//	t1->addMapGenerator(mg1);
//
//	transformations[reactantIndex2].push_back(transformation2);
//	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex1].size()-1);
//	t2->addMapGenerator(mg2);
//
//	return true;
//}



bool TransformationSet::addUnbindingTransform(TemplateMolecule *t, string bSiteName, TemplateMolecule *t2, string bSiteName2)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	TemplateMolecule *tToTransform = 0;
	if(t==0 && t2==0) {
		cerr<<"Error in transformation set! when creating unbinding transform!"<<endl;
		cerr<<"Both molecules you gave me are null!\n";
		return false;
	} else if(t2==0) {
		tToTransform = t;
	} else if(t==0) {
		tToTransform = t2;
	} else {
		// they are both real, so randomly pick t1
		tToTransform=t;

		//Check for symmetric unbinding
		bool isSymmetric = TemplateMolecule::checkSymmetryAroundBond(t,t2,bSiteName,bSiteName2);
		if( isSymmetric ) {
			hasSymUnbinding = true;
		}

	}

	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(tToTransform);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you might get this error."<<endl;
		return false;
	}

	// 2) Create a Transformation object to remember the information
	unsigned int cIndex = tToTransform->getMoleculeType()->getCompIndexFromName(bSiteName);
	Transformation *transformation = TransformationFactory::genUnbindingTransform(cIndex, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 4) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	tToTransform->addMapGenerator(mg);

	// 4) Create an empty transformation for the binding partner so that
	// connectivity can be inferred. Arvind Rasi Subramaniam
	if (t2 != 0) {
		unsigned int cIndex2 = t2->getMoleculeType()->getCompIndexFromName(bSiteName2);
		Transformation *transformation2 = TransformationFactory::genUnbindingTransform2(cIndex2, t2);
		transformations[reactantIndex].push_back(transformation2);
		MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
		t2->addMapGenerator(mg);
	}


	return true;
}





/*!
	Adds a delete rule to the given TemplateMolecule.
	@author Michael Sneddon
*/
bool TransformationSet::addDeleteMolecule(TemplateMolecule *t, int deletionType) {
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDeleteMolecule!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genRemoveMoleculeTransform(deletionType, t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}


/*!
	Adds a decrement population rule to the given TemplateMolecule.
	@author Justin Hogg
*/
bool TransformationSet::addDecrementPopulation(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }
	int reactantIndex = find(t);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set - addDecrementPopulation!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
		return false;
	}
	Transformation *transformation = TransformationFactory::genDecrementPopulationTransform(t);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	t->addMapGenerator(mg);
	return true;
}

bool TransformationSet::addAddSpecies( SpeciesCreator *sc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddSpeciesTransform * transformation = TransformationFactory::genAddSpeciesTransform( sc );

	// 3) Add the transformation object to the TransformationSet
	addSpeciesTransformations.push_back( transformation );

	// 3) No map generators needed for an add species!
	return true;
}


bool TransformationSet::addAddMolecule( MoleculeCreator *mc )
{
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	// We don't need a polymorphic transform because AddTransforms are handled separately!
	//  But we do need to call some methods specific to AddMoleculeTransform.
	//  So we're modified TransformationFactory to return the specific object type  --JUstin
	AddMoleculeTransform * transformation = TransformationFactory::genAddMoleculeTransform( mc, mc->getTemplateMolecule() );

	// 3) Add the transformation object to the TransformationSet
	addMoleculeTransformations.push_back( transformation );

	// 3) No map generators needed for an add molecule!
	return true;
}




int TransformationSet::find(TemplateMolecule *t)
{
	if(finalized) { cerr<<"TransformationSet cannot search for a templateMolecule once it has been finalized!"<<endl; exit(1); }
	int findIndex = -1;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(this->reactants[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	// also check add molecule templates
	for(unsigned int r=0; r<n_addmol; r++)  {
		if(this->addmol[r]->contains(t)) {
			if(findIndex==-1) {
				findIndex = r + n_reactants;
			}
			else {
				cerr<<"Found duplicate template molecule in two reaction lists!!  (in transformationSet)."<<endl;
				exit(1);
			}
		}
	}
	return findIndex;
}


bool TransformationSet::transform(MappingSet **mappingSets)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	/*
	 * NOTE: Check for "null conditions" was moved to ReactionClass::fire. This allows rejection of a reaction 
	 * prior to removing molecules from observables. A general method TransformationSet::checkMolecularity has
	 * been implemented to check for incorrect molecularity or reaction center conflicts. --Justin
	 */


	// addMolecule transforms are applied before other transforms so the molecules exist
	//  for potential modification by other transforms.
	int size = addMoleculeTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			addMoleculeTransformations.at(i)->apply_and_map( mappingSets[n_reactants+i]  );
		}
	}

	// apply addSpecies transforms, so we have all the molecules out there
	size = addSpeciesTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			addSpeciesTransformations.at(i)->apply(NULL,NULL);
		}
	}

	// loop over reactants and added molecules, apply transforms to each
	for(unsigned int r=0; r<getNmappingSets(); r++)
	{
		MappingSet *ms = mappingSets[r];
		for ( unsigned int t=0;  t<transformations[r].size();  t++ )
		{
			if( transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
			{	// handle deletions
				Molecule * mol = ms->get(t)->getMolecule();
				if ( transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
				{	// complex deletion: flag connected molecules for deletion
					mol->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT);
				}
				else
				{	// molecule deletion: flag this molecule for deletion
					deleteList.push_back( mol );
				}
			}
			else
			{	// handle other transforms
				transformations[r].at(t)->apply(ms->get(t), mappingSets);
			}
		}
	}


	//Each molecule that is on the delete list must be dealt with
	Molecule * mol;
	for( it = deleteList.begin(); it!=deleteList.end(); it++)
	{
		mol = *it;
		mol->getMoleculeType()->removeMoleculeFromRunningSystem(mol);
	}
	deleteList.clear();

	return true;
}


bool TransformationSet::checkMolecularity( MappingSet ** mappingSets )
{
	if ( n_reactants < 2 )
	{	// unimolecular, so there's nothing to check
		return true;
	}
	else if ( complex_bookkeeping )
	{	// verify that each reactant pattern points to a unique complex
		complex_ids.clear();
		for ( unsigned int ir = 0;  ir < n_reactants;  ++ir )
		{
			// skip populations
			if ( reactants[ir]->getMoleculeType()->isPopulationType() ) continue;

			complex_id = mappingSets[ir]->getComplexID();
			complex_id_iter = std::find( complex_ids.begin(), complex_ids.end(), complex_id );
			if ( complex_id_iter == complex_ids.end() )
			{
				complex_ids.push_back( complex_id );
			}
			else
			{   // two reactant patterns matched the same complex!
				return false;
			}
		}
		return true;
	}
	else if ( check_collisions )
	{	// we won't do a proper check for molecularity, but we should ensure that mappingSets
		//  point to non-overlapping reaction centers.
		for ( collision_pair_iter = collision_pairs.begin(); collision_pair_iter != collision_pairs.end(); ++collision_pair_iter )
		{
			if ( MappingSet::checkForCollisions( mappingSets[ (*collision_pair_iter).first  ],
					                             mappingSets[ (*collision_pair_iter).second ]  ) )
			{	// reaction centers overlap!
				return false;
			}
		}
		return true;
	}
	else
	{   // do nothing
		return true;
	}
}


bool TransformationSet::getListOfProducts(MappingSet **mappingSets,
		vector <Molecule *> &products, int traversalLimit, bool polymerFlag)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	vector <Molecule *>::iterator molIter;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		// if we are deleting the entire complex, we don't have to track molecules in this complex
		if (mappingSets[r]->hasSpeciesDeletionTransform()) continue;

		//cout<<"Traversing:"<<endl;
		//mappingSets[r]->get(0)->getMolecule()->printDetails();
		//mappingSets[r]->get(0)->getMolecule()->traverseBondedNeighborhood(products,traversalLimit);

		/*
		 * I thought that making sure we don't go over the same molecule multiple
		 * times would make the code faster - but this is rarely used for most rxn
		 * systems, so it is commented out for now.  But actually, we do have to check
		 * because if have the same molecule in here twice, then it can mess up our
		 * observable lists...  --michael */
		else
		{
			// For each of the molecules that we possibly affect, traverse the neighborhood
			// Q: Is it sufficient to just look at the first mapping?
			// A: It should be if the traversal limit is set high enough, at least for
			// all standard reactions.  I'm wondering now, though, if it is enough in
			// all cases where you would use the connected-to syntax.  I think so, but
			// someone should test it.  --michael 9Mar2011
			if (!polymerFlag) {
				Molecule * molecule = mappingSets[r]->get(0)->getMolecule();
				// is this molecule already on the product list?
				if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
				{	// Traverse neighbor and add molecules to list
					molecule->traverseBondedNeighborhood(products,traversalLimit);
					//molecule->traverseBondedNeighborhoodForUpdate(products,traversalLimit);
				}
				continue;
			}

			// modified by rasi to skip traversing bonds
			// these require high traversal limit that slows down the simulation
			// instead add all the molecules in the mapping set
			for (unsigned int i = 0; i < mappingSets[r]->getNumOfMappings(); i++) {
				Molecule * molecule = mappingSets[r]->get(i)->getMolecule();

				// This is for non-polymeric molecules like ribosomes
				if (molecule->getVisitedMolecule() == false &&
						molecule->getMoleculeType()->checkIfPolymer() == false) {
					molecule->getBondedProductsForNonpolymers(products,traversalLimit);
				}
				// Search for bonds within interaction distance
				// of sites that are being changed by the MappingSet
				// Add molecules bonded within this distance.
				// This is done to prevent searching across the full molecule, which is very
				// inefficient for simulating translation on an mRNA

				// Note that both polymer molecules that are mapped as well as non-polymeric
				// molecules that are bonded to polymers are searched. This is because the
				// all state changes of polymeric molecules are not present in the mappingSet.
				// If not, you can accidentally miss searching within the interaction distance
				// of all site changes.
				// Arvind Rasi Subramaniam
				molecule->traversePolymerNeighborhood(products, mappingSets[r]->get(i));
			}
		}
	}

	// reset visitation of molecules in products for next round
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
  		(*molIter)->setVisitedMolecule(false);

	// Next, find added molecules that are treated as populations.
	//  Populations molecules have to be removed from observables, then incremented,
	//  and then added back to the observables (Add molecules treated as particles are handled later)
	vector <AddMoleculeTransform *>::iterator addmol_iter;
	for ( addmol_iter = addMoleculeTransformations.begin();
			addmol_iter != addMoleculeTransformations.end();  ++addmol_iter )
	{
		// get molecule creator
		AddMoleculeTransform * addmol = *addmol_iter;
		if ( !(addmol->isPopulationType()) ) continue;

		// Get the population molecule pointer
		Molecule * molecule = addmol->get_population_pointer();

		// is this molecule already on the product list?
		if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
		{	// Add molecule to list
			products.push_back( molecule );
		}
	}

	//cout<<"All together, we have: "<<products.size()<<endl;
	return true;
}


Molecule * TransformationSet::getPopulationPointer( unsigned int r ) const
{
	return addMoleculeTransformations.at(r)->isPopulationType()
			? addMoleculeTransformations.at(r)->get_population_pointer()
		    : NULL;
}


bool TransformationSet::getListOfAddedMolecules(MappingSet **mappingSets, vector <Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	// Add new molecules (particle type) to the list of products
	vector <Molecule *>::iterator molIter;
	for (unsigned int r=n_reactants; r<getNmappingSets(); r++)
	{
		//For each of the molecules that we possibly affect, traverse the neighborhood
		// NOTE: in this instance, it's okay to only look at the first mapping
		Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

		// Skip populations
		if ( molecule->isPopulationType() ) continue;

		// Is the molecule already in the products list?  If not, add to list.
		if ( std::find( products.begin(), products.end(), molecule ) == products.end() )
		{	// Add molecule to list.
			products.push_back( molecule );
			// NOTE: we don't need to traverse neighbors. All new molecules will be put in this
			//  list separately and old molecules that bind to new molecules will be traversed elsewhere
		}
	}

	return true;
}


MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(!finalized) { cerr<<"TransformationSet cannot generate blank mapping if it is not finalized!"<<endl; exit(1); }
	if( reactantIndex>=getNmappingSets() ) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(mappingSetId, transformations[reactantIndex]);
}

void TransformationSet::finalize()
{
	//Be sure to add at least a blank transformation to every reactant if there is no transformation
	//specified so that we count the reactants even if we don't do anything to it.
	for(unsigned int r=0; r<getNmappingSets(); r++)  {
		if(transformations[r].size()==0) {
			transformations[r].push_back(TransformationFactory::genEmptyTransform());
			MapGenerator *mg = new MapGenerator(transformations[r].size()-1);
			getTemplateMolecule(r)->addMapGenerator(mg);
		}
	}

	// Determine if we need to do any reactant center overlap checks.
	// Currently we check a necessary (but not sufficient) condition for the
	//   possibility of reactant center overlap: are there a common molecule types in
	//   a pair of reactant templates. In the future, we could check a necessary and sufficent condition
	//   (e.g. pattern overlap) to avoid extra work.
	if ( (n_reactants>1)  &&  !complex_bookkeeping )
	{
		vector <TemplateMolecule *> tmList1;
		vector <TemplateMolecule *> tmList2;
		vector <TemplateMolecule *>::iterator tm_iter;

		vector <MoleculeType *> moltypes;
		vector <MoleculeType *>::iterator  found_iter;

		for ( unsigned int ir1 = 0;  ir1 < n_reactants;  ++ir1 )
		{
			// skip populations
			if ( reactants[ir1]->getMoleculeType()->isPopulationType() ) continue;

			tmList1.clear();
			moltypes.clear();

			// get all the template molecules in reactant pattern ir1
			TemplateMolecule::traverse(reactants[ir1], tmList1, false);

			// collect the molecule types included in the pattern
			for ( tm_iter = tmList1.begin(); tm_iter != tmList1.end();  ++tm_iter )
			{
				moltypes.push_back( (*tm_iter)->getMoleculeType() );
			}

			for ( unsigned  int ir2 = ir1 + 1;  ir2 < n_reactants;  ++ir2 )
			{
				// skip populations
				if ( reactants[ir2]->getMoleculeType()->isPopulationType() ) continue;

				tmList2.clear();

				// get all the template molecules in reactant pattern ir2
				TemplateMolecule::traverse(reactants[ir2], tmList2, false);

				// check if any moleculeTypes collide
				for ( tm_iter = tmList2.begin(); tm_iter != tmList2.end();  ++tm_iter )
				{
					found_iter = std::find( moltypes.begin(), moltypes.end(), (*tm_iter)->getMoleculeType() );
					if ( found_iter != moltypes.end() )
					{
						// we found a molecule type that is commont to patterns ir1 and ir2!
						check_collisions = true;
						collision_pairs.push_back( pair<int,int>(ir1,ir2) );
						break;
					}
				}
			}
		}

	}

	finalized = true;
}

bool TransformationSet::checkConnection(ReactionClass * rxn) {
	TemplateMolecule * t1;
	MoleculeType * mt1;
	Transformation * transfn;
	int c1;
	for(unsigned int r=0; r<n_reactants; r++) {
		for (unsigned int i=0; i<transformations[r].size(); i++) {
			transfn = transformations[r].at(i);
			t1 = transfn->getTemplateMolecule();
			if (!t1) continue;
			mt1 = t1->getMoleculeType();
			c1 = transfn->getComponentIndex();
			// If the moleculetype or component is not present in the other reaction,
			// it is not connected
			if (!rxn->areMoleculeTypeAndComponentPresent(mt1, c1)) continue;

			// If the TemplateMolecule is 'incompatible' with any of the reactants
			// or products, then the reaction is not connected
			if (!rxn->isTemplateCompatible(t1)) continue;
			// Both checks passed for one op so return true
			return true;
		}
	}
	// Do the same as above but now for the transformed product templates
	for(unsigned int r=0; r<n_reactants; r++) {
		for (unsigned int i=0; i<transformations[r].size(); i++) {
			transfn = transformations[r].at(i);
			// Use the transformed product here
			t1 = transfn->getTemplateMolecule();
			if (!t1) continue;
			t1 = t1->getMappedPartner();
			if (!t1) continue;
			mt1 = t1->getMoleculeType();
			c1 = transfn->getComponentIndex();
			// If the moleculetype or component is present in the other reaction,
			// it is not connected
			if (!rxn->areMoleculeTypeAndComponentPresent(mt1, c1)) continue;

			// If the TemplateMolecule is 'incompatible' with any of the reactants
			// or products, then the reaction is not connected
			if (!rxn->isTemplateCompatible(t1)) continue;
			// Both checks passed for one op so return true
			return true;
		}
	}
	// Now consider the newly added molecules where components don't matter, but
	// compatibility
	for (unsigned int i=0; i<addMoleculeTransformations.size(); i++) {
		transfn = addMoleculeTransformations.at(i);
		// Use the transformed product here
		t1 = transfn->getTemplateMolecule();
		if (!t1) continue;
		// If the TemplateMolecule is 'incompatible' with any of the reactants
		// or products, then the reaction is not connected
		if (!rxn->isTemplateCompatible(t1)) continue;
		// Both checks passed for one op so return true
		return true;
	}
	// Both checks did not pass for any reactant or product template, so not connected
	return false;
}
