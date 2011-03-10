
#include "transformationSet.hh"

using namespace NFcore;




list <Molecule *> TransformationSet::deleteList;
list <Molecule *> TransformationSet::updateAfterDeleteList;
list <Molecule *>::iterator TransformationSet::it;

TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates)
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

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [n_reactants];
	finalized = false;
}


TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates,
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


TemplateMolecule *
TransformationSet::getTemplateMolecule( unsigned int reactantIndex ) const
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
	Transformation *transformation = TransformationFactory::genStateChangeTransform(cIndex, finalStateValue);

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
	Transformation *transformation = TransformationFactory::genIncrementStateTransform(cIndex);

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
	Transformation *transformation = TransformationFactory::genDecrementStateTransform(cIndex);

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
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
	else
		transformation1 = TransformationFactory::genBindingTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex2].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}




bool TransformationSet::addBindingSeparateComplexTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2)
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
		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size()+1);
	else
		transformation1 = TransformationFactory::genBindingSeparateComplexTransform1(cIndex1, reactantIndex2, transformations[reactantIndex2].size());

	Transformation *transformation2 = TransformationFactory::genBindingTransform2(cIndex2);

	transformations[reactantIndex1].push_back(transformation1);
	MapGenerator *mg1 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t1->addMapGenerator(mg1);

	transformations[reactantIndex2].push_back(transformation2);
	MapGenerator *mg2 = new MapGenerator(transformations[reactantIndex1].size()-1);
	t2->addMapGenerator(mg2);

	return true;
}



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
	Transformation *transformation = TransformationFactory::genUnbindingTransform(cIndex);

	// 3) Add the transformation object to the TransformationSet
	transformations[reactantIndex].push_back(transformation);

	// 3) Create a MapGenerator object and add it to the templateMolecule
	MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	tToTransform->addMapGenerator(mg);

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
	Transformation *transformation = TransformationFactory::genRemoveMoleculeTransform(deletionType);

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
	Transformation *transformation = TransformationFactory::genDecrementPopulationTransform();

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
	AddMoleculeTransform * transformation = TransformationFactory::genAddMoleculeTransform( mc );

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

	// addMolecule transforms applied before other transforms so the molecules exist
	//  for potential modification by other transforms.  --Justin
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

		//If there is a species removal, we have to do this first
		if(ms->hasSpeciesDeletionTransform()) {
			for ( unsigned int t=0;  t<transformations[r].size();  t++ )
			{
				if( transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
				{
					Molecule * mol = ms->get(t)->getMolecule();
					if ( transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
					{   // handle species deletion
						mol->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT);
					}
				}
			}
		}
		for ( unsigned int t=0;  t<transformations[r].size();  t++ )
		{
			if( transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE )
			{
				// handle molecule deletion
				Molecule * mol = ms->get(t)->getMolecule();
				if ( !transformations[r].at(t)->getRemovalType()==(int)TransformationFactory::COMPLETE_SPECIES_REMOVAL )
				{
					deleteList.push_back( mol );
				}
			}
			else
			{
				//cout<<transformations[r].at(t)->getType()<<endl;
				//ms->printDetails();
				//cout<<"here"<<endl;
				transformations[r].at(t)->apply(ms->get(t),mappingSets);
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


	// Made some changes here:  we identify neighbor molecules that need to be updated following
	//  a delete molecule transform in the usual way (getListOfProducts), so it's no longer
	//  necessary to update reaction membership here. This change was made because (I think)
	//  we weren't updating observables associated with these neighbor molecules.
	//--Justin, 8Mar2010

	//Update anything we have to update because we deleted some of its connected friends
	// (needed when we delete some, but not all molecules, in a species)
	//for( it = updateAfterDeleteList.begin(); it != updateAfterDeleteList.end(); it++ )
	//{
	//	if(!(*it)->isAlive()) { continue; } // skip over molecules that we just removed...  don't actually need this check
	//  	(*it)->updateRxnMembership();
	//  	(*it)->updateTypeIIFunctions();
	//  	(*it)->updateDORRxnValues();
	//}
	//updateAfterDeleteList.clear();

	return true;
}


bool TransformationSet::checkMolecularity( MappingSet ** mappingSets )
{
	// no need to check if fewer than 2 reactants or complex bookkeeping is off
	if ( (n_reactants < 2)  ||  !complex_bookkeeping )
		return true;

	// otherwise we need to make sure each mappingSet is targeting a different complex
	complex_ids.clear();
	for ( unsigned int ir = 0;  ir < n_reactants;  ++ir )
	{
		insert_retval = complex_ids.insert( mappingSets[ir]->getComplexID() );
		if ( insert_retval.second==false ) return false;
	}

	return true;
}


bool TransformationSet::getListOfProducts(MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }
	bool isPresent = false;
	list <Molecule *>::iterator molIter;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		// if we are deleting this guy, it doesn't have to get updated...
		// NOTE: this only applies to Species deletion!!
		if(mappingSets[r]->hasSpeciesDeletionTransform())
		{
			continue;
		}
		//cout<<"Traversing:"<<endl;
		//mappingSets[r]->get(0)->getMolecule()->printDetails();
		//mappingSets[r]->get(0)->getMolecule()->traverseBondedNeighborhood(products,traversalLimit);

		/*
		 * I thought that making sure we don't go over the same molecule multiple
		 * times would make the code faster - but this is rarely used for most rxn
		 * systems, so it is commented out for now.  But actually, we do have to check
		 * because if have the same molecule in here twice, then it can mess up our
		 * observable lists... */
		else
		{
			// For each of the molecules that we possibly affect, traverse the neighborhood
			// TODO:Is it sufficient to just look at the first mapping???  --Justin, 8Mar2011
			// It should be if the traversal limit is set high enough, at least for
			// all standard reactions.  I'm wondering now, though, if it is enough in
			// all cases where you would use the connected-to syntax.  I think so, but
			// someone should test it.  --michael 9Mar2011
			Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

			isPresent=false;
			for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
				if((*molIter)==molecule) { isPresent = true; break;}
			}

			if(!isPresent)
			{
				molecule->traverseBondedNeighborhood(products,traversalLimit);
				//molecule->traverseBondedNeighborhoodForUpdate(products,traversalLimit);
			}
		}
	}

	// Next, find added molecules that are treated as populations.
	//  Populations molecules have to be removed from observables, then incremented,
	//  and then added back to the observables
	// (Add molecules treated as particles are handled later)
	isPresent = false;
	vector <AddMoleculeTransform *>::iterator addmol_iter;
	for ( addmol_iter = addMoleculeTransformations.begin();
			addmol_iter != addMoleculeTransformations.end();  ++addmol_iter )
	{
		// get molecule creator
		AddMoleculeTransform * addmol = *addmol_iter;
		if ( !(addmol->isPopulationType()) ) continue;

		// Get the population molecule pointer
		Molecule * molecule = addmol->get_population_pointer();

		// See if its already in the products list
		isPresent=false;
		for( molIter = products.begin(); molIter != products.end(); molIter++ )
		{
			if ( (*molIter)==molecule ) { isPresent = true; break; }
		}

		// Add to products list, if not already there..
		if(!isPresent) products.push_back( molecule );
	}

	//cout<<"All together, we have: "<<products.size()<<endl;
	return true;
}


bool TransformationSet::getListOfAddedMolecules(MappingSet **mappingSets, list <Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	// Add new molecules (particle type) to the list of products
	bool isPresent = false;
	list <Molecule *>::iterator molIter;
	for (unsigned int r=n_reactants; r<getNmappingSets(); r++)
	{
		//For each of the molecules that we possibly affect, traverse the neighborhood
		// NOTE: in this instance, it's okay to only look at the first mapping
		//  Molecule
		Molecule * molecule = mappingSets[r]->get(0)->getMolecule();

		// Is this a population?
		if ( molecule->isPopulationType() ) continue;

		isPresent=false;
		for( molIter = products.begin(); molIter != products.end(); molIter++ )
		{
			if ((*molIter)==molecule) { isPresent = true; break;}
		}

		if(!isPresent)
		{
			products.push_back( molecule );
			// Pretty sure we don't need to traverse neighbors --Justin
		}

	}
	//cout<<"All together, we have: "<<products.size()<<endl;
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
	finalized = true;
}
