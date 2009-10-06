
#include "transformationSet.hh"

using namespace NFcore;




list <Molecule *> TransformationSet::deleteList;

TransformationSet::TransformationSet(vector <TemplateMolecule *> reactantTemplates)
{
	this->hasSymUnbinding=false;
	this->hasSymBinding = false;

	//cout<<"creating transformationSet..."<<endl;
	//Remember our reactants
	this->n_reactants = reactantTemplates.size();
	this->reactants = new TemplateMolecule *[n_reactants];
	for(unsigned int r=0; r<n_reactants; r++)
		this->reactants[r] = reactantTemplates.at(r);

	//Set up our transformation vectors
	this->transformations = new vector <Transformation *> [n_reactants];
	finalized = false;
}


TransformationSet::~TransformationSet()
{
	for(unsigned int r=0; r<n_reactants; r++)  {
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

	delete [] transformations;
	delete [] reactants;
	this->n_reactants = 0;
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


	//Check for symmetric binding (should do a better check here!)
	if(t1->getMoleculeType()==t2->getMoleculeType())  {
		if(bSiteName1==bSiteName2) {
			hasSymBinding=true;
			//cout<<"**************************setting sym binding reaction!"<<endl;
		}
	}


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

		//Check for symmetric unbinding (should do a better check here!
		if(t->getMoleculeType()==t2->getMoleculeType())  {
			if(bSiteName==bSiteName2) {
				hasSymUnbinding=true;
				//cout<<"setting sym unbinding reaction!"<<endl;
			}
		}
	}

	// 1) Check that the template molecule is contained in one of the reactant templates we have
	int reactantIndex = find(tToTransform);
	if(reactantIndex==-1) {
		cerr<<"Couldn't find the template you gave me!  In transformation set!"<<endl;
		cerr<<"This might be caused if you declare that two molecules are connected, but you\n";
		cerr<<"don't provide how they are connected.  For instance: if you have declared \n";
		cerr<<" A(b).B(a),( instead of, say, A(b!1).B(a!1) ) you will get this error."<<endl;
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

bool TransformationSet::addAddMolecule(SpeciesCreator *sc) {
	if(finalized) { cerr<<"TransformationSet cannot add another transformation once it has been finalized!"<<endl; exit(1); }

	Transformation *transformation = TransformationFactory::genAddMoleculeTransform(sc);

	// 3) Add the transformation object to the TransformationSet
	addMoleculeTransformations.push_back(transformation);

	// 3) No map generators needed for an add molecule!
	//MapGenerator *mg = new MapGenerator(transformations[reactantIndex].size()-1);
	//t->addMapGenerator(mg);
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
	return findIndex;
}
bool TransformationSet::transform(MappingSet **mappingSets)
{
	if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }

	//cout<<"transforming!"<<endl;
	//list <Molecule *> deleteList;

	for(unsigned int r=0; r<n_reactants; r++)  {
		MappingSet *ms = mappingSets[r];
		for(unsigned int t=0; t<transformations[r].size(); t++)
		{
			if(transformations[r].at(t)->getType()==(int)TransformationFactory::REMOVE) {
				Mapping *m1 = ms->get(t);
				if(transformations[r].at(t)->getRemovalType()==TransformationFactory::COMPLETE_SPECIES_REMOVAL) {
					m1->getMolecule()->traverseBondedNeighborhood(deleteList,ReactionClass::NO_LIMIT);
				} else if (transformations[r].at(t)->getRemovalType()==TransformationFactory::DELETE_MOLECULES) {
					deleteList.push_back(m1->getMolecule());
				}


			} else {
				//cout<<"applying!"<<endl;
				//cout<<transformations[r].at(t)->getType()<<endl;
				transformations[r].at(t)->apply(ms->get(t),mappingSets);
			}
		}
	}

	if(deleteList.size()>0) {
		list <Molecule *>::iterator it;

		//Each molecule that is on the list must be dealt with
		for(it = deleteList.begin(); it!=deleteList.end(); it++) {
			(*it)->getMoleculeType()->removeMoleculeFromRunningSystem((*it));
		}
		deleteList.clear();
	}


	int size = addMoleculeTransformations.size();
	if(size>0) {
		for(int i=0; i<size; i++) {
			addMoleculeTransformations.at(i)->apply(NULL,NULL);
		}
	}
	return true;
}
bool TransformationSet::getListOfProducts(MappingSet **mappingSets, list<Molecule *> &products, int traversalLimit)
{
	//if(!finalized) { cerr<<"TransformationSet cannot apply a transform if it is not finalized!"<<endl; exit(1); }
	bool isPresent = false;
	list <Molecule *>::iterator molIter;
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(mappingSets[r]->hasDeletionTransform()) continue;  //if we are deleting this guy, it doesn't have to get updated

		//cout<<"Traversing:"<<endl;
		//mappingSets[r]->get(0)->getMolecule()->printDetails();
		//mappingSets[r]->get(0)->getMolecule()->traverseBondedNeighborhood(products,traversalLimit);

		/*
		 * I thought that making sure we don't go over the same molecule multiple
		 * times would make the code faster - but this is rarely used for most rxn
		 * systems, so it is commented out for now.  But actually, we do have to check
		 * because if have the same molecule in here twice, then it can mess up our
		 * observable lists... */
		//For each of the molecules that we possibly affect, traverse the neighborhood
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
	//cout<<"All together, we have: "<<products.size()<<endl;
	return true;
}


MappingSet *TransformationSet::generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId)
{
	if(!finalized) { cerr<<"TransformationSet cannot generate blank mapping if it is not finalized!"<<endl; exit(1); }
	if(reactantIndex>=n_reactants) {
		cerr<<"Gave me (a transformation Set) a reactant index that was too high!"<<endl;
		exit(1);
	}
	return new MappingSet(mappingSetId, transformations[reactantIndex]);
}

void TransformationSet::finalize()
{
	//Be sure to add at least a blank transformation to every reactant if there is no transformation
	//specified so that we count the reactants even if we don't do anything to it.
	for(unsigned int r=0; r<n_reactants; r++)  {
		if(transformations[r].size()==0) {
			transformations[r].push_back(TransformationFactory::genEmptyTransform());
			MapGenerator *mg = new MapGenerator(transformations[r].size()-1);
			getTemplateMolecule(r)->addMapGenerator(mg);
		}
	}
	finalized = true;
}
