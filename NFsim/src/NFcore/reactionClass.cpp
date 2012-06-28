


#include "NFcore.hh"


using namespace std;
using namespace NFcore;



ReactionClass::ReactionClass(string name, double baseRate, string baseRateParameterName, TransformationSet *transformationSet, System *s)
{
	//cout<<"\n\ncreating reaction "<<name<<endl;
	this->system=s;
	this-> tagged = false;

	totalRateFlag=false;
	isDimerStyle=false;
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->baseRate = baseRate;
	this->baseRateParameterName=baseRateParameterName;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;


	//Set up the template molecules from the transformationSet
	this->n_reactants   = transformationSet->getNreactants();
	this->n_mappingsets = transformationSet->getNmappingSets();
	this->reactantTemplates = new TemplateMolecule *[n_reactants];
	vector <TemplateMolecule*> tmList;
	vector <int> hasMapGenerator;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		//The main reactant should be the one that is getting modified...
		//In other words, we select the reactant that has at least one map generator, and
		//to minimize mistakes, with the least sym sites...
		TemplateMolecule *curTemplate = transformationSet->getTemplateMolecule(r);
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//First, single out all the templates that have at least one map generator
		for(unsigned int i=0; i<tmList.size(); i++) {
			//cout<<"looking at:"<<endl;
			//tmList.at(i)->printDetails();

			if(tmList.at(i)->getN_mapGenerators()>0) {
				hasMapGenerator.push_back(i);
			}
		}

		//Find the one with the least sym comp bonds...
		int minSymSites = 999999;
		for(unsigned int k=0; k<hasMapGenerator.size(); k++) {
			if(tmList.at(hasMapGenerator.at(k))->getN_symCompBonds()<minSymSites) {
				curTemplate = tmList.at(hasMapGenerator.at(k));
				minSymSites = curTemplate->getN_symCompBonds();
			}
		}

		reactantTemplates[r] = curTemplate;
		tmList.clear(); hasMapGenerator.clear();
	}
	mappingSet = new MappingSet *[n_mappingsets];



	/* create blank mappingSets for the added molecules. These will be used
	 * to hold mappings to added molecules, which is useful for rules that create
	 *  molecules and then perform other transformations.  --Justin, 1Mar2011
	 */
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		mappingSet[r] = transformationSet->generateBlankMappingSet(r,0);
	}



	//Here, if we identify that there are disjoint sets in this pattern, from
	//the connected-to syntax, then we have to flag the ones that we actually
	//have to traverse down...
	for(unsigned int r=0; r<n_reactants; r++)
	{
		tmList.clear();

		// Get the connected set of molecules
		TemplateMolecule *curTemplate = reactantTemplates[r];
		TemplateMolecule::traverse(curTemplate,tmList,TemplateMolecule::FIND_ALL);

		//Label the unique sets, and only continue if we have more than one set
		vector <vector <TemplateMolecule *> > sets;
		vector <int> uniqueSetId;
		int setCount = TemplateMolecule::getNumDisjointSets(tmList,sets,uniqueSetId);
		if(setCount<=1) continue;


		//count the number of map generators (rxn centers) in each set
		vector <int> numMapGenerators;
		for(int s=0; s<setCount; s++) { numMapGenerators.push_back(0); }

		int curTemplateSetId = -1;
		for(unsigned int t=0;t<tmList.size();t++) {
			if(tmList.at(t)==curTemplate) {
				curTemplateSetId = uniqueSetId.at(t);
			}
			int n_maps = numMapGenerators.at(uniqueSetId.at(t));
			numMapGenerators.at(uniqueSetId.at(t)) = n_maps+tmList.at(t)->getN_mapGenerators();
		}

		// Debug output
		//cout<<"found "<<setCount<<" unique sets."<<endl;
		//cout<<"found that reactant molecule head is in set: "<<curTemplateSetId<<endl;
		//for(int s=0; s<setCount; s++) {
		//	cout<<"set: "<<s<<" has "<<numMapGenerators.at(s)<<" map generators."<<endl;
		//}
		//for(unsigned int i=0; i<tmList.size(); i++) {
		//	cout<<"looking at:"<<endl;
		//	tmList.at(i)->printDetails();
		//}


		//Lets rearrange the connected-to elements so that the one head is listed as
		//connected to all other molecules.  This will better suit our needs.

		// first, clear out the old connections
		for(unsigned int i=0; i<tmList.size(); i++) {
			tmList.at(i)->clearConnectedTo();
		}

		// add back the connections, but always through the head template
		int rxnCenterSets = 1;
		int curSet=0;
		for(unsigned int i=0; i<uniqueSetId.size(); i++) {
			if(uniqueSetId.at(i)==curTemplateSetId) {
				if(curSet==curTemplateSetId) curSet++;
				continue;
			}
			if(uniqueSetId.at(i)==curSet) {
				bool otherHasRxnCenter = false;
				if(numMapGenerators.at(curSet)>0) {
					otherHasRxnCenter=true;
					rxnCenterSets++;
				}
				TemplateMolecule *otherTemplate = tmList.at(i);
				int ctIndex1=curTemplate->getN_connectedTo();
				int ctIndex2=otherTemplate->getN_connectedTo();
				curTemplate->addConnectedTo(otherTemplate,ctIndex2,otherHasRxnCenter);
				otherTemplate->addConnectedTo(curTemplate,ctIndex1);
				curSet++;
			}
		}

		if(rxnCenterSets>2) {
			cout.flush();
			cerr<<"\n\n   Error in Reaction Rule: "<<name<<endl;
			cerr<<"   You created a reaction with a pattern that includes the connected-to\n";
			cerr<<"   syntax (ie: A().B()).  You included 3 or more disjoint sets of molecules\n";
			cerr<<"   where there are more than 2 sets with rxn centers.  This may work ok, \n";
			cerr<<"   but you really shouldn't ever do something this crazy, so I'm just going\n";
			cerr<<"   to stop you now.  Goodbye.\n"<<endl;
			exit(1);
		}


		//cout<<"++++++++++++++++"<<endl;
		//for(unsigned int i=0; i<tmList.size(); i++) {
		//	tmList.at(i)->printDetails();
		//}


		//Finally, clear out the data structures.
		for(unsigned int i=0; i<sets.size(); i++) sets.at(i).clear();
		sets.clear(); uniqueSetId.clear();
		numMapGenerators.clear();
	}


	//Check here to see if we have molecule types that are the same across different reactants
	//Because if so, we will give a warning
	if(n_reactants>2) cerr<<"Warning!! You created a reaction ("<< name <<") that has more than 2 reactants.  This has not been extensively tested!"<<endl;

	if(n_reactants==2)
	{	//If the reactants are of the same type, then we have to make a few special considerations
		if(reactantTemplates[0]->getMoleculeType()==reactantTemplates[1]->getMoleculeType())
		{
			cout<<endl;
			cout<<"Warning! You have a binding rxn (" << name << ") that allows a moleculeType to bind another of the same type."<<endl;
			cout<<"Make sure that is correct, because this can potentially make long polymers or large aggregates."<<endl;
			cout<<endl;
		}
	}

	if ( this->transformationSet->usingSymmetryFactor() )
	{	// new general method for handling reaction center symmetry
		baseRate *= this->transformationSet->getSymmetryFactor();
	}
	else
	{	// old method for handling symmetric binding and unbinding
		if(n_reactants==2)
		{
			//If the binding is symmetric
			if(transformationSet->hasSymBindingTransform()) {
				cout<<endl;
				cout<<"Warning! You have an binding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;

				cout<<endl;
				baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
		if(n_reactants==1)
		{
			if(transformationSet->hasSymUnbindingTransform())
			{
				cout<<endl;
				cout<<"Warning! You have an unbinding rxn (" << name << ") that is symmetric."<<endl;
				cout<<"Make sure that is correct."<<endl;
				cout<<endl;
				baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
				isDimerStyle=true;
			}
		}
	}


	onTheFlyObservables=true;


	// check for population type reactants
	isPopulationType = new bool[n_reactants];
	for( unsigned int i=0; i < n_reactants; ++i )
	{
		isPopulationType[i] = reactantTemplates[i]->getMoleculeType()->isPopulationType();
	}


	// calculate discrete count corrections for symmetric population reactants
	//  e.g. number of reactant pairs = A*(A-1)/2.  Note that the factor of two
	//  is part of the symmetry factor above.
	identicalPopCountCorrection = new int[n_reactants];
	for ( int i=0; i < (int)n_reactants; ++i )
	{
		identicalPopCountCorrection[i] = 0;
		if ( isPopulationType[i] )
		{
			for ( int j=i-1; j >= 0; --j )
			{
				if ( reactantTemplates[i]->getMoleculeType() == reactantTemplates[j]->getMoleculeType() )
				{
					identicalPopCountCorrection[i] = identicalPopCountCorrection[j] + 1;
					break;
				}
			}
		}
	}
}



ReactionClass::~ReactionClass()
{
	delete [] reactantTemplates;
	delete transformationSet;
	for ( unsigned int r = n_reactants; r < n_mappingsets; ++r )
	{
		delete mappingSet[r];
	}

	delete [] mappingSet;
	delete [] isPopulationType;
	delete [] identicalPopCountCorrection;
}


void ReactionClass::setBaseRate(double newBaseRate,string newBaseRateName) {
	if ( this->transformationSet->usingSymmetryFactor() )
	{	this->baseRate = this->transformationSet->getSymmetryFactor() * newBaseRate;   }
	else if (isDimerStyle)
	{	this->baseRate = 0.5 * newBaseRate;   }
	else
	{	this->baseRate = newBaseRate;   }

	this->baseRateParameterName = newBaseRateName;
	update_a();
};


void ReactionClass::resetBaseRateFromSystemParamter() {

	if(!this->baseRateParameterName.empty()) {
		if ( transformationSet->usingSymmetryFactor() ) {
			this->baseRate = transformationSet->getSymmetryFactor() * system->getParameter(this->baseRateParameterName);
		}
		else {
			// TODO: do we need to handle DimerStyle here?? --Justin
			this->baseRate=system->getParameter(this->baseRateParameterName);
		}
		this->update_a();
	}

}



void ReactionClass::printDetails() const {
	cout<< name <<"  (id="<<this->rxnId<<", baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -|"<< this->getReactantCount(r)<<" mappings|\t";
		cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
		//cout<<"head: "<<endl; this->reactantTemplates[r]->printDetails(cout);
		//reactantTemplates[r]->printDetails();
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
	cout<<"\n";
}


void ReactionClass::fire(double random_A_number) {
	//cout<<endl<<">FIRE "<<getName()<<endl;
	fireCounter++;


	// First randomly pick the reactants to fire by selecting the MappingSets
	this->pickMappingSets(random_A_number);


	// Check reactants for correct molecularity:
	if ( ! transformationSet->checkMolecularity(mappingSet) ) {
		// wrong molecularity!  this is a NULL event
		++(System::NULL_EVENT_COUNTER);
		return;
	}


	// output something if the reaction was tagged
	if(tagged) {
		cout<<"#RT "<<this->rxnId<<" "<<this->system->getCurrentTime();
		for(unsigned int k=0; k<n_reactants; k++) {
			cout<<" [";
			for(unsigned int p=0; p<mappingSet[k]->getNumOfMappings();p++) {
				Molecule *mForTag = mappingSet[k]->get(p)->getMolecule();
				cout<<" "<<mForTag->getMoleculeTypeName()<<mForTag->getUniqueID();
			}
			cout<<" ]";
		}
		cout<<endl;
	}


	// Generate the set of possible products that we need to update
	// (excluding new molecules, we'll get those later --Justin)
	this->transformationSet->getListOfProducts(mappingSet,products,traversalLimit);


	// display product molecules for debugging..
	//for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
	//	cout<<">>molecule: "<<(*molIter)->getMoleculeTypeName()<<endl;
	//	(*molIter)->printDetails();
	//	cout<<"<<"<<endl;
	//}


	// Loop through the products (excluding added molecules) and remove from observables
	if (this->onTheFlyObservables) {

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ )
			(*molIter)->removeFromObservables();

		// species observables..
		if(system->getNumOfSpeciesObs()>0) {
			// we can find reactant complexes by following mappingSets to target molecules
			int matches = 0;
			Complex * c;
			for ( unsigned int k=0; k<transformationSet->getNreactants(); k++) {
				// get complexID and check if we've already updated that complex
				int complexId = mappingSet[k]->get(0)->getMolecule()->getComplexID();
				if ( std::find( updatedComplexes.begin(), updatedComplexes.end(), complexId ) == updatedComplexes.end() ) {
					// complex has not been updated, so do it now.
					updatedComplexes.push_back(complexId);
					c = mappingSet[k]->get(0)->getMolecule()->getComplex();
					for(int i=0; i<system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						system->getSpeciesObs(i)->straightSubtract(matches);
					}
				}
			}

			// grab added molecules that are represented as populations and remove from observables
			for ( int k=0; k<transformationSet->getNumOfAddMoleculeTransforms(); k++)
			{
				Molecule * addmol = transformationSet->getPopulationPointer((unsigned int)k);
				if ( addmol == NULL ) continue;

				// get complexID and check if we've already updated that complex
				int complexId = addmol->getComplexID();
				if ( std::find( updatedComplexes.begin(), updatedComplexes.end(), complexId ) == updatedComplexes.end() ) {
					// complex has not been updated, so do it now.
					updatedComplexes.push_back(complexId);
					c = addmol->getComplex();
					for (int i=0; i < system->getNumOfSpeciesObs(); i++) {
						matches = system->getSpeciesObs(i)->isObservable(c);
						system->getSpeciesObs(i)->straightSubtract(matches);
					}
				}
			}
			updatedComplexes.clear();
		}
	}


	// Through the MappingSet, transform all the molecules as neccessary
	//  This will also create new molecules, as required.  As a side effect,
	//  deleted molecules will be removed from observables.
	this->transformationSet->transform(this->mappingSet);


	// Add newly created molecules to the list of products
	this->transformationSet->getListOfAddedMolecules(mappingSet,products,traversalLimit);


	// if complex bookkeeping is on, find all product complexes
	// (this is useful for updating Species Observables and TypeII functions, so keep the info handy).
	// NOTE: this is a brute force approach: check complex of each molecule. there may be a more
	//  elegant way to do this, but it's tricky to get it right.
	if (system->isUsingComplex()) {
		Complex * complex;
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			// get complexID and check if we've already updated that complex
			complex = (*molIter)->getComplex();
			if ( std::find( productComplexes.begin(), productComplexes.end(), complex ) == productComplexes.end() )
				productComplexes.push_back(complex);
		}
	}


	// If we're handling observables on the fly, tell each molecule to add itself to observables.
	if (onTheFlyObservables) {

		// molecule observables..
		for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
			// skip dead molecules
			if ( ! (*molIter)->isAlive() ) continue;
			(*molIter)->addToObservables();
		}

		// species observables..
		if (system->getNumOfSpeciesObs()>0) {
			Complex * c;
			int matches;
			// we can assume that complex bookkeeping is enabled..
			for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter ) {
				// update all species observables for this complex
				c = *complexIter;
				matches = 0;
				for ( int i=0; i < system->getNumOfSpeciesObs(); i++ ) {
					matches = system->getSpeciesObs(i)->isObservable(c);
					system->getSpeciesObs(i)->straightAdd(matches);
				}
			}

			// NOTE: we don't need to handle added population types separately since they are
			//  among the product molecules
		}
	}


	// Now update reaction membership, functions, and update any DOR Groups
	//  also, gather a list of typeII dependencies that will require updating
	typeII_products.clear();
	for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
		Molecule * mol = *molIter;
		MoleculeType * mt = mol->getMoleculeType();

		// If this moleculeType has typeII dependencies, add it to the list
		// (do this for alive and dead molecules, since molecule deletion may influence
		//    the value of a local function)
		if ( mt->getNumOfTypeIIFunctions() > 0 ) {
			if ( std::find( typeII_products.begin(), typeII_products.end(), mt ) == typeII_products.end() )
				typeII_products.push_back( mt );
		}

		//Update this molcule's reaction membership
		//  NOTE: as a side-effect, DORreactions that depend on molecule-scoped local functions
		//   (typeI relationship) will be updated as long as UTL is set appropriately.
		if ( mol->isAlive() )
			mol->updateRxnMembership();
	}


	// update complex-scoped local functions for typeII dependencies
	// NOTE: as a side-effect, dependent DOR reactions (via typeI molecule dependencies) will be updated
	if (system->getEvaluateComplexScopedLocalFunctions()) {
		// for each typeII product molecule, update all dependent local functions
		if (system->isUsingComplex()) {
			// this is the easy way: update all typeI molecules on each complex
			for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
				MoleculeType * mt = *typeII_iter;
				for (int i=0; i < mt->getNumOfTypeIIFunctions(); i++) {
					for ( complexIter = productComplexes.begin(); complexIter != productComplexes.end(); ++complexIter )
						mt->getTypeIILocalFunction(i)->evaluateOn( *complexIter );
				}
			}
		}
		else {
			// this is the hard way: find a representative molecule from each connected set
			//  and evaluate TypeII functions on that representative.
			list <Molecule *> allMols;
			Molecule * mol;
			for ( molIter = products.begin(); molIter != products.end(); molIter++ ) {
				mol = *molIter;
				if ( std::find( allMols.begin(), allMols.end(), mol ) == allMols.end() ) {
					// remember everything connected to this molecule
					//  (so we don't evaluate this connected set multiple times)
					mol->traverseBondedNeighborhood( allMols, ReactionClass::NO_LIMIT );
					// evaluate typeII local functions on this connected set
					for ( typeII_iter = typeII_products.begin(); typeII_iter != typeII_products.end(); ++typeII_iter ) {
						MoleculeType * mt = *typeII_iter;
						for (int i=0; i<mt->getNumOfTypeIIFunctions(); i++)
							mt->getTypeIILocalFunction(i)->evaluateOn( mol, LocalFunction::SPECIES );
					}
				}
			}
		}
	} // done updating complex-scoped local functions


	// display final product molecules for debugging..
	//for( molIter = products.begin(); molIter != products.end(); molIter++ ) {
	//	cout<<">>molecule: "<<(*molIter)->getMoleculeTypeName()<<endl;
	// 	(*molIter)->printDetails();
	//  	cout<<"<<"<<endl;
	//}


	//Tidy up
	products.clear();
	productComplexes.clear();
}











