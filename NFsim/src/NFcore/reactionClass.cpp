


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


	//cout<<"good, very good."<<endl;
	//exit(0);



	//Check here to see if we have molecule types that are the same across different reactants
	//Because if so, we will give a warning
	if(n_reactants>2) cerr<<"Warning!! You created a reaction ("<< name <<") that has more than 2 reactants.  This has not been extensively tested!"<<endl;
	if(n_reactants==2)
	{
		//If the reactants are of the same type, then we have to make a few special considerations
		if(reactantTemplates[0]->getMoleculeType()==reactantTemplates[1]->getMoleculeType())
		{
			cout<<endl;
			cout<<"Warning! You have a binding rxn (" << name << ") that allows a moleculeType to bind another of the same type."<<endl;
			cout<<"Make sure that is correct, because this can potentially make long polymers or large aggregates."<<endl;
			cout<<endl;
		}
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
	if(n_reactants==1) {
		if(transformationSet->hasSymUnbindingTransform()) {
			cout<<endl;
			cout<<"Warning! You have an unbinding rxn (" << name << ") that is symmetric."<<endl;
			cout<<"Make sure that is correct."<<endl;
			cout<<endl;
			baseRate = baseRate*0.5;  //We have to correct the rate to get the proper factor
			isDimerStyle=true;
		}
	}
	onTheFlyObservables=true;
}


ReactionClass::~ReactionClass()
{
	delete [] reactantTemplates;
	delete transformationSet;
	delete [] mappingSet;

}



void ReactionClass::resetBaseRateFromSystemParamter() {

	if(!this->baseRateParameterName.empty()) {
		this->baseRate=system->getParameter(this->baseRateParameterName);
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


void ReactionClass::fire(double random_A_number)
{
	fireCounter++; //Remember that we fired

	//if(getName()=="Rule2") {
	//	this->printFullDetails();
	//	this->system->printAllObservableCounts(0);
	//	this->system->printAllReactions();
	//	cout<<"---"<<endl; cout<<"found: "<<products.size()<<" products."<<endl;
	//	exit(1);
	//}

	//cout<<"\n\n-----------------------\nfiring: "<<name<<endl;;
	//this->system->printAllObservableCounts(0);
	//this->system->printAllReactions();


	// First randomly pick the reactants to fire by selecting the MappingSets
	this->pickMappingSets(random_A_number);


	// Check reactants for correct molecularity
	// (This is more general than the previous checks for intra- and inter-complex
	//   binding, and is compliant with the BNGL definition of molecularity.
	//   However, we still aren't checking for correct product molecularity,
	//   which is a harder problem.)
	if ( !(transformationSet->checkMolecularity( mappingSet )) )
	{
		// wrong molecularity!  this is a NULL event
		++(System::NULL_EVENT_COUNTER);
		return;
	}



	// output something if the reaction was tagged
	if(tagged)
	{
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


	//	if(getName()=="Rule2") {
	//		cout<<"---"<<endl;cout<<"-------found: "<<products.size()<<" products."<<endl;
	//		for( molIter = products.begin(); molIter != products.end(); molIter++ )
	//			(*molIter)->printDetails();
	//		cout<<system->getObservableByName("Lig_free")->getCount()<<"/"<<system->getObservableByName("Lig_tot")->getCount()<<endl;
	//	}

	//cout<<"found: "<<products.size()<<" products."<<endl;


	// Loop through the products and remove them from their observables
	if(this->onTheFlyObservables)
	{
		// remove products from molecule observables
		for ( molIter = products.begin(); molIter != products.end(); molIter++ )
		{
			//cout<<"Removing: "<<(*molIter)->getMoleculeTypeName()<<"_"<<(*molIter)->getUniqueID()<<endl;
			//(*molIter)->printDetails();
			(*molIter)->removeFromObservables();
		}

		// Remove complexes containing products from species oservables
		if(system->getNumOfSpeciesObs()>0) {
			bool found = false;
			// we can assume that complex bookkeeping is on, and that each reactant
			// is in a separate (and single) complex
			int matches = 0;
			for (int k=0; k<transformationSet->getNumOfReactants(); k++)
			{
				// get complexID and check if we've already updated that complex
				int complexId = mappingSet[k]->get(0)->getMolecule()->getComplexID();
				found = false;
				for(unsigned int k2=0; k2<updatedComplexes.size(); k2++) {
					if(updatedComplexes.at(k2)==complexId) {
						found = true;
						break;
					}
				}
				// if we already handled this, go to the next product
				if(found) continue;

				// if we didn't handle this, remember that we're handling it now..
				updatedComplexes.push_back(complexId);

				// update species observables for this complex
				Complex *c = mappingSet[k]->get(0)->getMolecule()->getComplex();
				for(int i=0; i<system->getNumOfSpeciesObs(); i++) {
					matches = system->getSpeciesObs(i)->isObservable(c);
					for(int j=0; j<matches; j++) system->getSpeciesObs(i)->straightSubtract();
				}
			}
			updatedComplexes.clear();
		}

	}


	// Through the MappingSet, transform all the molecules as neccessary
	//  This will also create new molecules, as required
	this->transformationSet->transform(this->mappingSet);


	// Add newly created molecules to the list of products
	this->transformationSet->getListOfAddedMolecules(mappingSet,products,traversalLimit);


	// If we're handling observables on the fly, tell each molecule to add
	//  itself back into observable counts.
	//  NOTE: this will also take care of adding new molecules into the observables. --Justin, 8Mar2011
	if(onTheFlyObservables)
	{
		// add product molecules to molecule observables
		for( molIter = products.begin(); molIter != products.end(); molIter++ )
		{
			// skip dead molecules
			if ( !(*molIter)->isAlive() ) continue;
			(*molIter)->addToObservables();
		}

		// add complexes containing products into species observables
		if(system->getNumOfSpeciesObs()>0)
		{
			// we can assume that complex bookkeeping is on, and that each reactant
			bool found = false;
			for( molIter = products.begin(); molIter != products.end(); molIter++ )
			{
				// skip dead molecules
				if ( !(*molIter)->isAlive() ) continue;

				// get complexID and check if we've already updated that complex
				int complexId = (*molIter)->getComplexID();
				found = false;
				for (unsigned int k=0; k<updatedComplexes.size(); k++) {
					if(updatedComplexes.at(k)==complexId) {
						found = true;
						break;
					}
				}
				// if we already handled this, go to the next product
				if(found) continue;

				// if we didn't handle this, remember that we're handling it now..
				updatedComplexes.push_back(complexId);

				// update species observables for this complex
				Complex *c = (*molIter)->getComplex();
				int matches = 0;
				for ( int i=0; i<system->getNumOfSpeciesObs(); i++ ) {
					matches = system->getSpeciesObs(i)->isObservable(c);
					for(int j=0; j<matches; j++) system->getSpeciesObs(i)->straightAdd();
				}
			}
			updatedComplexes.clear();
		}
	}

	//if(getName()=="Rule2") {
	//	cout<<"---"<<endl;cout<<"-------found: "<<products.size()<<" products."<<endl;
	//	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	//				(*molIter)->printDetails();
	//
	//	cout<<system->getObservableByName("Lig_free")->getCount()<<"/"<<system->getObservableByName("Lig_tot")->getCount()<<endl;
	//}


	// Now update reaction membership, functions, and update any DOR Groups
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	{
		// skip over dead molecules (NOTE: we do need this now.  --Justin, 8Mar2011
		if ( !(*molIter)->isAlive() ) continue;
	  	(*molIter)->updateRxnMembership();
	  	(*molIter)->updateTypeIIFunctions();
	  	(*molIter)->updateDORRxnValues();
	  	//(*molIter)->printDetails();
	}


	//Molecule::printMoleculeList(products);
	//this->printFullDetails();
	//this->system->printAllObservableCounts(0);
	//this->system->printAllReactions();
	//exit(1);



	//cout<<",  everything done"<<endl;
	//Tidy up
	products.clear();
	
}











