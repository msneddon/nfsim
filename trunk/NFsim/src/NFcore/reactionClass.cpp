


#include "NFcore.hh"


using namespace std;
using namespace NFcore;



ReactionClass::ReactionClass(string name, double baseRate, TransformationSet *transformationSet)
{
	//cout<<"\n\ncreating reaction "<<name<<endl;

	isDimerStyle=false;
	//Setup the basic properties of this reactionClass
	this->name = name;
	this->baseRate = baseRate;
	this->fireCounter = 0;
	this->a = 0;
	this->traversalLimit = ReactionClass::NO_LIMIT;
	this->transformationSet = transformationSet;

	//Set up the template molecules from the transformationSet
	this->n_reactants = transformationSet->getNreactants();
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
	mappingSet = new MappingSet *[n_reactants];





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



void ReactionClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
		//reactantTemplates[r]->printDetails();
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}


void ReactionClass::fire(double random_A_number)
{
	fireCounter++; //Remember that we fired
	//cout<<"\n\n-----------------------\nfiring: "<<name<<endl;;

	//First randomly pick the reactants to fire by selecting the MappingSets
	pickMappingSets(random_A_number);

	//Generate the set of possible products that we need to update
	this->transformationSet->getListOfProducts(mappingSet,products,traversalLimit);

	//cout<<"found: "<<products.size()<<" products."<<endl;
	//Loop through the products and remove them from thier observables
	//cout<<"------------------------------------------"<<endl;
	if(this->onTheFlyObservables) {
		for( molIter = products.begin(); molIter != products.end(); molIter++ )
		{
			//cout<<"Removing: "<<(*molIter)->getMoleculeTypeName()<<"_"<<(*molIter)->getUniqueID()<<endl;
			//(*molIter)->printDetails();
			(*molIter)->removeFromObservables();
		}
	}

	//cout<<",  obs removed";

	//cout<<"before: "<<endl;
	//Molecule::printMoleculeList(products);

	//Through the MappingSet, transform all the molecules as neccessary
	this->transformationSet->transform(this->mappingSet);

//	cout<<",  transformed updated";

	//Tell each molecule in the list of products to add itself back into
	//the counts of observables and update its class lists, and update any DOR Groups
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	{
		if(onTheFlyObservables) (*molIter)->addToObservables();
	}
//	cout<<",  obs updated";

//	cout<<"after:"<<endl;
	for( molIter = products.begin(); molIter != products.end(); molIter++ )
	{
	  	(*molIter)->updateRxnMembership();
	  	(*molIter)->updateTypeIIFunctions();
	  	(*molIter)->updateDORRxnValues();
	 // 	(*molIter)->printDetails();
	}
	//Molecule::printMoleculeList(products);

	//	cout<<",  everything done"<<endl;
	//Tidy up
	products.clear();
}











