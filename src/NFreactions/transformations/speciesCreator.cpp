

#include "speciesCreator.hh"

using namespace NFcore;


SpeciesCreator::SpeciesCreator(
		vector <MoleculeType *> &productMoleculeTypes,
		vector < vector <int> > &stateInformation,
		vector < vector <int> > &bindingSiteInformation)
{
	try {
		this->n_molecules = productMoleculeTypes.size();
		cout << "\t\t\tNew species Creator!! num of molecules in species: "<<n_molecules<<endl;

		//Copy over the molecules
		n_molecules = productMoleculeTypes.size();
		moleculeTypes = new MoleculeType *[n_molecules];
		for(unsigned int m=0; m<n_molecules; m++)
		{
			moleculeTypes[m]=productMoleculeTypes.at(m);
			cout << "\t\t\t\t" << m << ": " << moleculeTypes[m]->getName() << endl;
		}

		newMoleculeCreations = new Molecule *[n_molecules];

		//Save the configuration of states we have to set
		//nd prefix stands for non-default
		n_ndStates = stateInformation.at(0).size();
		ndStateMolecule = new int [n_ndStates];
		ndStateIndex = new int [n_ndStates];
		ndStateValue = new int [n_ndStates];

		for(unsigned int s=0; s<n_ndStates; s++)
		{
			ndStateMolecule[s]=stateInformation.at(0).at(s);
			ndStateIndex[s]=stateInformation.at(1).at(s);
			ndStateValue[s]=stateInformation.at(2).at(s);

			cout<<"state"<<s<<" : "<<ndStateMolecule[s]<<" "<< ndStateIndex[s]<<" "<<ndStateValue[s]<<endl;
		}

		//Save the bonds that we have to create
		n_bonds = bindingSiteInformation.at(0).size();
		bMolecule1 = new int [n_bonds];
		bMolecule2 = new int [n_bonds];
		bSite1 = new int [n_bonds];
		bSite2 = new int [n_bonds];

		for(unsigned int b=0; b<n_bonds; b++)
		{
			bMolecule1[b]=bindingSiteInformation.at(0).at(b);
			bMolecule2[b]=bindingSiteInformation.at(2).at(b);
			bSite1[b]=bindingSiteInformation.at(1).at(b);
			bSite2[b]=bindingSiteInformation.at(3).at(b);

			cout<<"bonds"<<b<<" : "<<bMolecule1[b]<<","<<bSite1[b]<<"  to  "<<bMolecule2[b]<<","<<bSite2[b]<<endl;


		}

	} catch(std::exception& err){
		cout<<"Error when creating a new SpeciesCreator object: SpeciesCreator was not properly creaated... quitting."<<endl;
		err.what();
		exit(1);
	}
}

SpeciesCreator::SpeciesCreator(vector <TemplateMolecule *> &templates)
{

	cerr<<"Calling an old and nonfunctional SpeciesCreator constructor.  Quitting."<<endl;
	exit(1);

//	try {
//
//		//This method of initializing is more labor intensive at this end, but makes more sense if you
//		//are initializing a model directly from code
//		this->n_molecules = templates.size();
//
//
//		//First, some vectors to dynamically keep track of the states and bonds to create
//		//This would normally just be input from the xml, but in this case, we need to generate this
//		//information from the set of templates
//		const int proMolTypeIndex = 0;
//		const int stateIndex = 1, stateValue=2;
//		const int bSiteIndex = 1, partnerMolTypeIndex = 2, partnerBsiteIndex = 3;
//		vector <MoleculeType *> productMoleculeTypes;
//		vector < vector <int> > stateInformation;
//		vector < vector <int> > bindingSiteInformation;
//
//		//First, set up the stateInformation vector
//		{
//			//v1 maps index into productMoleculetypes, v2 maps the stateIndex
//			//v3 maps the state value
//			vector <int> v1,v2,v3;
//			stateInformation.push_back(v1);
//			stateInformation.push_back(v2);
//			stateInformation.push_back(v3);
//		}
//
//		//And here we go with the bindingSiteInformation vector
//		{
//			//v1 maps index into productMoleculetypes, v2 maps the bindingSiteIndex
//			//v3 maps the partner's index into productMoleculeTypes, and v4 maps the
//			vector <int> v1,v2,v3, v4;
//			bindingSiteInformation.push_back(v1);
//			bindingSiteInformation.push_back(v2);
//			bindingSiteInformation.push_back(v3);
//			bindingSiteInformation.push_back(v4);
//		}
//
//
//
//		//@TODO Fix this species creator!
//		//Loop through each template and remember the states and binding sites we have to set
//		for(unsigned int m=0; m<n_molecules; m++) {
//
////			TemplateMolecule *tm = templates.at(m);
////			for(unsigned int s=0; s<tm->get(); s++) {
////				stateInformation.at(proMolTypeIndex).push_back(m);
////				stateInformation.at(stateIndex).push_back(tm->getStateIndex(s));
////				stateInformation.at(stateValue).push_back(tm->getStateValue(s));
////			}
//
////			for(unsigned int b=0; b<tm->getNumBindingSites(); b++) {
////
////				//Find the molecule we are bound to
////				TemplateMolecule *tm2 = tm->getBondedTemplateMolecule(b);
////				unsigned int m2 = 0;
////				for(m2=m+1; m2<n_molecules; m2++) {
////					if(tm2==tm->getBondedTemplateMolecule(m2)) {
////						//Remember the information, only if we found the other partner
////						//This ensures we only add this bond once
////						bindingSiteInformation.at(proMolTypeIndex).push_back(m);
////						bindingSiteInformation.at(bSiteIndex).push_back(tm->getBindingSiteIndex(b));
////						bindingSiteInformation.at(partnerMolTypeIndex).push_back(m2);
////						bindingSiteInformation.at(partnerBsiteIndex).push_back(tm->getBindingSiteIndexOfBondedTemplate(b));
////					}
////				}
////			}
//		}
//
//
//		//Now, we just have to copy the information in the vectors over to arrays to keep our memory effeciency and speed good
//		//when we run the simulation
//
//		//Copy over the molecules
//		moleculeTypes = new MoleculeType *[n_molecules];
//		for(unsigned int m=0; m<n_molecules; m++)
//		{
//			moleculeTypes[m]=templates.at(m)->getMoleculeType();
//		}
//
//		newMoleculeCreations = new Molecule *[n_molecules];
//
//		//Save the configuration of states we have to set
//		//nd prefix stands for non-default
//		n_ndStates = stateInformation.at(0).size();
//		ndStateMolecule = new int [n_ndStates];
//		ndStateIndex = new int [n_ndStates];
//		ndStateValue = new int [n_ndStates];
//
//		for(unsigned int s=0; s<n_ndStates; s++)
//		{
//			ndStateMolecule[s]=stateInformation.at(0).at(s);
//			ndStateIndex[s]=stateInformation.at(1).at(s);
//			ndStateValue[s]=stateInformation.at(2).at(s);
//		}
//
//		//Save the bonds that we have to create
//		n_bonds = bindingSiteInformation.at(0).size();
//		bMolecule1 = new int [n_bonds];
//		bMolecule2 = new int [n_bonds];
//		bSite1 = new int [n_bonds];
//		bSite2 = new int [n_bonds];
//
//		for(unsigned int b=0; b<n_bonds; b++)
//		{
//			bMolecule1[b]=bindingSiteInformation.at(0).at(b);
//			bMolecule2[b]=bindingSiteInformation.at(1).at(b);
//			bSite1[b]=bindingSiteInformation.at(2).at(b);
//			bSite2[b]=bindingSiteInformation.at(3).at(b);
//		}
//
//	} catch(std::exception& err){
//		cout<<"Error when creating a new SpeciesCreator object: SpeciesCreator was not properly creaated... quitting."<<endl;
//		err.what();
//		exit(1);
//	}
}



SpeciesCreator::~SpeciesCreator()
{
	delete [] bMolecule1;
	delete [] bMolecule2;
	delete [] bSite1;
	delete [] bSite2;

	delete [] ndStateMolecule;
	delete [] ndStateIndex;
	delete [] ndStateValue;

	delete [] newMoleculeCreations;
	delete [] moleculeTypes;


}

void SpeciesCreator::create()
{
	//Create the molecules in this new species
	for(unsigned int m=0; m<n_molecules; m++)
	{
		newMoleculeCreations[m] = moleculeTypes[m]->genDefaultMolecule();
	}

	//Set the component state values correctly
	for(int c=0; c<n_ndStates; c++)
	{
		newMoleculeCreations[ndStateMolecule[c]]->setComponentState(ndStateIndex[c],ndStateValue[c]);
	}

	//Make the bonds
	for(unsigned int b=0; b<n_bonds; b++)
	{
		Molecule *m1 = newMoleculeCreations[bMolecule1[b]];
		Molecule *m2 = newMoleculeCreations[bMolecule2[b]];
		Molecule::bind(m1,bSite1[b],m2,bSite2[b]);
	}

	//Create all the molecules
	//Molecule *m = mt->genDefaultMolecule();

	//Set all the states etc...

	//Prep the molecules and enter them into the simulation
	for(unsigned int m=0; m<n_molecules; m++)
	{
		moleculeTypes[m]->addMoleculeToRunningSystem(newMoleculeCreations[m]);
		//cout<<"Created:"<<endl;
		//newMoleculeCreations[m]->printDetails();
		//cout<<endl;
	}
}
