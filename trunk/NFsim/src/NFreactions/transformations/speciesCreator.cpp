

#include "speciesCreator.hh"

using namespace NFcore;


SpeciesCreator::SpeciesCreator(
		vector <MoleculeType *> &productMoleculeTypes,
		vector < vector <int> > &stateInformation,
		vector < vector <int> > &bindingSiteInformation)
{
//	try {
	this->n_molecules = productMoleculeTypes.size();
	cout<<"New species Creator!! num of molecules in species: "<<n_molecules<<endl;
	
	//Copy over the molecules
	n_molecules = productMoleculeTypes.size();
	moleculeTypes = new MoleculeType *[n_molecules];
	for(unsigned int m=0; m<n_molecules; m++)
	{
		moleculeTypes[m]=productMoleculeTypes.at(m);
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
		bMolecule2[b]=bindingSiteInformation.at(1).at(b);
		bSite1[b]=bindingSiteInformation.at(2).at(b);
		bSite2[b]=bindingSiteInformation.at(3).at(b);
	} 
	
	
//	} catch(std::exception& err){
//		cout<<"in here..."<<endl;
//		err.what();
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
	
	//Set the state values correctly
	for(unsigned int s=0; s<n_ndStates; s++)
	{
		newMoleculeCreations[ndStateMolecule[s]]->setState(ndStateIndex[s],ndStateValue[s]);
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
	}
}