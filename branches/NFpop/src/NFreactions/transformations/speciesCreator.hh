#ifndef SPECIESCREATOR_HH_
#define SPECIESCREATOR_HH_


#include "../NFreactions.hh"

using namespace std;

namespace NFcore
{
	class MoleculeType;
	class TemplateMolecule;

	class SpeciesCreator
	{
		public:
			SpeciesCreator(
					vector <MoleculeType *> &productMoleculeTypes,
					vector < vector <int> > &stateInformation,
					vector < vector <int> > &bindingSiteInformation);
			SpeciesCreator(vector <TemplateMolecule *> &templates);
			~SpeciesCreator();
			
			void create();
		
		protected:
			
			unsigned int n_molecules;
			MoleculeType ** moleculeTypes;
			
			Molecule **newMoleculeCreations;
			

			//Save the configuration of states we have to set
			//nd prefix stands for non-default 
			unsigned int n_ndStates;
			int *ndStateMolecule;
			int *ndStateIndex;
			int *ndStateValue;
			
			
			//Save the bonds that we have to create
			unsigned int n_bonds;
			int *bMolecule1;
			int *bMolecule2;
			int *bSite1;
			int *bSite2;
	};




}
























#endif /*SPECIESCREATOR_HH_*/
