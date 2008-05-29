#ifndef SPECIES_HH_
#define SPECIES_HH_





namespace NFcore {

	class SpeciesList {
		
		
		create();
		remove();
		at();
		size();
		
		
		vector <Species *> species;
	};




	class Species {
	
		public:
		
			Species(Molecule *m);
			~Species();
			
			bool compare(Species *s);
			
			void add();
			void subtract();
			
	
		protected:
	
			int n_molecules;
			MoleculeType ** moleculeTypes;
			
			int n_states;
			int **states;
			
			int n_bonds;
			int **bonds;
			
			
			//To keep track of state index
			static const int M_INDEX=0;
			static const int STATE_INDEX=1;
			static const int STATE_VALUE=2;
			
			
			//To keep all the numbers straight when we
			//keep track of data in the bonds array
			static const int M_INDEX_1=0;
			static const int M_INDEX_2=1;
			static const int B_SITE_INDEX_1=2;
			static const int B_SITE_INDEX_2=3;
			
			
			
			
			
			
			
			
			static int SpeciesCount=0;
	};




};























#endif /*SPECIES_HH_*/
