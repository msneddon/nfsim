#ifndef MOLECULELISTS_HH_
#define MOLECULELISTS_HH_


#include "../NFcore.hh"

namespace NFcore
{
	class MoleculeType;


	class MoleculeList 
	{
		
		
		public:
			MoleculeList(MoleculeType *mt, int init_capacity);
			~MoleculeList();
		
			int size() const;
			Molecule *at(int index) const;
			
			int create(Molecule *&m);
			void remove(int listId);
			void removeLast();
			
			void printDetails();
		
		protected:
			int n_molecules;
			int lastAllocated;
			int capacity;
			
			MoleculeType *mt;
			Molecule **mArray;
			int *molPos;
		
	};


}





#endif /*MOLECULELISTS_HH_*/
