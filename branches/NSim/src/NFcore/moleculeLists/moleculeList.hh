#ifndef MOLECULELISTS_HH_
#define MOLECULELISTS_HH_


#include "../NFcore.hh"

namespace NFcore
{
	class MoleculeType;
	class Molecule;

	//!  Keeps a set of molecules neatly for a MoleculeType
	/*!
	   MoleculeTypes have to store a link to all Molecule objects that currently
	   exist.  This is done using this list.  This list sets aside a certain number
	   of pointers for molecules of a particular type so that adding and deleting
	   molecules from the system is very fast.  The code enforces that all existing
	   Molecules are maintained with this list again for maximum effeciency and memory
	   management.
	   @author Michael Sneddon
	*/
	class MoleculeList
	{
		public:

			/*!
				Create a new MoleculeList that stores the given type of molecule with
				an initial capacity (that is dynamically expanded if you add more than
				the given capacity of molecules).
			*/
			MoleculeList(MoleculeType *mt, int init_capacity, int finalCapacity);

			/*!
				Deletes the list and all Molecule objects associated with this list.  A
				Molecule cannot exist outside of a list, so deleting this deletes all
				member molecules.
			*/
			~MoleculeList();

			/*!
				Returns the number of molecules that are on the list
			*/
			int size() const;

			/*!
				Returns a pointer to the Molecule with the given index.
			*/
			NFcore::Molecule *at(int index) const;

			/*!
				Allocates a spot for a new Molecule in the system.  The given
				pointer is assigned to the new Molecule (so it should be passed
				in as a null pointer) and the method returns the index this
				Molecule was assigned.
			*/
			int create(Molecule *&m);

			/*!
				Removes a Molecule from this list with the given ID.
			*/
			void remove(int listId, Molecule *m);

			/*!
				Removes the very last Molecule on the list (which is the one you just
				created if you haven't touched anything)
			*/
			void removeLast();

			/*!
				Print out some (maybe too much) diagnostic debug information about the list.
			*/
			void printDetails();

			static const int NO_LIMIT = -1;

		protected:

			/*! The number of Molecule objects currently on the list */
			int n_molecules;

			/*! The number of Molecule objects currently on the list */
			int lastAllocated;

			/*! The maximum number of Molecules that can be added before the list resizes itself */
			int capacity;

			/*! The maximum number of Molecules that can be added, period. */
			int finalCapacity;

			/*! Keeps track of the type of molecule stored by this list */
			MoleculeType *mt;

			/*! The actual array of Molecules that are stored */
			Molecule **mArray;

			/*! Allows the list to map index values of Molecules to the index values in the list array  */
			int *molPos;
	};


}





#endif /*MOLECULELISTS_HH_*/
