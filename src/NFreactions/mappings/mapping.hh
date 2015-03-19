




#ifndef MAPPING_HH_
#define MAPPING_HH_

#include "../NFreactions.hh"

namespace NFcore
{
	//!  Keeps a pointer to a molecule and remembers a single component to act on
	/*!
	 	This is the low level unit that is used by reactions to identify the particular
	 	component of a molecule to update or change when a transformation is
	 	required.  It is maintained by a MappingSet which is just a collection of
	 	Mapping objects.  The "type" of mapping is identified in the same way as the
	 	"type" of transformation.  See the Transformation class for a list of those
	 	types.  The "index" simply is the index of the binding site or state of a molecule.
	 	It is possible that some "types" do not require an "index" (as in the case of
	 	adding or deleting an entire molecule).
	    @author Michael Sneddon
	 */
	class Mapping
	{
		public:

			/*!
			 	Creates a new Mapping object with the specified "type" (as defined by the
			 	types in the Transformation class) and a specified "index" which is the
			 	binding site or state index (in a MoleculeType) that this Mapping will map
			 	onto.  The type and index does not (and should not) change during a
			 	simulation.
			    @author Michael Sneddon
			 */
			Mapping(unsigned int type, int index);

			/*!
			 	Destroys this mapping, but does not delete the molecule it is mapped to.
			    @author Michael Sneddon
			 */
			~Mapping();

			/*!
			 	Returns the type of mapping this is.  See Transformation for the possible
			 	types that can be assigned.
			    @author Michael Sneddon
			 */
			unsigned int getType() const;

			/*!
			 	Returns the index into either the binding site or state of the Molecule
			 	this Mapping points to.
			    @author Michael Sneddon
			 */
			int getIndex() const;

			/*!
			 	Once this Mapping has been mapped to a molecule, you can get a pointer to that
			 	Molecule with this function.
			    @author Michael Sneddon
			 */
			Molecule * getMolecule() const;

			/*!
			 	This clears the Mapping by forgetting the Molecule that this maps to.  It is good
			 	practice to call this function after the Mapping becomes unused (as is done by
			 	ReactantList) but it is not always necessary as long as you have another way to
			 	keep track of mappings (as ReactionList does).
			    @author Michael Sneddon
			 */
			void clear();

			/*!
			 	This sets the Mapping to actually point to a particular Molecule in the system.  This
			 	Molecule is remembered until a new Molecule is set or the clear() function is called.
			    @author Michael Sneddon
			 */
			bool setMolecule(Molecule *m);


			static void clone(Mapping *original, Mapping *newClone);


			void printDetails() const;
			void printDetails(ostream &o) const;


		protected:

			/*!
				The "type" of mapping this is, taken from the same list of types in
				the Transformation class.
			*/
			unsigned int type;

			/*!
				The index of the binding site or state that this Mapping keeps track of.
			*/
			int index;

			/*!
				A pointer to the Molecule that this Mapping keeps track of.  It can only
				keep track of one Molecule at a time.
			*/
			Molecule * m;
	};



	inline
	Molecule * NFcore::Mapping::getMolecule() const
	{
		//Make sure the Molecule points somewhere.  For effeciency, this check
		//can be removed.
	//	if(m==NULL) {
	//		cout<<"Trying to get a molecule from a null mapping (in class Mapping)!! Quitting!"<<endl;
	//		exit(1);
	//	}
		return m;
	}
}









#endif /*MAPPING_HH_*/
