#ifndef MAPPINGSET_HH_
#define MAPPINGSET_HH_





#include "../NFreactions.hh"

using namespace std;

namespace NFcore
{
	//Needs a forward declaration here
	class Mapping;
	class Transformation;


	//!  Keeps a list of mappings needed by a particular list of reactants.
	/*!
	 	MappingSet objects manage a list of simple Mappings onto Molecules all
	 	within a single 
	    @author Michael Sneddon
	 */
	class MappingSet
	{
		public:
			
			/*!
			 	Creates a new MappingSet.  This MappingSet constructor will take care
			 	of populating itself with Mapping objects based on the list of 
			 	transformations that are specified in the given vector.  TransformationSet
			 	objects know how to make MappingSet objects which are immediately made
			 	when a ReactantList (or ReactantTree) initializes or expands its
			 	capacity.  They are never destroyed until the simulation is complete.
				@author Michael Sneddon
			*/
			MappingSet(unsigned int id, vector <Transformation *> &transformations);
			
			/*!
			 	Destroys this MappingSet and any Mapping contained by this MappingSet
				@author Michael Sneddon
			*/
			~MappingSet();
			
			
			/*!
			 	Assigns a particular Mapping in this set at the specified index to
			 	a given Molecule.  MapGenerator objects use this function to Map molecules.
			    @author Michael Sneddon
			 */
			bool set(unsigned int mappingIndex, Molecule *m);
			
			
			/*!
			 	Retrieves a particular Mapping at the specified index.  This is used
			 	when we need to actually transform a given Molecule involved with a rule
			 	that fired.
			    @author Michael Sneddon
			 */
			Mapping *get(unsigned int mappingIndex) { return mappings[mappingIndex]; }
			
			/*!
			 	This clears the Mapping objects belonging to this MappingSet.  It is generally
			 	good practice to do this once this MappingSet is no longer being used, but this
			 	function is not required for any real functionality.  Note that currently this
			 	function does nothing (for effeciency).  To make it do something, look at the
			 	file mappingSet.cpp and use the (slightly slower) yet less error prone version
			 	of these functions.
			    @author Michael Sneddon
			 */
			void clear() { };
			
			
			/*!
			 	Returns the Id of this MappingSet.  The Id of the MappingSet is used to track the
			 	MappingSet as it exists in a ReactantList or ReactantTree.  This Id should be unique
			 	within the list or tree, but is not unique within the system.
			    @author Michael Sneddon
			 */
			unsigned int getId() const { return id; };
			
		protected:
			
			/*!
				The Id of this MappingSet assigned by the ReactantList or ReactantTree in which
				this MappingSet lives.
			*/
			unsigned int id;
			
			/*!
				Keeps track of whether or not this MapppingSet is being used.  This is merely for
				error checking and serves no other useful purpose.
			*/
			bool isSet;
			
			/*!
				Keeps track of the number of Mapping objects contained in this MappingSet.
			*/
			unsigned int n_mappings;
			
			/*!
				An array of pointers to Mapping objects.  This is where the actual Mappings are stored.
			*/
			Mapping ** mappings;
	};

}




#endif /*MAPPINGSET_HH_*/
