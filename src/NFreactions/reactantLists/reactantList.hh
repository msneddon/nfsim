#ifndef REACTANTLIST_HH_
#define REACTANTLIST_HH_




#include "../NFreactions.hh"

namespace NFcore
{
	//Forward Declarations
	class TransformationSet;
	class MappingSet;


	//!  Maintains a list of MappingSets needed by ReactionClass
	/*!
	  This is essentially a specialized vector implementation that allows a ReactionClass
	  to easily keep track of all the MappingSet objects that can be involved in the reaction.
	  All MappingSet objects stored in this class are created once and reused throughout the
	  course of the simulation.  This prevents new MappingSet objects to be created and destroyed.
	  This class has the ability to automatically expand its capacity when extra mappings are
	  needed, so there is no need for the user to manage these details.  This class allows for
	  constant time removal, insertion, and random selection of MappingSets, while std::vector
	  requires linear time removal.  We gain this speedup because the ordering of the list
	  is unimportant.  To use this class, call the pushNextAvailableMappingSet to get a pointer
	  to the available MappingSet.  Pass this MappingSet to a TemplateMolecule to actually
	  map it onto molecules.  If you don't end up using this MappingSet, call popLastMappingSet
	  immediately.  To remove MappingSets, call removeMappingSet with the ID of the MappingSet
	  you want to remove.  You can get this Id from the MappingSet object and Molecule objects
	  also keep a vector of the MappingSet objects that point to it.
	    @author Michael Sneddon
	 */
	class ReactantList : public ReactantContainer
	{
		public:

			/*!
				Creates a new empty ReactantList with the given initial capacity (default is 50).  This capacity
				should roughly be set to the number of mappings you expect this list to have.  A reactantList must
				also be told what its reactionIndex is in the reaction and the TransformationSet of the
				reaction so that it can populate itself with MappingSets.
			 */
			ReactantList(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);


			/*!
				Deconstructor that deletes all mappingSets associated with this list.
			 */
			virtual ~ReactantList();


			/*!
				Returns the number of mappingSets that have been added to this list
			 */
			virtual int size() const { return n_mappingSets; };

			/*!
				Returns the sum population of all mappingSets that have been added to this list
			 */
			virtual int getPopulation() const;


			/*!
				Adds a new MappingSet to this list and returns a pointer to the new mapping set for you
				to map (usually by comparing to some template molecule).  Warning: even if you don't use
				this mapping set, it will be counted until you pop it! (see popLastMappingSet()).
			 */
			virtual MappingSet * pushNextAvailableMappingSet();

			/*!
				Removes the very last mappingSet that was added to the list.
			 */
			virtual void popLastMappingSet();

			/*!
				Removes the mapping set with the specified mappingSetId.  Be careful here: make sure the mapping
				set is actually on the list before trying to remove or else this will give you an error!
			 */
			virtual void removeMappingSet(unsigned int mappingSetId);

			/*!
				Randomly selects a MappingSet from the list of available MappingSets.
			 */
			void pickRandom(MappingSet *&ms);

			/*!
				Randomly selects a MappingSet from the population weighted list of available MappingSets.
			 */
			void pickRandomFromPopulation(MappingSet *&ms);


			virtual MappingSet * getMappingSet(unsigned int mappingSetId) const;


			/*!
				Outputs basic details about this list - used only for debugging.
			 */
			virtual void printDetails() const;

		protected:

			/*! Maintains the number of mappingSets on this list */
			int n_mappingSets;

			/*! The total capacity that this list can hold */
			int capacity;

			/*! The transformation set of the ReactionClass that owns this list */
			TransformationSet *ts;

			/*! The index of the reactant that this list maintains */
			unsigned int reactantIndex;

			/*! The array that maps MappingSet Ids onto the location in the list
			    that the MappingSet is sitting in. */
			unsigned int * msPositionMap;

			/*! The actual array that stores a list of pointers to MappingSet objects */
			MappingSet **mappingSets;
	};
}
















#endif /*REACTANTLIST_HH_*/
