

#ifndef REACTANTTREE_HH_
#define REACTANTTREE_HH_


#include "../NFreactions.hh"

namespace NFcore
{
	//Forward Declarations
	class TransformationSet;
	class MappingSet;
	class ReactantContainer;


	//!  Maintains a tree of MappingSets needed by Distribution of Rates Reactions
	/*!
	 *  This is one of the more complex classes in NFsim.  It is written in order to handle
	 *  distribution of rates reactions.  These types of reactions have many reactants, but
	 *  each reactant can participate with a different rate.  The rate is determined by the
	 *  set of molecules that the reactant is connected to.  Therefore, we need a tree (or we could have also
	 *  used logarithmic classes) to efficiently select the next reactant weighted by
	 *  each of its propensities.
	 *
	 *  A key example of this is in the chemotaxis system.  Clusters of receptors influence
	 *  the rate of CheA autophosphorylation.
	 *
	 *
	 */
	class ReactantTree : public ReactantContainer {


		public:

			/*!
				Create a reactant tree, based on the given transformation set (so we can create mappings
				at the given reactant index.  The tree will be initialized such that it can immediately
				hold the number of objects specified by the init_capacity.
			*/
			ReactantTree(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);


			/*!
				Standard deconstructor for the tree
			*/
			virtual ~ReactantTree();

			/*!
				Returns the number of mappingSets that have been added to this tree
			*/
			virtual int size() const { return n_mappingSets; };

			/*!
				Adds a new MappingSet to this tree and returns a pointer to the new mapping set for you
				to map (usually by comparing to some template molecule).  Warning: even if you don't use
				this mapping set, it will be counted until you pop it or remove it! (see removeMappingSet()).  This has
				a special behavior in the ReactantTree.  This method merely gets you a pointer to the next
				available mappingset: if you decide to keep the mapping, you have to confirm it. Confirming
				the mapping actually places it in the tree.  If you don't do this, it won't ever be selected
				by the tree.
			*/
			virtual MappingSet * pushNextAvailableMappingSet();


			/*!
				Once you have pushed on the next available mappingSet, you have to confirm that push.  In
				the event that mappingsets are cloned, this will confirm all mappingsets that have been
				cloned off of the given mapping set.
			*/
			void confirmPush(int mappingSetId, double rateFactor);



			/*!
				Removes the very last mappingSet that was added to the list.  You should actually now
				use the removeMappingSet function instead, just to be careful.
			*/
			virtual void popLastMappingSet();

			/*!
				Removes the mapping set with the specified mappingSetId.  Be careful here: make sure the mapping
				set is actually on the list before trying to remove or else this will give you an error!
			*/
			virtual void removeMappingSet(unsigned int mappingSetId);


			/*!
				This is the key utility of a reactant tree: given a random value, we can select the next reactant
				to fire with this function that weights the selection of any particular reactant by
			*/
			void pickReactantFromValue(MappingSet *&ms, double value, double baseRate);

			/*!
				When a local function value changes, it must update the value in the reactant tree.  This
				method allows you to update values without changing the mappingSet membership of this tree.
			 */
			void updateValue(unsigned int mappingSetId, double newRateFactor);

			/*!
				Returns a MappingSet so that a DOR can evaluate a local function on it.
			 */
			virtual MappingSet * getMappingSet(unsigned int mappingSetId) const;



			/*!
				Outputs basic details about this tree - used only for debugging, because
				it really does print out everything
			*/
			virtual void printDetails() const;




			/*!
				Returns the combined rate factor sum of this tree, which is needed by
				the DOR reactionclass in order to properly update its propensity
			*/
			double getRateFactorSum() const { return leftRateFactorSum[0]; };


			/*!
				Returns the depth of the tree (which you shouldn't ever really need...)
			*/
			int getDepthOfTree() const { return treeDepth; };


		protected:

			/*!
				Removes a mapping set from the tree only.  This method should only be
				used by the tree!
			*/
			void removeFromTreeOnly(int msTreeArrayPosition, unsigned int mappingSetId);

			/*!
				If we try to add more than this tree can handle, we have to expand it.  Because
				the tree is a complete binary tree, this will necessarily double its capacity,
				or enlarge it more depending on the given new capacity.
			*/
			void expandTree(int newCapacity);


			TransformationSet *ts;       //Keeps track of the set of transformations
			unsigned int reactantIndex;  //the index of the tree

			//Basic tree parameters and constants
			int maxElementCount;
			int treeDepth;
			int numOfNodes;

			//The tree stored as 3 arrays indexed as:
			// Children are at 2i and 2i+1
			// Parent is at i/2 (using integer division)
			// root is at index 1
			// index 0 is empty and always zero
			double * leftRateFactorSum;
			int * leftElementCount;
			int * rightElementCount;

			//The actual list of mappingSets, stored as a list
			MappingSet ** mappingSets;



			/////  below are three arrays that are necessary to quickly grab a particular
			////   mappingSet from the tree or vice versa

			//Given a mappingSet Id, this tells us what position it is in
			//in the mappingSets array
			int * msPositionMap;

			//Given a mappingSet Id, this tells us what position it is in
			//in the actual tree
			int * msTreePositionMap;

			//Given a tree position index, this tells us the mappingSet Id
			//which we can use to get our mappingSet out of the mappingSets array
			int * reverseMsTreePositionMap;


			//The number of mappingSets currently set
			int n_mappingSets;

			//The index of the first molecule in the tree
			unsigned int firstMappingTreeIndex;
	};
}







#endif /* REACTANTTREE_HH_ */
