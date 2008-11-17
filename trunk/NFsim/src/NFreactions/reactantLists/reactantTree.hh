

#ifndef REACTANTTREE_HH_
#define REACTANTTREE_HH_


#include "../NFreactions.hh"

namespace NFcore
{
	//Forward Declarations
	class TransformationSet;
	class MappingSet;


	//!  Maintains a tree of MappingSets needed by Distribution of Rates Reactions
	/*!
	 *  This is one of the more complex classes in NFsim.  It is written in order to handle
	 *  distribution of rates reactions.  These types of reactions have many reactants, but
	 *  each reactant can participate with a different rate.  The rate is determined by the
	 *  module that the reactant is in.  Therefore, we need a tree (or we could have also
	 *  used logarithmic classes) to efficiently select the next reactant weighted by
	 *  each of its propensities.
	 *
	 *  A key example of this is in the chemotaxis system.  Clusters of receptors influence
	 *  the rate of CheA autophosphorylation.
	 *
	 *
	 */
	class ReactantTree {


		public:
			ReactantTree(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);
			~ReactantTree();

			/*!
				Returns the number of mappingSets that have been added to this tree
			*/
			int size() const { return n_mappingSets; };

			/*!
				Adds a new MappingSet to this tree and returns a pointer to the new mapping set for you
				to map (usually by comparing to some template molecule).  Warning: even if you don't use
				this mapping set, it will be counted until you pop it! (see popLastMappingSet()).  This has
				a special behavior in the ReactantTree.  This method merely gets you a pointer to the next
				available mappingset: if you decide to keep the mapping, you have to confirm it. Confirming
				the mapping actually places it in the tree.  If you don't do this, it won't ever be selected
				by the tree.
			*/
			MappingSet * pushNextAvailableMappingSet();


			void confirmPush(int mappingSetId, double rateFactor);



			/*!
				Removes the very last mappingSet that was added to the list.
			*/
			void popLastMappingSet();

			/*!
				Removes the mapping set with the specified mappingSetId.  Be careful here: make sure the mapping
				set is actually on the list before trying to remove or else this will give you an error!
			*/
			void removeMappingSet(unsigned int mappingSetId);

			/*!
				Randomly selects a MappingSet from the list of available MappingSets.
			*/
			void pickRandom(MappingSet *&ms);

			/*!
				This is the key utility of a reactant tree: given a random value, we can select the next reactant
				to fire with this function that weights the selection of any particular reactant by
			*/
			void getReactantFromValue(MappingSet *&ms, double value, double baseRate);

			/*!

			 */
			void updateValue(unsigned int mappingSetId, double newRateFactor);



			/*!
				Outputs basic details about this list - used only for debugging.
			*/
			void printDetails();




			void expandTree(int newCapacity);



			//Once you insert into this tree, you will get the position you were
			//inserted into.  This will never change until you are removed
			//int insert(Molecule * m, double rateFactor);
			//void remove(Molecule * m, unsigned int rxnListIndex);



			//void updateValue(Molecule * m, unsigned int rxnListIndex, double newRateFactor);
			//Molecule * getReactantFromValue(double value, double baseRate) const;


			//int getNumOfMolecules() const { return numOfMolecules; };
			double getRateFactorSum() const { return leftRateFactorSum[0]; };
			int getDepthOfTree() const { return treeDepth; };
			//void printDetails() const;


		protected:

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
			//which we can use to get our mappingSet out of the mappingSets
			//array
			int * reverseMsTreePositionMap;


			//bool hasOpenPush;


			//The number of mappingSets currently set
			int n_mappingSets;

			//The index of the first molecule in the tree
			unsigned int firstMappingTreeIndex;
	};
}







#endif /* REACTANTTREE_HH_ */
