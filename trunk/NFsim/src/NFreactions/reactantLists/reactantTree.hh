#ifndef REACTANTTREE_HH_
#define REACTANTTREE_HH_



#include "../NFreactions.hh"


namespace NFcore
{
	class ReactantTree {
		
		
		public:
			
			ReactantTree(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);
			~ReactantTree();
			
			int size() const { return n_mappingSets; };
			double getRateFactorSum() const { return leftRateFactorSum[0]; };
			
			
			
			
			/*
			 * Order of operations in DOR Reactions
			 * 1) use pushNextAvailableMappingSet to get a MappingSet, with the rateFactor from the reactant molecule
			 * 2) use compare to see if we should keep it
			 * 3) if(!match) call removeMappingSet
			 * 
			 * 
			 */
			
			
			MappingSet * pushNextAvailableMappingSet(double rateFactor);
			
			
			void removeMappingSet(unsigned int mappingSetId);
			void updateValue(unsigned int mappingSetId, double newRateFactor);
			void pickReactantFromValue(MappingSet *&ms, double value, double baseRate) const;

			
			
			
			int getDepthOfTree() const { return treeDepth; };
			void printDetails() const;
		
		
		protected:
		
			//Basic tree parameters and constants
			unsigned int capacity;
			unsigned int treeDepth;
			unsigned int treeFirstMoleculeIndex;
			unsigned int numOfNodes;
			
			//The tree stored as 3 arrays indexed as:
			// Children are at 2i and 2i+1
			// Parent is at i/2 (using integer division)
			// root is at index 1
			// index 0 is empty and always zero
			double * leftRateFactorSum;
			unsigned int * leftElementCount;
			unsigned int * rightElementCount;
			
			//The actual list of MappingSets
			MappingSet ** mappingSets;
		
			
			
			
			
			//The number of molecules currently in the tree
			int n_mappingSets;
			
			/*! The transformation set of the ReactionClass that owns this list */
			TransformationSet *ts;
			
			/*! The index of the reactant that this list maintains */
			unsigned int reactantIndex;
			
			//The index of the first molecule in the tree
			unsigned int firstMoleculeTreeIndex;
		
	};
}


#endif /*REACTANTTREE_HH_*/
