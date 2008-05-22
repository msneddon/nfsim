//
///*
//
//class NFlist 
//{
//	public:
//
//			/*!
//				Creates a new empty ReactantList with the given initial capacity (default is 50).  This capacity
//				should roughly be set to the number of mappings you expect this list to have.  A reactantList must
//				also be told what its reactionIndex is in the reaction and the TransformationSet of the
//				reaction so that it can populate itself with MappingSets.
//			 */
//			ReactantList(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);
//			
//			
//			/*!
//				Deconstructor that deletes all mappingSets associated with this list.
//			 */
//			~ReactantList();
//		
//		
//			/*!
//				Returns the number of mappingSets that have been added to this list
//			 */
//			int size() const;
//			
//			/*!
//				Adds a new MappingSet to this list and returns a pointer to the new mapping set for you
//				to map (usually by comparing to some template molecule).  Warning: even if you don't use
//				this mapping set, it will be counted until you pop it! (see popLastMappingSet()).
//			 */
//			MappingSet * pushNextAvailableMappingSet();
//			
//			/*!
//				Removes the very last mappingSet that was added to the list.
//			 */
//			void popLastMappingSet();
//			
//			/*!
//				Removes the mapping set with the specified mappingSetId.  Be careful here: make sure the mapping
//				set is actually on the list before trying to remove or else this will give you an error!
//			 */
//			void removeMappingSet(unsigned int mappingSetId);
//			
//			/*!
//				Randomly selects a MappingSet from the list of available MappingSets.
//			 */
//			void pickRandom(MappingSet *&ms);
//			
//			
//			/*!
//				Outputs basic details about this list - used only for debugging.
//			 */
//			void printDetails();
//		
//		protected:
//			
//			/*! Maintains the number of mappingSets on this list */
//			int n_mappingSets;
//			
//			/*! The total capacity that this list can hold */
//			int capacity;
//			
//			/*! The transformation set of the ReactionClass that owns this list */
//			TransformationSet *ts;
//			
//			/*! The index of the reactant that this list maintains */
//			unsigned int reactantIndex;
//			
//			/*! The array that maps MappingSet Ids onto the location in the list
//			    that the MappingSet is sitting in. */
//			unsigned int * msPositionMap;
//			
//			/*! The actual array that stores a list of pointers to MappingSet objects */
//			MappingSet **mappingSets;
//	};
//	
//	
	
	