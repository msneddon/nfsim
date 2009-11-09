/*
 * netgen.hh
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#ifndef NETGEN_HH_
#define NETGEN_HH_


#include "../NFcore/NFcore.hh"


using namespace std;

namespace NFcore
{
	// forward declarations for cyclic dependencies
	//class ComplexList;
	//class MappingSetsIterator;
	//class Netgen;
	//class Reaction;
	//class ReactionList;
    //class ReactionClassIterator;


	/* This class implements a unidirectional reaction object.
	 *  (e.g. A + B -> C)  The key members are lists of reactants,
	 *  products and a ratelaw.
	 */
	class Reaction
	{
		public:
			Reaction  ( );
			~Reaction ( );

			vector <Complex *> reactants;
			vector <Complex *> products;

		protected:

		private:

	};


	/*  A container class for Reactions.
	 *  TODO: define methods.
	 */
	class ReactionList
	{
		public:
			ReactionList  ( );
			~ReactionList ( );

			// add reaction to list
			bool addReaction( Reaction * rxn );

		protected:
			vector <Reaction *> reaction_list;

		private:

	};


	/* A container class for Complex objects.
	 *  Key features:
	 *  > Maintains a map from complex labels to complex pointers.
	 *  > Implements a canonical labeling method for complex objects.
	 */
	class ComplexList
	{
		public:
			ComplexList  ( );
			~ComplexList ( );

		protected:
			// Add a complex to the list.
			// returns true if Complex is not already in list, false otherwise.
			bool addComplexToList ( Complex* c );

			// Get pointer to complex by label.
			Complex * getComplexByLabel ( string& label );

			// Get canonical label for a complex.
			string getComplexLabel ( Complex* c );

			// Compare complexes by label.
			// e.g. return label(c1) <=> label(c2)
			// valid return values are among {-1,0,1}
			int compare( Complex* c1, Complex* c2 );

			// Maps complex labels to complex pointers:
			map <string, Complex *>  label_map;
			// A list of complex pointers:
			vector <Complex *>       complex_list;

		private:
			// some iterators that are frequently used.
			map <string, Complex *>::iterator labelMapIter;
			vector <Complex *>::iterator      complexIter;
			vector <Molecule *>::iterator     molIter;
	};


	/* During network generation, each ReactionClass object has an
	 *  associated MappingSetsIterator that efficiently iterates over
	 *  all possible combinations of reactants.  Call method
	 *  'getNextMappingSets' to load a vector with new MappingSets.
	 *  Returns true if a new mappingSets were loaded and false otherwise.
	 */
	class MappingSetsIterator
	{
		public:
			MappingSetsIterator  ( ReactionClass * rc_ );
			~MappingSetsIterator ( );

			bool getNextMappingSets ( vector <MappingSet*> * next_set );

		protected:
			ReactionClass * rc;

		private:

	};


	/* Implements an iterator object for ReactionClasses.
	 *  Call 'getNextReactionClass' to load the next reactionClass
	 *  pointer.  Returns true if a ReactionClass is returned,
	 *  false otherwise.
	 */
	class ReactionClassIterator
	{
		public:
			ReactionClassIterator  ( );
			~ReactionClassIterator ( );

			void setReactionClassList ( vector <ReactionClass *> * rc_list );
			bool getNextReactionClass ( ReactionClass* rc );

		protected:
			vector <ReactionClass *> * reactionClassList;

		private:


	};


	/* The main class for Network Generation.
	 *  This is a wrapper around a system object that performs network
	 *  generation.  Call "network_generation" to populate the reaction
	 *  network.  The class contains members that hold the list of
	 *  complexes and reactions.
	 */
	class Netgen
	{
		public:
			Netgen  (System * sys_ );
			~Netgen ( );

			void generate_network( );

		protected:
			System                * sys;
			ReactionClassIterator rc_iter;
			ComplexList           complex_list;
			ReactionList		  reaction_list;

		private:

	};

}

#endif /* NETGEN_HH_ */
