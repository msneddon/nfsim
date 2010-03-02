/*
 * netgen.hh
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#ifndef NETGEN_HH_
#define NETGEN_HH_


#include "../NFcore/NFcore.hh"
#include "../NFreactions/reactantLists/reactantList.hh"
#include "../NFreactions/reactions/reaction.hh"


using namespace std;

namespace NFcore
{
	// forward declarations for cyclic dependencies
	//class ComplexList;
	//class MatchSetIter;
	//class Netgen;
	//class Reaction;
	//class ReactionList;
    //class ReactionClassIter;


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
	 *  associated MatchSetIter that efficiently iterates over
	 *  all possible combinations of reactant matches.  Call method
	 *  'nextMatchSet' to load a vector with new MappingSets.
	 *  Returns true if a new mappingSets were loaded and false otherwise.
	 */
	class MatchSetIter
	{
		public:
			MatchSetIter  ( ReactionClass * _rc );
			~MatchSetIter ( );

			void reset ( );
			void update ( );
			bool nextMatchSet ( vector <MappingSet*> & match_set );

			BasicRxnClass            * rc;
			unsigned int             n_reactants;

		protected:

			vector <unsigned int>    curr_set;
			vector <ReactantList *>  reactantLists;
			bool                     more_sets;

			// advance iterator
			void advance ();


		private:
			 vector <unsigned int>::iterator     curr_set_iter;
			 vector <ReactantList *>::iterator   reactantList_iter;
			 vector <MappingSet *>::iterator     mappingSet_iter;

	};


	/* Implements an iterator object for ReactionClasses.
	 *  Call 'nextReactionClass' to load the MatchSetIter
	 *  corresponding to the next reactionClass. Returns true if a
	 *  a new MatchSetIter is returned, false otherwise.
	 */
	class ReactionClassIter
	{
		public:
			ReactionClassIter  ( );
			~ReactionClassIter ( );

			void setSystem ( vector <ReactionClass *> * _reactionClasses );
			void reset ( );
			MatchSetIter * nextReactionClass ( );

		protected:
			vector <ReactionClass*>            * reactionClasses;
			vector <MatchSetIter*>               matchSetIters;
			vector <MatchSetIter*>::iterator     matchSetIters_iter;

		private:
		    vector<ReactionClass*>::iterator     rc_iter;
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

			System                * sys;

			ComplexList           complex_list;
			ReactionList		  reaction_list;

			ReactionClassIter     rc_iter;
			MatchSetIter          * match_set_iter;

			vector <MappingSet *>            match_set;
			vector <MappingSet *>::iterator  match_iter;

		protected:
			vector <Molecule *>              product_molecules;

		private:
			 // containers for copying species targets
		     list <Molecule *>                   mol_list;
			 map <Molecule *, Molecule *>        mol_map;

			 // handy iterators for copying species
			 map <Molecule *, Molecule *>::iterator  mol_map_iter;
			 list <Molecule *>::iterator             mol_iter;
	};

}

#endif /* NETGEN_HH_ */
