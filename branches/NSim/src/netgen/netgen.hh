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
