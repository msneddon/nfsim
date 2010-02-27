/*
 * netgen.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: justin
 */

#include "netgen.hh"

using namespace std;
using namespace NFcore;


// Constructor
Netgen::Netgen ( System * _sys )
{
	sys = _sys;
	rc_iter.setSystem( &(sys->allReactions) );
}

// Destructor
Netgen::~Netgen ( )
{

}


// Network generation method
void Netgen::generate_network ( )
{
	cout << ">>> Generate network!" << endl;

	// network generation iteration loop
	//while ( true )
	//{
		// ReactionClass loop
	    size_t ii = 0;
		while ( match_set_iter = rc_iter.nextReactionClass() )
		{
			cout << "reactionClass " << ii << endl;  ++ii;

			// MatchSet loop
			size_t jj = 0;
			while ( match_set_iter->nextMatchSet( match_set ) )
			{
				cout << "matchSet " << jj << endl;  ++jj;

				// loop over MappingSets
				for( match_iter = match_set.begin(); match_iter < match_set.end(); match_iter++ )
				{
					// get first Map in the set
					Mapping * map = ((*match_iter)->mappings)[0];
					// get the target molecule reference
					Molecule * mol = map->m;
					// get the complex ID for the target molecule
					Complex * complex = mol->getComplex();

					vector <Molecule *> molvec_orig = new vector <Molecule *>;
					vector <Molecule *> molvec_copy = new vector <Molecule *>;
					map <Molecule *, Molecule *> mol_map = new map <Molecule *, Molecule *>;

					// get molecules in complex
					for ( c->molIter = (c->complexMembers).begin;
						    c->molIter < (c->complexMembers).end;
					           (c->molIter)++  )
					{
						// get molecule pointer
						Molecule * mol_orig = *(c->molIter);
						// copy molecule pointer to "original" vector
						molvec_orig->push_back( mol );

						// get molecule type
						MoleculeType * mol_type = mol->getMoleculeType();
						// make copy of molecule
						Molecule * mol_copy = mol_type->genDefaultMolecule();
						// add copy to vector
						molvec_copy->push_back( mol_copy );

						// add map from original to copy
						mol_map.insert ( pair<Molecule *,Molecule *>(mol_orig,mol_copy) );

						// copy component states
						// TODO: what about equivalent classes of components?
						size_t cIndex;
						for ( cIndex = 0; cIndex < mol_type->getNumOfComponents(); cIndex++ )
						{
							mol_copy->setComponentState( cIndex, mol_orig->getComponentState(cIndex) );
						}
					}

					// now that all molecules in complex are copied, add bonds
					vector <Molecule *>::iterator mol_iter = new vector <Molecule *>::iterator;
					for ( mol_iter = molvec_orig->begin(); mol_iter < molvec_orig->end(); mol_iter++ )
					{
				        Molecule * mol = *mol_iter;
						MoleculeType * mol_type = mol->getMoleculeType();

						// loop over components to find bonds
						size_t cIndex;
				        for ( cIndex = 0; cIndex < mol_type->getNumOfComponents(); cIndex++ )
						{
				        	if ( mol->isBindingSiteOpen(cIndex) ) {  continue;  }

				        	Molecule * bond_mol = mol->getBondedMolecule( cIndex );
				        	int bond_index = mol->getComponentIndexOfBond( cIndex );

							//bind(Molecule *m1, int cIndex1, Molecule *m2, int cIndex2);
							//unbind(Molecule *m1, int bSiteIndex);
						}
					}

					// re-target MappingSet to the copy

				}

				// load match_set into the reactionClass
				match_set_iter->rc->set_match( match_set );

				// apply rule to match_set, and get product molecules
				molecule_vec.clear();
				match_set_iter->rc->apply( molecule_vec );


				// create copy of targets
				// Complex * Molecule::getComplex()
				// re-map the mappingSets onto the copy (how?)

				// fire reaction rule

				// add products to ComplexList (how?),
				//     delete products from instantiated moleculeLists if complex already was found

				// construct reaction

				// add reaction to ReactionList

			}
		}

	    //rc_iter.reset();
	//}
}
