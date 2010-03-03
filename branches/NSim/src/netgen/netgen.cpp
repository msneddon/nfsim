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
	// continue iterations until no new species are generated (or max_iters is exceeded)
	//while ( new_species )
	//{
		// ReactionClass loop:
	    // rc_iter is an object that iterates over a list of reactionClasses,
        //  but the return value is not the reactionClass itself, but rather an object
	    //  that iterates over reactant combinations, i.e. match set iterator.
	    size_t iReactionClass = 0;
		while ( (match_set_iter = rc_iter.nextReactionClass()) )
		{
			cout << "reactionClass " << iReactionClass << endl;

			// get reactionClass ptr
			ReactionClass * rc = match_set_iter->rc;

			// get transformation set for this rule
			TransformationSet * ts = rc->transformationSet;

			// MatchSet loop
			// call "nextMatchSet" method from the match_set iterator to get the next
			// set of reactants.  The reactants, or more specifically the MappingSets from
			// the reactant patterns into the matching species, are loaded into a vector
			// called "match_set"
			size_t iMatchSet = 0;
			//match_set_iter->reset();
			while ( match_set_iter->nextMatchSet( match_set ) )
			{
				cout << "matchSet " << iMatchSet << endl;

				// Now that we have a set of Maps into reactants, we need to identify the
				// corresponding complexes and create copies.  New maps are created that
				// pointed into the copied complexes (as induced by the original maps).
				// The rule transformation is applied to the copy, so the original complexes
				// are preserved. The products of the rule are captured and analyzed to determine
				// the specific reaction (i.e. specific reactants -> specific products). If novel
				// complexes are created, they are added to the ComplexList


				// iterate over the matching reactants in match_set:
				// Copy the complex and create a new MappingSet that is pointed to the copy.
				size_t iReactant = 0;
				for( match_iter = match_set.begin(); match_iter < match_set.end(); match_iter++ )
				{
					cout << "reactant " << iReactant << endl;

					// get mappingSet ptr
					MappingSet * mapset = *match_iter;

					// create blank mappingSets (to hold new map)
					MappingSet * new_mapset = ts->generateBlankMappingSet( iReactant, 0 );

					// replace the original map on match_set vector with the new map
					*match_iter = new_mapset;

					// find the complex that the mapset is pointing to:
					Mapping * map_orig = (mapset->mappings)[0];     // get first Map in the set
					Molecule * mol_orig = map_orig->m;              // get the target molecule reference
					Complex * complex = mol_orig->getComplex();     // get the complex ID for the target molecule


					// clear out containers
					mol_list.clear();        // mol_list is a list of pointers to molecules in the complex.
					mol_map.clear();         // mol_map is a map from the original molecules to the copies.

					// get molecules in complex and make copies
					list <Molecule *> & complex_members = complex->complexMembers;
					for ( mol_iter = complex_members.begin();
						    mol_iter != complex_members.end();
					           mol_iter++  )
					{
						// get molecule pointer
						Molecule * mol_orig = *mol_iter;
						//debug
						mol_orig->printDetails();

						// put molecule pointer on "original" vector
						mol_list.push_back( mol_orig );

						// get molecule type, and ask the type object to make a new molecule
						MoleculeType * mol_type = mol_orig->getMoleculeType();
						Molecule * mol_copy = mol_type->genDefaultMolecule();

						// add map from original to copy
						mol_map.insert ( pair<Molecule *,Molecule *>(mol_orig,mol_copy) );

						// copy component states
						// TODO: what about equivalent classes of components?
						int iComponent;
						for ( iComponent = 0; iComponent < mol_type->getNumOfComponents(); iComponent++ )
						{
							mol_copy->setComponentState( iComponent, mol_orig->getComponentState(iComponent) );
						}
					}

					// now that all molecules in complex are copied, add bonds
					for ( mol_iter = mol_list.begin(); mol_iter != mol_list.end(); mol_iter++ )
					{
				        Molecule * mol_orig = *mol_iter;
						MoleculeType * mol_type = mol_orig->getMoleculeType();

						// get pointer to copy of this molecule
						Molecule * mol_copy = (mol_map.find( mol_orig ))->second;

						// loop over components to find bonds
						int iComponent;
				        for ( iComponent = 0; iComponent < mol_type->getNumOfComponents(); iComponent++ )
						{
				        	// see if this component is bonded
				        	if ( mol_orig->isBindingSiteOpen(iComponent) ) {  continue;  }

				        	// if so, get pointer to binding partner
				        	Molecule * bond_mol_orig = mol_orig->getBondedMolecule( iComponent );
				        	// get pointer to copy of binding partner
				        	Molecule * bond_mol_copy = (mol_map.find( bond_mol_orig ))->second;

				        	// get component index of the bond on the binding partner
				        	int bond_index = mol_orig->getComponentIndexOfBond( iComponent );

				            // check if copy is already bonded
				        	if ( mol_copy->isBindingSiteBonded(iComponent) ) {  continue;  }

				        	// add bond between copied molecules
				        	Molecule::bind( mol_copy, iComponent, bond_mol_copy, bond_index );
						}
					}

					// point new MappingSet to the copy
					// TODO: handle cloned mappingSets
					for ( size_t iMap = 0; iMap < mapset->n_mappings; iMap++ )
					{

						Mapping * map     = (mapset->mappings)[iMap];
						Mapping * new_map = (new_mapset->mappings)[iMap];

						// get molecule that map points to
						Molecule * mol_orig_ptr = map->getMolecule();
						// get pointer to the copy of molecule
						Molecule * mol_copy_ptr = (mol_map.find( mol_orig_ptr ))->second;

                        // point the new map to the copy molecule and component
						new_map->setMolecule( mol_copy_ptr );
					}

					++iReactant;
				}

				// load match_set into the reactionClass
				// recall that match_set now holds Mappings to the copied reactants
				match_set_iter->rc->set_match( match_set );

				// apply rule to match_set, and get product molecules
				product_molecules.clear();
				match_set_iter->rc->apply( product_molecules );

				cout << "product molecules: " << endl;
                vector <Molecule *>::iterator molvec_iter;
                for ( molvec_iter = product_molecules.begin(); molvec_iter < product_molecules.end(); molvec_iter++ )
                {
                	Molecule * mol = *molvec_iter;
                	mol->printDetails();
                }


				// add products to ComplexList,
				//     delete products from instantiated moleculeLists if complex already was found

				// construct reaction

				// add reaction to ReactionList

                // delete MappingSet objects in the match_set vector (no longer needed since the products
                //  no longer match the reactant patterns, in general).
                for( match_iter = match_set.begin();
						match_iter < match_set.end(); match_iter++ )
				{   delete *match_iter;   }

				++iMatchSet;
			}

			++iReactionClass;
		}

		// track new species created and formally add them to the sysem (i.e. set alive)

		// reset reactionClass iterator to the first reaction
	    //rc_iter.reset();
	//}
}
