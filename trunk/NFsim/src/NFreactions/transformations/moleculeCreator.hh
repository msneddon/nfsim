/*
 * moleculeCreator.hh
 *
 *  Created on: Mar 4, 2011
 *      Author: justin
 */

#ifndef MOLECULECREATOR_HH_
#define MOLECULECREATOR_HH_


#include "../NFreactions.hh"

using namespace std;

namespace NFcore
{

	class MoleculeType;
	class TemplateMolecule;

	class MoleculeCreator
	{
		public:
			MoleculeCreator( TemplateMolecule *          _template_molecule,
					         MoleculeType *              _molecule_type,
			                 vector < pair<int,int> > &  _component_states    );

			~MoleculeCreator( );

			// create a molecule (don't get pointer)
			void create( );

			// create a molecule and get pointer
			Molecule * create_molecule( );

			// get the template molecule
			TemplateMolecule * getTemplateMolecule () const { return template_molecule; };

			// check if this is a population type
			bool isPopulationType () const { return population_type; }

			// get pointer to population molecule (NULL if molecule is a particle)
			Molecule * get_population_pointer() const;

		protected:

			// for population types:
			//  try to find map to existing population molecule
			Molecule * map_molecule();

			// pointer to newly created molecule object
			Molecule *         molecule_object;

			// is this a population or a regular molecule type?
			bool               population_type;

			// template that matches created molecule
			TemplateMolecule * template_molecule;

			// type of molecule to create
			MoleculeType *      molecule_type;

			// non-default component index / state index pairs
			vector < pair<int,int> >  component_states;

		private:
			vector < pair<int,int> >::iterator   comp_iter;
	};

}
#endif /* MOLECULECREATOR_HH_ */
