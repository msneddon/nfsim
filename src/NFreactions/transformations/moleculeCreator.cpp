/*
 * moleculeCreator.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: justin
 */

#include "moleculeCreator.hh"

using namespace NFcore;

MoleculeCreator::MoleculeCreator(
		TemplateMolecule           * _template_molecule,
		MoleculeType               * _molecule_type,
		vector < pair <int,int> >  & _component_states   )
{
	template_molecule = _template_molecule;
	molecule_type     = _molecule_type;
	component_states  = _component_states;

	// is this a population type molecule?
	population_type   = molecule_type->isPopulationType();

	if ( population_type )
	{
		// Handle a population molecule...
		// try to find this molecule among those initiatize in seed species block:
		molecule_object = map_molecule();
	}

}



MoleculeCreator::~MoleculeCreator()
{

}


void
MoleculeCreator::create()
{

	if ( isPopulationType() )
	{
		// increment population
		molecule_object->incrementPopulation();
	}
	else
	{
		//Create the molecule
		molecule_object = molecule_type->genDefaultMolecule();

		//Set the component state values correctly
		for( comp_iter = component_states.begin();  comp_iter != component_states.end();  ++comp_iter )
		{
			molecule_object->setComponentState( (*comp_iter).first, (*comp_iter).second );
		}

		//Prep the molecule and enterinto the simulation
		molecule_type->addMoleculeToRunningSystemButDontUpdate( molecule_object );
	}

}


Molecule *
MoleculeCreator::create_molecule()
{

	if ( isPopulationType() )
	{
		// increment population
		molecule_object->incrementPopulation();
	}
	else
	{
		//Create the molecule
		molecule_object = molecule_type->genDefaultMolecule();

		//Set the component state values correctly
		for( comp_iter = component_states.begin();  comp_iter != component_states.end();  ++comp_iter )
		{
			molecule_object->setComponentState( (*comp_iter).first, (*comp_iter).second );
		}

		//Prep the molecule and enterinto the simulation
		molecule_type->addMoleculeToRunningSystemButDontUpdate( molecule_object );
	}

	// return a pointer to the new molecule, yeah!
	return molecule_object;
}

Molecule *
MoleculeCreator::get_population_pointer() const
{
	if ( population_type )
	{
		return molecule_object;
	}
	else
	{
		return NULL;
	}
}

Molecule *
MoleculeCreator::map_molecule()
{
	Molecule * mol;
	Molecule * found_mol = NULL;
	int        im, n_mol;

	n_mol = molecule_type->getMoleculeCount();
	for ( im = 0; im < n_mol; ++im )
	{
		mol = molecule_type->getMolecule( im );
		if ( mol->isAlive() )
		{
			if ( template_molecule->compare( mol ) )
			{
				if ( found_mol != NULL )
				{
					cerr << "!! While initializing MoleculeCreator: found multiple populations that match template!\n"
					     << "Not sure what to do... Quitting." << endl;
					exit(1);

				}
				// save this molecule pointer
				found_mol = mol;
			}
		}
	}

	if ( found_mol == NULL )
	{
		cerr << "!! While initializing MoleculeCreator: could not find molecule matching template!\n"
			 << "(HINT: Try adding the molecule to seed species block in model file.)\n"
		     << "Not sure what to do... Quitting." << endl;
		exit(1);
	}

	return found_mol;
}
