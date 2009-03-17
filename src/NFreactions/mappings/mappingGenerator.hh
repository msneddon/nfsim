#ifndef MAPPINGGENERATOR_HH_
#define MAPPINGGENERATOR_HH_






	
#include "../NFreactions.hh"

namespace NFcore
{
	

	//!  Knows how to assign mappings in a MappingSet to a particular Molecule
	/*!
	   MapGenerators are simple objects given to TemplateMolecules (so are maintained and
	   deleted by TemplateMolecule objects) that assign Mapping objects to particular
	   Molecule objects.  This allows us to Map specific components of a reactant as
	   we compare the potential molecules to the TemplateMolecule.  (see the TemplateMolecule
	   compare() functions).  If there is a match, the Template Molecule uses this to assign
	   a mapping.  This class is useful because it knows which Mapping to use in the
	   MappingSet provided in the map() function.  This is necessary information as we always
	   reuse Mapping and MappingSet objects for effeciency.
	    @author Michael Sneddon
	 */
	class MapGenerator
	{
		public:
			
			/*!
			   Creates a MapGenerator (generally created in a TransformationSet when
			   we are creating Transformations) that remembers the index of a Mapping
			   in a MappingSet so that we always reuse the same Mapping.
			    @author Michael Sneddon
			 */
			MapGenerator(unsigned int mappingIndex);
			
			/*!
			   Deconstructor that does nothing significant.
			    @author Michael Sneddon
			 */
			~MapGenerator() {};
			
			/*!
			   The one key function of a MapGenerator - this takes a MappingSet and
			   finds the particular Mapping that this MapGenerator needs to create
			   and sets that Mapping to point to the given molecule.  This is called
			   in the TemplateMolecule compare() function.
			    @author Michael Sneddon
			 */
			bool map(MappingSet *mappingSet, Molecule *molecule);
		
		protected:
			/*!
			   Keeps track of the index of the Mapping in the MappingSet that this
			   MapGenerator must assign.
			 */
			unsigned int mappingIndex;
	};

}







#endif /*MAPPINGGENERATOR_HH_*/
