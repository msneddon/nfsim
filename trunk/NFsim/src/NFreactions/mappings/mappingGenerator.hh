#ifndef MAPPINGGENERATOR_HH_
#define MAPPINGGENERATOR_HH_






	
#include "../NFreactions.hh"

namespace NFcore
{
	

	//!  Knows how to assign mappings in a MappingSet to a particular Molecule
	/*!
	    @author Michael Sneddon
	 */
	class MapGenerator
	{
		public:
			MapGenerator(unsigned int mappingIndex);
			~MapGenerator();
			bool map(MappingSet *mappingSet, Molecule *molecule);
			
		protected:
			unsigned int mappingIndex;
	};

}







#endif /*MAPPINGGENERATOR_HH_*/
