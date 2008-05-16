#ifndef MAPPINGSET_HH_
#define MAPPINGSET_HH_





#include "../NFreactions.hh"





using namespace std;



namespace NFcore
{

	class Mapping;


	//!  Keeps a list of mappings
	/*!
	    @author Michael Sneddon
	 */
	class MappingSet
	{
		public:
			MappingSet(unsigned int id, vector <Transformation *> &transformations);
			~MappingSet();
			
			
			/*!
			    @author Michael Sneddon
			 */
			bool set(unsigned int mappingIndex, Molecule *m);
			
			
			/*!
			    @author Michael Sneddon
			 */
			Mapping *get(unsigned int mappingIndex) { return mappings[mappingIndex]; }
			
			/*!
			    @author Michael Sneddon
			 */
			void clear() { };
			
			
			/*!
			    @author Michael Sneddon
			 */
			unsigned int getId() const { return id; };
			
		protected:
			unsigned int id;
			
			bool isSet;
			unsigned int n_mappings;
			Mapping ** mappings;
	};


}




#endif /*MAPPINGSET_HH_*/
