








#ifndef MAPPING_HH_
#define MAPPING_HH_

#include "../NFreactions.hh"

namespace NFcore
{

//!  Keeps a pointer to a molecule and remembers a single component to act on
	/*!
	    @author Michael Sneddon
	 */
	class Mapping
	{
		public:
			
			/*!
			    @author Michael Sneddon
			 */
			Mapping(unsigned int type, unsigned int index);
			
			/*!
			    @author Michael Sneddon
			 */
			~Mapping();
			
			/*!
			    @author Michael Sneddon
			 */
			unsigned int getType() const;
			
			/*!
			    @author Michael Sneddon
			 */
			unsigned int getIndex() const;
			
			/*!
			    @author Michael Sneddon
			 */
			Molecule * getMolecule() const;
			
			/*!
			    @author Michael Sneddon
			 */
			void clear();
			
			/*!
			    @author Michael Sneddon
			 */
			bool setMolecule(Molecule *m);

		protected:
			unsigned int type;
			unsigned int index;
			Molecule * m;
	};

}









#endif /*MAPPING_HH_*/
