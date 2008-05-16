#ifndef TRANSFORMATION_HH_
#define TRANSFORMATION_HH_









#include "../NFreactions.hh"

namespace NFcore
{

	
	class Transformation
	{
		public:
			~Transformation();
			unsigned int getType() const { return type; };
			unsigned int getStateOrSiteIndex() const { return stateORsiteIndex; };
			int getNewStateValue() const { return newStateValue; };
			unsigned int getPartnerReactantIndex() const { return otherReactantIndex; };
			unsigned int getPartnerMappingIndex() const { return otherMappingIndex; };
			
			static Transformation * genStateChangeTransform(unsigned int stateIndex, int newStateValue);
			static Transformation * genBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex);
			static Transformation * genBindingTransform2(unsigned int bSiteIndex);
			static Transformation * genUnbindingTransform(unsigned int bSiteIndex);
			static Transformation * genAddMoleculeTransform();
			static Transformation * genRemoveMoleculeTransform();
			static Transformation * genEmptyTransform();
			
			//Keep track of types
			static const unsigned int STATE_CHANGE = 0;
			static const unsigned int BINDING = 1;
			static const unsigned int UNBINDING = 2;
			static const unsigned int REMOVE = 3;
			static const unsigned int ADD = 4;
			static const unsigned int SKIP = 5; //used to mark a transformation placeholder that we don't actually want to execute, as
						          
			
		protected:
			Transformation();
			
			unsigned int type;
			unsigned int stateORsiteIndex;
			
			//For a state change transformation
			int newStateValue;
			
			//For a binding reaction
			unsigned int otherReactantIndex;
			unsigned int otherMappingIndex;
			
			//For a creation of a new molecule
			//Species s
	};
}







#endif /*TRANSFORMATION_HH_*/
