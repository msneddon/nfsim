#ifndef TRANSFORMATION_HH_
#define TRANSFORMATION_HH_




#include "../mapping/mapping.hh"
#include "../../NFcore/NFcore.hh"

#include <vector>




namespace NFcore {

	//Forward declarations
	class TemplateMapping;
	class MappingSet;
	class Mapping;
	class Transformation;
	class System;
	class MoleculeType;
	class Molecule;
	class TemplateMolecule; 
	class ReactionClass;
	
	
	

	class Transformation {
		
		public:
			virtual ~Transformation();
			
			
			virtual void transform(Mapping **m) {};
			virtual void printDetails();
			
			int getIndex() const {return t_index;};
			int getNumOfMappings() const {return n_mappings;};
			int getType() const {return type;};
			const char * getTypeName() const;
			
			
			
			
			
			const static unsigned int STATE_CHANGE = 0;
			const static unsigned int BINDING = 1;
			const static unsigned int UNBINDING = 2;
			
			const static unsigned int INCREMENT_STATE = 3;
			const static unsigned int DECREMENT_STATE = 4;
			
			
			static Transformation * genStateChangeTransform(TemplateMolecule *tm, const char *stateName, int stateValue, ReactionClass *r);
			static Transformation * genBindingTransform(TemplateMolecule *tm1, TemplateMolecule *tm2, const char *bSiteName1, const char *bSiteName2, ReactionClass *r);
			static Transformation * genUnbindingTransform(TemplateMolecule *tm, const char *bSiteName, ReactionClass *r);
			
			
		protected:
			Transformation(ReactionClass *r, unsigned int type, unsigned int n_mappings);
			
			
			
			TemplateMapping ** templateMappings;
			
			unsigned int t_index;
			unsigned int n_mappings;
			unsigned int type;
			
		private:
			static unsigned int currentIndex;
	};

	
	
	class StateChangeTransform : public Transformation {
		
		public:
			StateChangeTransform(TemplateMapping *tm, int newStateValue, ReactionClass *r);
			virtual ~StateChangeTransform() {};
			
			virtual void transform(Mapping **m);

		protected:
			int newStateValue;
		
	};
	
	
	
	class BindingTransform : public Transformation {
			
			public:
				BindingTransform(TemplateMapping *tm1, TemplateMapping *tm2, ReactionClass *r);
				virtual ~BindingTransform() {};
				
				virtual void transform(Mapping **m);
	};
	
	
	class UnbindingTransform : public Transformation {
				
			public:
				UnbindingTransform(TemplateMapping *tm, ReactionClass *r);
				virtual ~UnbindingTransform() {};
				virtual void transform(Mapping **m);
	};



}




















#endif /*TRANSFORMATION_HH_*/
