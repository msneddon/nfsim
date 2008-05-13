#ifndef MAPPING_HH_
#define MAPPING_HH_


#include "../../NFcore/NFcore.hh"
#include "../transformation/transformation.hh"

#include <vector>
#include <map>
using namespace std;


namespace NFcore {

	//Forward Declarations
	class TemplateMapping;
	class MappingSet;
	class Mapping;
	class Transformation;
	
	class System;
	class MoleculeType;
	class Molecule;
	class TemplateMolecule; 
	class ReactionClass;
	class Observable;			
	

	class TemplateMapping {
		
		public:
			
			TemplateMapping(TemplateMolecule *tm, unsigned int type, const char *name);
			TemplateMapping(TemplateMolecule *tm, unsigned int type, unsigned int index);
			~TemplateMapping() {}
			
			
			Mapping *createNewMapping(Molecule *m);
			void setTransformation(Transformation *t);
			MoleculeType * getMoleculeType() const { return moleculeType; };
			unsigned int getMappingType() const { return type; };
			unsigned int getMappingIndex() const { return index; };
			
			void printDetails();
			
		protected:
			
			//The type of molecule we are mapping to
			MoleculeType *moleculeType;
			
			//The transformation this mapping is associated with
			Transformation *transformation;
			
			// This gives whether we are mapping to a state or to a binding site
			unsigned int type;
			
			// Stores the index of the state / binding site
			unsigned int index;
			
		
	};

	
	class MappingSet {
			
			public:
				MappingSet();
				~MappingSet();
				
				void add(Mapping * m);
				void clear();
				void start();
				Mapping * getNextMapping();
				int getNumOfMappings();
				Mapping * getMapping(unsigned int index);
				
				
				static void transform(vector <MappingSet *> &mappingSets);
				static void getPossibleProducts(vector <MappingSet *> &mappingSets, list <Molecule *> &products, unsigned int traversalLimit);
				
				void printDetails();
				
			protected:
				vector <Mapping *> mappings;
	};
	
	

	class Mapping {
		
		public:
			Mapping(Molecule *m, Transformation *t, unsigned int type, unsigned int index);
			~Mapping() { molecule=NULL; transformation=NULL; };
			
			Molecule * getMolecule() const {return molecule;};
			Transformation * getTransformation() const {return transformation;};
			unsigned int getType() const {return type;};
			unsigned int getIndex() const {return index;};
			void printDetails();
			
			static const unsigned int BOND = 0;
			static const unsigned int STATE = 1;
			
		protected:
			Molecule *molecule;
			Transformation *transformation;
			unsigned int type;
			unsigned int index;
	};





}













#endif /*MAPPING_HH_*/
