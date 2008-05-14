#ifndef NFREACTIONS_HH_
#define NFREACTIONS_HH_

#include <vector>
#include <string>

#include "../NFcore/NFcore.hh"


#include "ReactantList.hh"




using namespace std;



namespace NFcore
{
	//Forward Declarations
	class TemplateMolecule;
	class Molecule;
	class ReactionClass;



	void test();
	void test_simple();
	
	
	class Transformation;
	class TransformationSet;
	class MapGenerator;
	class MappingSet;
	class Mapping;
	class ReactantList;


	
	class TransformationSet
	{
		public:
			TransformationSet(vector <TemplateMolecule *> reactantTemplates);
			~TransformationSet();
			
			bool addStateChangeTransform(TemplateMolecule *t, string stateName, int finalStateValue);
			
			
			bool addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2);
			bool addUnbindingTransform(TemplateMolecule *t, string bSiteName);
			
			
			bool addDeleteMolecule() { return false; };
			bool addAddMolecule() { return false; };
			
			
			bool transform(MappingSet **mappingSets);
			MappingSet *generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId);
		
			
			//Used for initializations in ReactionClass
			TemplateMolecule * getTemplateMolecule(unsigned int reactantIndex) const { return reactants[reactantIndex]; };
			unsigned int getNreactants() const {return n_reactants; }; 
			void finalize();
			bool isFinalized() const { return finalized; };
			bool getListOfProducts(MappingSet **mappingSets, list<Molecule *> &products, int traversalLimit);
			
		protected:
			bool finalized;
			
			int find(TemplateMolecule *t);
			
			unsigned int n_reactants;
			TemplateMolecule ** reactants;
			
			
			vector <Transformation *> *transformations;
	};
	
	
	
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
	
	
	
	//!  Keeps a list of mappings
	/*!
	    @author Michael Sneddon
	 */
	class MappingSet
	{
		public:
			MappingSet(unsigned int id, vector <Transformation *> &transformations);
			~MappingSet();
			
			bool set(unsigned int mappingIndex, Molecule *m);
			Mapping *get(unsigned int mappingIndex);
			bool clear();
			
			unsigned int getId() const { return id; };
			
			
		protected:
			unsigned int id;
			
			bool isSet;
			unsigned int n_mappings;
			Mapping ** mappings;
	};

	
	
	
	
	//!  Keeps a pointer to a molecule and remembers a component to act on
	/*!
	    @author Michael Sneddon
	 */
	class Mapping
	{
		public:
			Mapping(unsigned int type, unsigned int index);
			~Mapping();
			
			unsigned int getType() const;
			unsigned int getIndex() const;
			Molecule * getMolecule() const;
			
			void clear();
			bool setMolecule(Molecule *m);

		protected:
			unsigned int type;
			unsigned int index;
			Molecule * m;
	};
	
}

#endif /*NFREACTIONS_HH_*/
