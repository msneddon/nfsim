#ifndef NFREACTIONS_HH_
#define NFREACTIONS_HH_

#include <vector>
#include <string>

#include "../NFcore/NFcore.hh"
#include "../NFtest/simple_system/simple_system.hh"





using namespace std;



namespace NFcore
{
	//Forward Declarations
	class TemplateMolecule;
	class Molecule;




	void test();
	
	
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
			
			
			
			bool addBindingTransform() { return false; };
			bool addUnbindingTransform() { return false; };
			bool addDeleteMolecule() { return false; };
			bool addAddMolecule() { return false; };
			
			
			bool transform(MappingSet **mappingSets);
			MappingSet *generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId);
		
		protected:
			
			
			int find(TemplateMolecule *t);
			
			unsigned int n_reactants;
			TemplateMolecule ** reactants;
			
			
			vector <Transformation *> *transformations;
	};
	
	
	
	class Transformation
	{
		public:
			
			Transformation();
			~Transformation();
			unsigned int getType() const { return type; };
			int getNewStateValue() const { return newStateValue; };
			unsigned int getBiPartnerReactantIndex() const { return reactant; };
			unsigned int getBiPartnerMappingIndex() const { return mappingIndex; };
			
			static Transformation * genStateChangeTransform(int newStateValue);
			static Transformation * genBindingTransform1(unsigned int reactant, unsigned int mappingIndex);
			static Transformation * genBindingTransform2();
			static Transformation * genUnbindingTransform();
			static Transformation * genAddMoleculeTransform();
			static Transformation * genRemoveMoleculeTransform();
			
			//Keep track of types
			static const unsigned int STATE_CHANGE = 0;
			static const unsigned int BINDING = 1;
			static const unsigned int UNBINDING = 2;
			static const unsigned int REMOVE = 3;
			static const unsigned int ADD = 4;
			static const unsigned int SKIP = 5; //used to mark a transformation placeholder that we don't actually want to execute, as
						          
			
		protected:
			unsigned int type;
			
			//For a state change transformation
			int newStateValue;
			
			//For a binding reaction
			unsigned int reactant;
			unsigned int mappingIndex;
			
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
			MappingSet(vector <Transformation *> &transformations);
			~MappingSet();
			
			bool set(unsigned int mappingIndex, Molecule *m);
			Mapping *get(unsigned int mappingIndex);
			bool clear();
			
			
		protected:
			unsigned int id;
			
			bool isSet;
			unsigned int n_mappings;
			Mapping * mappings;
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
			
			
			const static unsigned int STATE =  0;
			const static unsigned int BSITE =  1;
			const static unsigned int REMOVAL =  2;
			const static unsigned int ADDITION =  3;
			

		protected:
			unsigned int type;
			unsigned int index;
			Molecule * m;
	};
	
	
	
	//!  Maintains a list of Mapping & Molecule objects needed by ReactionClass
		/*!
		  This is essentially a specialized vector implementation that allows a ReactionClass 
		  to easily keep track of all the Molecule objects that can be involved in the reaction.  It
		  also has to maintain Mapping objects into those molecules (so that transformations can
		  easily be made through the Transformation class).  This class automatically expands 
		  its capacity when extra mappings are added.  The advantage over the traditional vector 
		  class is that this class allows for near constant time removal and insertion of elements, 
		  while vectors require linear time removal.  We can gain this speedup because the indexing 
		  of the reactant list is unimportant.
		    @author Michael Sneddon
		 */
		class ReactantList
		{
			
			
			public:
				//!  Default Constructor
				/*!
					  Creates a new empty ReactantList with the given initial capacity.  This capacity
					  should roughly be set to the number of mappings you expect this list to have.
				 */
				ReactantList(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity);
				~ReactantList();
			
			
				unsigned int size() const;
				
				
				MappingSet * pushNextAvailableMappingSet();
				void popLastMappingSet();
				void removeMappingSet(unsigned int mappingSetId);
			
				
				
				
				
				void printDetails();
			
			protected:
				
				
				unsigned int n_mappingSets;
				unsigned int capacity;
				TransformationSet *ts;
				unsigned int reactantIndex;
				
				unsigned int * msPositionMap;
				MappingSet **mappingSets;
		};
		
	
	
}

#endif /*NFREACTIONS_HH_*/
