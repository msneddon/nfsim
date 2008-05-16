#ifndef TRANSFORMATIONSET_HH_
#define TRANSFORMATIONSET_HH_


#include "../NFreactions.hh"


namespace NFcore
{

	class Transformation;
	class mappingSets;
	
	class TemplateMolecule;
	class Molecule;


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


}












#endif /*TRANSFORMATIONSET_HH_*/
