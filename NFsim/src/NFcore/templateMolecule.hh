
#ifndef TEMPLATEMOLECULE_HH_
#define TEMPLATEMOLECULE_HH_


#include "NFcore.hh"

namespace NFcore
{

	class MoleculeType;
	class MapGenerator;
	class Molecule;
	class MappingSet;


	//!  Used for matching Molecule objects to the given pattern
	/*!
	    TemplateMolecules are regular expression like objects needed to identify specific
	    configurations of connected Molecules.  Individual TemplateMolecules are derived
	    from a particular MoleculeType and inherit their set of components from their parent
	    MoleculeType.  TemplateMolecules can be connected to other TemplateMolecules through
	    component bonds forming the regular expression pattern.  An individual Molecule matches
	    an individual TemplateMolecule if they are of the same MoleculeType and their
	    component bonds and state values match.  Any component state not explicitly specified
	    in a TemplateMolecule is treated as a wild-card and will always match the corresponding
	    component of a Molecule.  For more complex patterns, an entire connected set of Molecules
	    is matched to a connected set of TemplateMolecules through a recursive algorithm that
	    checks for graph isomorphism between the two sets.  The worst case performance of the
	    recursive matching algorithm is proportional to the number of connected TemplateMolecules
	    in the pattern.  However, the average performance is much better because a match is rejected
	    as soon as a single difference in component states or molecule connectivity is found.
        @author Michael Sneddon
	 */
	class TemplateMolecule {
	public:
		TemplateMolecule(MoleculeType * parentMoleculeType);
		~TemplateMolecule();

		/* get functions */
		MoleculeType * getMoleculeType() const { return parentMoleculeType; };
		string getMoleculeTypeName() const;

		unsigned int getNumCompareStates() const { return stateIndex.size(); };
		unsigned int getStateIndex(int state) const { return stateIndex.at(state); };
		unsigned int getStateValue(int state) const { return stateValue.at(state); };

		unsigned int getNumBindingSites() const { return bonds.size(); };
		bool isBindingSiteOpen(int bIndex) const;// { return bonds.at(bIndex)->isOpen(); };
		bool isBindingSiteBonded(int bIndex) const; //{ return bonds.at(bIndex)->isBonded(); };
		TemplateMolecule * getBondedTemplateMolecule(int bIndex) const;
		unsigned int getTemplateBsiteIndexFromMoleculeBsiteIndex(int molBsiteIndex);

		int getBindingSiteIndex(int bIndex) const;
		int getBindingSiteIndexOfBondedTemplate(int bIndex) const ;


		//int compareAll(TemplateMolecule *tm1, TemplateMolecule *tm2);

		/* set functions */
		void setHasVisited(int bSite);


		// add constraints....
		int addEmptyBindingSite(string bSiteName);

		void addOccupiedBindingSite(string bSiteName);
		static void bind(TemplateMolecule *t1, int bSiteIndex1, TemplateMolecule *t2, int bSiteIndex2);
		static void bind(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2);

		void addStateValue(int cIndex, int stateValue);
		void addStateValue(string cName, int stateValue);
		void addStateValue(string cName, string stateValue);
		void addNotStateValue(string stateName, int notStateValue);



		void clear() {
			this->matchMolecule = 0;
			for(unsigned int i=0; i<hasVisitedBond.size(); i++) hasVisitedBond.at(i) = false;
			hasVisited=false; };


		/* the primary function and purpose of a template molecule
			   is to compare itself to an instance of a molecule */
		bool compare(Molecule * m);
		static bool compareBreadthFirst(TemplateMolecule *tm, Molecule *m);

		//general function used to get a list of all the template molecules connected to
		//this template molecule through bonds.
		static void traverse(TemplateMolecule *tempMol, vector <TemplateMolecule *> &tmList);

		/*
		 * Used to check if a particular state value matches - used when parsing an xml file and
		 * generating reaction transformations.
		 */
		bool isStateValue(string stateName, int stateValue);
		bool isBonded(string bSiteName);


	//	void addTemplateMapping(TemplateMapping *tm);
	//	bool compare(Molecule * m, MappingSet *mappingSet);



		void printDetails() const;


		static const int NO_CONSTRAINT = -1;
		static const int IS_EMPTY = 0;
		static const int IS_OCCUPIED=1;
		void addSymComponentConstraint(string name, int bondState, int stateConstraint);

		///////////////////////////////////////////////////////////////////
		void addMapGenerator(MapGenerator *mg);
		bool compare(Molecule *m, MappingSet *ms);
		bool contains(TemplateMolecule *tempMol);


		Molecule * matchMolecule;
		vector <bool> hasVisitedBond; //Change this to array for slight performance gain...


		/////////////////////////////////////////////////////////
		bool hasVisited;

	protected:


		int addEmptyBindingSite(int bSiteIndex);

		MoleculeType *parentMoleculeType;  // ptr to indicate type of molecule
		vector <int> stateIndex; // saves index into state array of MoleculeType
		vector <int> stateValue; // saves the value we need to have in the state array

		vector <int> notStateIndex;
		vector <int> notStateValue;

		vector <int> bSiteIndex; // saves index into bSite array in MoleculeType
		vector <TemplateMolecule *> bonds; // tells us a binding site
		vector <int> bSiteIndexOfBond;

		vector <int> sitesThatMustBeOccupied; //



		vector <MapGenerator *> mapGenerators;
		vector <MapGenerator *>::iterator mgIter;
		vector <int>::iterator intVecIter;


		static queue <TemplateMolecule *> tmq;
		static queue <Molecule *> mq;
		static list <TemplateMolecule *> tml;
		static queue <int> d;
		static list <TemplateMolecule *>::iterator tmIter;




		bool checkBasicSymConstraints(Molecule *m);



		bool hasSymmetricConstraint;
		int nComponents;
		bool * isComponentMapped;
		bool * componentIsAlwaysMapped;

		vector <string> symComponentNames;
		vector <int> symComponentBondState;
		vector <int> symComponentConstraintValue;

		vector <TemplateMolecule *> symBondedTo;

		// @TODO add not constraint value on sym states


		//vector <string> symEmptyBindingSite;
		//vector <string> symOccupiedBindingSite;

		//vector <string> symStateConstraint;
		//vector <int> symStateConstraintValue;

		//vector <string> symNotStateConstraint;
		//vector <int> symNotStateConstraintValue;

	};

}
#endif
