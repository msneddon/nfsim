
#ifndef TEMPLATEMOLECULE_HH_
#define TEMPLATEMOLECULE_HH_


#include "NFcore.hh"

namespace NFcore
{

	class MoleculeType;
	class MapGenerator;
	class Molecule;
	class MappingSet;
	class ReactantContainer;

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
		TemplateMolecule(MoleculeType * moleculeType);
		~TemplateMolecule();


		/* get functions */
		MoleculeType *getMoleculeType() const {return moleculeType;};
		string getMoleculeTypeName() const;

		int getN_symComps() const {return n_symComps; };
		int getN_symCompBonds() const {
			int symCompBondCounter=0;
			for(int i=0; i<n_symComps; i++) {
				if(symBondPartner[i]!=0) symCompBondCounter++;
			}
			return symCompBondCounter;
		}
		int getN_mapGenerators() const { return n_mapGenerators; }
		int getN_connectedTo() const { return n_connectedTo; };

		/* functions that allow you to set constraints */
		void addEmptyComponent(string cName);
		void addBoundComponent(string cName);
		void addComponentConstraint(string cName, string stateName);
		void addComponentConstraint(string cName, int stateValue);
		void addComponentExclusion(string cName, string stateName);
		void addComponentExclusion(string cName, int stateValue);
		void addBond(string thisBsiteName,TemplateMolecule *t2, string bSiteName2);

		/* Methods for adding a disjoint component to a template pattern
		 *  e.g.  X.Y
		 */
		void addConnectedTo(TemplateMolecule *t2, int otherConToIndex);
		void addConnectedTo(TemplateMolecule *t2, int otherConToIndex,bool otherHasRxnCenter);
		void clearConnectedTo();


		/* functions that allow you to set constraints for symmetric sites */
		const static int EMPTY=0;
		const static int OCCUPIED=1;
		const static int NO_CONSTRAINT=-1;
		void addSymCompConstraint(string cName, string uniqueId,
				int bondState,int stateConstraint);
		void addSymBond(string thisBsiteName, string thisCompId,
				TemplateMolecule *t2, string bSiteName2);

		/* static function for binding two templates together */
		static void bind(TemplateMolecule *t1, string bSiteName1, string compId1,
				TemplateMolecule *t2, string bSiteName2, string compId2);

		/* functions that provide mapping capabilities */
		void addMapGenerator(MapGenerator *mg);

		/* functions that are needed to perform TemplateMolecule operations */
		bool contains(TemplateMolecule *tempMol);

		const static bool FIND_ALL = false;
		const static bool SKIP_CONNECTED_TO = true;
		static void traverse(TemplateMolecule *tempMol, vector <TemplateMolecule *> &tmList, bool skipConnectedTo);

		/* searches the list of template molecules and identifies the number of disjoint
		   sets, and also returns the mapping onto those sets*/
		static int getNumDisjointSets(vector < TemplateMolecule * > &tMolecules,
				vector <vector <TemplateMolecule *> > &sets,
				vector <int> &uniqueSetId);

		/* functions that are needed to match to a molecule instance */
		bool compare(Molecule *m);
		bool compare(Molecule *m, ReactantContainer *rc, MappingSet *ms,bool holdMolClearToEnd=false);
		void clear();
		void clearTemplateOnly();
		bool tryToMap(Molecule *toMap, string toMapComponent,
				Molecule *mappedFrom, string mappedFromComponent);
		bool isSymMapValid();


		//////////////////////////////////////////////////////////////////////////////////////////////
		//returns false if they are not symmetric, or true if they are
		static bool checkSymmetry(TemplateMolecule *tm1, TemplateMolecule *tm2, string bSite1, string bSite2);
		static bool checkSymmetryAroundBond(TemplateMolecule *tm1, TemplateMolecule *tm2, string bSite1, string bSite2);



        /* functions that handle output for debugging and error messages */
		void printErrorAndExit(string message);
		void printDetails();
		void printDetails(ostream &o);


		string getPatternString();
		void printPattern();
		void printPattern(ostream &o);



	protected:

		static int TotalTemplateMoleculeCount;

		MoleculeType *moleculeType;
		int uniqueTemplateID;

		// Handling of transformations
		int n_mapGenerators;
		MapGenerator **mapGenerators;


		///////////////////////////////
		////  There are two classes of things we have to match that must
		////  be handled separately...
		////  1) unique components
		////  2) symmetric components


		// Which of the unique components must be empty (no bonds)
		int n_emptyComps;
		int *emptyComps;

		// Which of the unique components must be occupied (bonded to something, something
		// that is not specified)
		int n_occupiedComps;
		int *occupiedComps;

		// State value constraints
		int n_compStateConstraint;
		int *compStateConstraint_Comp; //index of the constrained component
		int *compStateConstraint_Constraint; //the constrained value

		// State value exclusions (state != exclusion)
		int n_compStateExclusion;
		int *compStateExclusion_Comp;
		int *compStateExclusion_Exclusion;

		// The set of connections that a particular site is connected to
		int n_bonds;
		int *bondComp;
		string *bondCompName;
		TemplateMolecule **bondPartner;
		string *bondPartnerCompName; //used if nonsymmetric bond is connected to partner symmetric site
		int *bondPartnerCompIndex; //used if nonsymmetric bond is connected to partner nonsymmetric site else =-1
		bool *hasVisitedBond;


		//This stores disjoint sets, in other words, this Template is
		//connected to some other Template via .. the dot operator "."
		int n_connectedTo;
		TemplateMolecule ** connectedTo;
		bool *hasTraversedDownConnectedTo;
		int *otherTemplateConnectedToIndex;
		bool *connectedToHasRxnCenter;


		//////////  Handling symmetric components
		int n_symComps;
		string *symCompName;
		string *symCompUniqueId; //Used to match up a particular component when creating bonds
		int *symCompStateConstraint;
		int *symCompBoundState;  //either Empty (0), Occupied (1), or No constraint(2)
		TemplateMolecule **symBondPartner; //the bound template, if this component is bound
		string *symBondPartnerCompName;
		int *symBondPartnerCompIndex;
		vector < vector <int> > canBeMappedTo; //might want to change this to a 2d array for memory/speed?
		bool *hasTraversedDownSym;

		//Used when matching to a given molecule
		int n_totalComps;
		bool *isSymCompMapped;
		bool *compIsAlwaysMapped;

		Molecule *matchMolecule;
		bool hasVisitedThis;


		//For depth first traversals on a template molecule
		static queue <TemplateMolecule *> q;
		static queue <int> d;
		static vector <TemplateMolecule *>::iterator tmVecIter;
		static list <TemplateMolecule *>::iterator tmIter;

	};

}


#endif
