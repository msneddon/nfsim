#ifndef NFRATECHANGEREACTIONS_HH_
#define NFRATECHANGEREACTIONS_HH_

#include "../../NFcore/NFcore.hh"



class DORrxn;
class NFReactantTree;

using namespace NFcore;

////////////////////////////////////////////////////////////////////
///////////////  Generalized DOR reaction


class DORrxn : public ReactionClass {
	
	public:
	
		/* DOR rxn constructor and deconstructor */
		DORrxn( char * name, 
				int n_reactants, 
				TemplateMolecule ** reactantTemplates,
				int DORreactantIndex,
				char * DORgroupName,
				int DORgroupValueIndex,
				double baseRate);
		virtual ~DORrxn();
								
		virtual void init();
		virtual double update_a();
		
		virtual Molecule ** pickReactants(double randA_Value);
		virtual void transformReactants(Molecule ** reactants, int nReactants);
		
		virtual int addToReactantList(Molecule *m, int reactantIndex);
		virtual void removeFromReactantList(Molecule *m, int reactantIndex, int rxnListIndex);
		
		virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
		
		virtual void prepareForSimulation();
		
		virtual void printDetails() const;
		
	protected:
		int stateIndex;
		int newStateValue;
		
		int DORreactantIndex;
		char * DORgroupName;
		int DORgroupValueIndex;
		
		NFReactantTree * dorReactants;
		
		double * rateFactorSum;
		double ** rateFactors;
};


/* Observable Dependent Rxn 
 * 
 *  - A reaction whose rate depends on the value of an observable
 * 
 * 
 * */

class ObsDrxn : public ReactionClass {
	
	public:
	
		/* DOR rxn constructor and deconstructor */
		ObsDrxn( char * name, 
				int n_reactants, 
				TemplateMolecule ** reactantTemplates,
				int n_observables,
				Observable ** obs);
		virtual ~ObsDrxn();
								
		virtual double update_a();
		
		virtual void transformReactants(Molecule ** reactants, int nReactants);
		virtual void printDetails() const;
		
	protected:
		int n_observables;
		Observable ** obs;
};



class NFReactantTree {
	
	
	public:
		NFReactantTree(unsigned long int maxElementCount);
		~NFReactantTree();
		
		//Once you insert into this tree, you will get the position you were
		//inserted into.  This will never change until you are removed
		int insert(Molecule * m, double rateFactor);
		void remove(Molecule * m, unsigned int rxnListIndex);
		void updateValue(Molecule * m, unsigned int rxnListIndex, double newRateFactor);
		
		Molecule * getReactantFromValue(double value, double baseRate) const;
		
		
		int getNumOfMolecules() const { return numOfMolecules; };
		double getRateFactorSum() const { return leftRateFactorSum[0]; };
		
		int getDepthOfTree() const { return treeDepth; };
		void printDetails() const;
	
	
	protected:
	
		//Basic tree parameters and constants
		unsigned int maxElementCount;
		unsigned int treeDepth;
		unsigned int treeFirstMoleculeIndex;
		unsigned int numOfNodes;
		
		//The tree stored as 3 arrays indexed as:
		// Children are at 2i and 2i+1
		// Parent is at i/2 (using integer division)
		// root is at index 1
		// index 0 is empty and always zero
		double * leftRateFactorSum;
		unsigned int * leftElementCount;
		unsigned int * rightElementCount;
		
		//The actual list of reactants
		Molecule ** reactants;
	
		//The number of molecules currently in the tree
		unsigned int numOfMolecules;
		
		//The index of the first molecule in the tree
		unsigned int firstMoleculeTreeIndex;
	
};








#endif /*NFRATECHANGEREACTIONS_HH_*/
