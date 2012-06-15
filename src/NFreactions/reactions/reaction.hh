#ifndef BASICREACTIONS_HH_
#define BASICREACTIONS_HH_


#include "../NFreactions.hh"



using namespace std;

namespace NFcore
{

	class BasicRxnClass : public ReactionClass {
		public:
			BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s);
			virtual ~BasicRxnClass();

			virtual void init();
			virtual void prepareForSimulation();
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printFullDetails() const;

		protected:
			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;

			ReactantList *rl;
			MappingSet *ms;
	};


	class FunctionalRxnClass : public BasicRxnClass {

		public:
			FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s);
			FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s);

			virtual ~FunctionalRxnClass();

			virtual double update_a();
			virtual void printDetails() const;

		protected:
			GlobalFunction *gf;
			CompositeFunction *cf;
	};

	class MMRxnClass : public BasicRxnClass {

		public:
			MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet, System *s);
			virtual ~MMRxnClass();

			virtual double update_a();
			virtual void printDetails() const;

		protected:
			double Km;
			double kcat;
			double sFree;
	};


	class DORRxnClass : public ReactionClass {
		public:
			DORRxnClass(
					string name,
					double baseRate,
					string baseRateName,
					TransformationSet *transformationSet,
					CompositeFunction *function,
					vector <string> &lfArgumentPointerNameList,
					System *s);
			virtual ~DORRxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a();

			virtual int getDORreactantPosition() const { return DORreactantIndex; };

			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual int getReactantCount(unsigned int reactantIndex) const;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const;

			virtual void printDetails() const;
			virtual void printFullDetails() const {};

			void directAddForDebugging(Molecule *m);
			void printTreeForDebugging();

			static void test1(System *s);

		protected:

			virtual double evaluateLocalFunctions(MappingSet *ms);

			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;
			ReactantTree *reactantTree;

			MappingSet *ms;

			CompositeFunction *cf;

			//Parameters to keep track of local functions
			int DORreactantIndex;

			int n_argMolecules;
			int * argIndexIntoMappingSet;
			Molecule ** argMappedMolecule;
			int * argScope;

			//vector <int> argIndexIntoMappingSet;

			//vector <LocalFunction *> lfList;
			//vector <int> indexIntoMappingSet;
			//vector <double> localFunctionValue;

	};

}







#endif /*BASICREACTIONS_HH_*/
