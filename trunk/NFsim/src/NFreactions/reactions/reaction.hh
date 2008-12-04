#ifndef BASICREACTIONS_HH_
#define BASICREACTIONS_HH_


#include "../NFreactions.hh"




using namespace std;

namespace NFcore
{

	class BasicRxnClass : public ReactionClass {
		public:
			BasicRxnClass(string name, double baseRate, TransformationSet *transformationSet);
			virtual ~BasicRxnClass();

			virtual void init();
			virtual void prepareForSimulation();
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos);
			virtual double update_a() {
					a = 1;

					for(unsigned int i=0; i<n_reactants; i++)
						a*=reactantLists[i]->size();

					a*=baseRate;
					return a;
			}


			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual unsigned int getReactantCount(unsigned int reactantIndex) const;

			virtual void printFullDetails() const;

		protected:
			virtual void pickMappingSets(double randNumber) const;

			ReactantList **reactantLists;

			ReactantList *rl;
			MappingSet *ms;
	};



	inline
	bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
		{
			//First a bit of error checking, that you should skip unless we are debugging...
		//	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
		//	{
		//		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		//		exit(1);
		//	}


			//Get the specified reactantList
			rl = reactantLists[reactantPos];

			//Check if the molecule is in this list
			int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
			//cout<<"got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;


			if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
			{
				if(!reactantTemplates[reactantPos]->compare(m)) {
				//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
					rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
					m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
				}

			} else {
				//Try to map it!
				ms = rl->pushNextAvailableMappingSet();
				if(!reactantTemplates[reactantPos]->compare(m,ms)) {
					rl->popLastMappingSet();
					//we just pushed, then popped, so we a has not changed...
				} else {
					m->setRxnListMappingId(rxnIndex,ms->getId());
				}
			}

			return true;
		}



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
			MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet );
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
					TransformationSet *transformationSet,
					vector <LocalFunction *> &lfList,
					vector <string> &lfArgumentPointerNameList);
			virtual ~DORRxnClass();

			virtual void init();
			virtual void prepareForSimulation() {};
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
			virtual void remove(Molecule *m, unsigned int reactantPos) {cout<<"calling remove in DOR?"<<endl; exit(0);};
			virtual double update_a();


			virtual int getDORreactantPosition() const { return DORreactantIndex; };



			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual unsigned int getReactantCount(unsigned int reactantIndex) const;

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

			//Parameters to keep track of local functions
			int DORreactantIndex;
			vector <LocalFunction *> lfList;
			vector <int> indexIntoMappingSet;
			vector <double> localFunctionValue;

//			LocalFunction
	};



}







#endif /*BASICREACTIONS_HH_*/
