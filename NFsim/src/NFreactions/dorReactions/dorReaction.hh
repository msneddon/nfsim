//#ifndef DORREACTIONS_HH_
//#define DORREACTIONS_HH_
//
//#include "../NFreactions.hh"
//
//
//using namespace std;
//
//namespace NFcore
//{
//
//	class DORrxnClass : public ReactionClass {
//		public:
//			DORrxnClass(
//					string name, 
//					double baseRate, 
//					TransformationSet *transformationSet, 
//					unsigned int DORreactantIndex, 
//					string DORgroupName, 
//					int DORgroupValueIndex);
//			virtual ~DORrxnClass();
//			
//			
//			//The main virtual functions that must be implemented in all implementing classes
//			virtual void init(); //called when the reaction is added to the system
//			virtual void prepareForSimulation(); //called once everything is added to the system
//			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
//			virtual void remove(Molecule *m, unsigned int reactantPos);
//			
//			virtual double update_a();
//			
//			
//			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
//			
//			
//			
//			virtual unsigned int getReactantCount(unsigned int reactantIndex) const;
//			virtual void printFullDetails() const;
//			
//		protected:
//			virtual void pickMappingSets(double random_A_number) const;
//			
//			int stateIndex;
//			int newStateValue;
//			
//			int DORreactantIndex;
//			string DORgroupName;
//			int DORgroupValueIndex;
//			
//			ReactantTree * reactantTree;  
//			ReactantList ** reactantLists;
//			
//			double * rateFactorSum;
//			double ** rateFactors;
//	};
//	
//};
//
//
//
//
//
//
//
//
//
//
//
//
//
//#endif /*DORREACTIONS_HH_*/
