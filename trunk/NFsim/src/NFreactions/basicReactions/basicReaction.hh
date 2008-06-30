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
						a*=reactantLists.at(i)->size();
					
					a*=baseRate;
					return a;
			}
						
						
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
			virtual unsigned int getReactantCount(unsigned int reactantIndex) const;
						
			virtual void printFullDetails() const;
		
		protected:
			virtual void pickMappingSets(double randNumber) const;
			vector <ReactantList *> reactantLists;
	};
}







#endif /*BASICREACTIONS_HH_*/
