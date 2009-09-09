/*
 * reactionSelector.hh
 *
 *  Created on: Jul 23, 2009
 *      Author: msneddon
 */

#ifndef REACTIONSELECTOR_HH_
#define REACTIONSELECTOR_HH_



#include "../NFcore.hh"
//#include <vector>


using namespace std;


namespace NFcore
{
	//Forward Declarations
	class ReactionClass;


	// Abstract Interface Class for the Reaction Selection Algorithm
	class ReactionSelector {

		public:

			//Initializations and basic functionality
			ReactionSelector() {};
			virtual ~ReactionSelector() {};

			virtual double refactorPropensities() = 0;


			virtual double update(ReactionClass *r,double oldA, double newA) = 0;
			virtual double getNextReactionClass(ReactionClass *&rc) = 0;
			virtual double getAtot() = 0;

	};


	class DirectSelector : public ReactionSelector {

		public:
			//Initializations and basic functionality
			DirectSelector(vector <ReactionClass *> &rxns);
			virtual ~DirectSelector();

			virtual double refactorPropensities();


			virtual double update(ReactionClass *r,double oldA, double newA);
			virtual double getNextReactionClass(ReactionClass *&rc);
			virtual double getAtot();


		protected:
			double Atot;
			int n_reactions;
			ReactionClass ** reactionClassList;

	};


	class LogClassSelector : public ReactionSelector {

		public:
			//Initializations and basic functionality
			LogClassSelector(vector <ReactionClass *> &rxns);
			virtual ~LogClassSelector();

			virtual double refactorPropensities();


			virtual double update(ReactionClass *r,double oldA, double newA);
			virtual double getNextReactionClass(ReactionClass *&rc);
			virtual double getAtot();


		protected:

			int calculateClass(double a);

			void place(ReactionClass *r,int logClass,double a);

			void setLogClassToActive(int logClass);
			void setLogClassToInactive(int logClass);


			// The absolute maximum log class
			int maxClassLimit;
			int minClassLimit;
			int totalLogClassCount;

			//the current range of log classes that we have to loop over
			int currentHighClass;
			int currentLowClass;


			// A 2d array of the logClasses
			ReactionClass *** logClassList;

			// A 1d array of the logClasses sizes
			int *logClassSize;
			int *logClassCapacity;


			//An array to keep track of the log classes with positive values
			int *activeLogClasses;
			bool *isLogClassActive;
			int n_activeLogClasses;



			// A 1d array of the propensity sum of the logClass
			double *logClassPropensity;



			int *mapRxnIdToLogClass;
			int *mapRxnIdToLogClassPosition;




			double Atot;
			int n_reactions;
			ReactionClass ** reactionClassList;

	};



}








#endif /* REACTIONSELECTOR_HH_ */
