#ifndef NG_ENVIRONMENT_HH_
#define NG_ENVIRONMENT_HH_

#include <math.h>

#include "../cell/cell.hh"


class AgentCell;

	class Environment
	{
		public:

			Environment();
			virtual ~Environment() {};

			virtual double getLigConc(double xPos, double yPos, double zPos, double time);

			virtual void getStartPosition(int cellNumber, double pos[3]);
			//virtual void tryToMove(double lastPosition[3], double position[3], double direction[3], double up[3]);

			virtual void tryToMove(
								double p[3],
								double p2[3],
								double u[3],
								double up[3]);

			virtual bool isInCap(AgentCell *c) { return false; };

		protected:

	};

	class ConstantEnvironment : public Environment
	{
		public:

			ConstantEnvironment(double ligandConc);
			~ConstantEnvironment() {};

			virtual double getLigConc(double xPos, double yPos, double zPos, double time);
			virtual void tryToMove(
					double p[3],
					double p2[3],
					double u[3],
					double up[3]);
			void setLigConc(double ligandConc) {this->ligandConc=ligandConc;};

		protected:
			double ligandConc;
	};




	class LinearEnvironment : public Environment
		{
			public:
	               LinearEnvironment() {
	                    this->slope=pow(10,-6.5);
	                    this->intercept=0;//1e-6;
	               };
				LinearEnvironment(double slope, double intercept) {
	                    this->slope=slope;
	                    this->intercept=intercept;
	               };
				~LinearEnvironment() {};

				virtual double getLigConc(double xPos, double yPos, double zPos, double time) {

	                    double conc = slope*zPos+intercept;
	                    if(conc<0) return 0;
	                    return conc;
	               };

				virtual void getStartPosition(int cellNumber, double pos[3]);

				virtual void tryToMove(
						double p[3],
						double p2[3],
						double u[3],
						double up[3]) {};


			protected:
	               double slope;
	               double intercept;
		};



#endif /*NG_ENVIRONMENT_HH_*/

