#ifndef ENVIRONMENT_HH_
#define ENVIRONMENT_HH_


#include <math.h>

#include "../cell/cell.hh"



namespace NG_swim {
	
	class Cell;
	
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

	

	
	
//	
//	class LinearGradientEnvironment : public Environment
//	{
//		public:
//			
//			LinearGradientEnvironment(double ligandConc);
//			~LinearGradientEnvironment() {};
//			
//			virtual double getLigConc(double xPos, double yPos, double zPos, double time);
//		
//		protected:
//		
//	};
//	
//	class ExponentialGradientEnvironment : public Environment
//	{
//		public:
//			
//			ExponentialGradientEnvironment(double ligandConc);
//			~ExponentialGradientEnvironment() {};
//			
//			virtual double getLigConc(double xPos, double yPos, double zPos, double time);
//		
//		protected:
//	};

	
	
	

	
	
	
	
	class CapillaryEnvironment : public Environment
	{
		public:
			CapillaryEnvironment();
			CapillaryEnvironment(
								double mouthZcoordinate,
								double initialCapillaryConc,
								double diffusionConstant,
								double radiusOfCapillary
								);
			~CapillaryEnvironment() {};
			
			virtual double getLigConc(double xPos, double yPos, double zPos, double time);
			
			virtual void getStartPosition(int cellNumber, double pos[3]);
			
			void checkPos();
			
			
			virtual void tryToMove(
					double p[3],
					double p2[3],
					double u[3],
					double up[3]);
			
			bool updatePosition(
					double p[3],
					double p2[3],
					double u[3],
					double up[3]);
					
			
		protected:
			
			double initialCapillaryConc;
			double diffusionConstant;
			
			
			//Needed parameters
			double Zt;	//Z coordinate of the top of the tube
			double Nt[3];
			double Zb;		//Z coordinate of the base of the tube
			double Nb[3];
			double Zc;   //Z coordinate of the mouth of the capillary
			double Re;	//Radius of the entire environment
			double Rc;	//Radius of the capillary
			
			static const int x=0,y=1,z=2;
	};
	
	
	

	
	
}


#endif /*ENVIRONMENT_HH_*/
