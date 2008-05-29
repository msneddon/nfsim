#ifndef CELL_HH_
#define CELL_HH_


#include <math.h>
#include <fstream>

#include "../../../../NFcore/NFcore.hh"
#include "../../../../NFutil/NFutil.hh"
#include "../../ng_param/ng_param.hh"
#include "../../ng_creator/ng_creator.hh"
//#include "../../ng/ng.hh"
#include "../../util/chemotaxisUtil.hh"


#include "../environment/environment.hh"





namespace NG_swim {

class Environment;



/** The basic unit of the NG chemotaxis simulation.
 * This is the core class needed to run simulations of swimming cells in a
 * 3D environment.  It is basically a wrapper for a NG network free simulator
 * object (@see ng.hh) that has additional functions that allow it to swim
 * according to the correct E.Coli statistics.  It allows interaction with
 * the given environment, but does not yet have the capability to interact
 * with other cells yet.
 */
class Cell
{
	public:
		
		Cell(	
				Environment *e,
				double speed,
				double rotDiffusionConstant,
				double startTime,
				const char * outputDirectoryPath,
				NGparam &p
				);
		
		~Cell();
		
		double stepTo(double endTime, double dt);
		void equilibriate(double duration);
		
		void turnOnDirOutput() { shouldOutputDirection=true; };
		void turnOffDirOutput() { shouldOutputDirection=false; };
		
		//Functions to get current position and direction
		double getXposition() const { return pos[X]; };
		double getYposition() const { return pos[Y]; };
		double getZposition() const { return pos[Z]; };
		double getXdirection() const { return dir[X]; };
		double getYdirection() const { return dir[Y]; };
		double getZdirection() const { return dir[Z]; };
		
		//Constants that can be used
		static const int X = 0;
		static const int Y = 1;
		static const int Z = 2;
		
		static const int TUMBLE = 0;
		static const int SWIM = 1;
		void swimToNewPosition(double elapsedTime);
	
	protected:
		
		
		void changeDirRandom();
		void changeDirDistribution();
		void changeDirRotDiffusion(double elapsedTime);
		void moveToNewPosition(double elapsedTime);
		
		double speed;
		double rotDiffusionConstant;
		
		
		//Identify the cell with a number
		int cellId;
	
		//Keep track of curent 3d position, direction, and local orientation
		// (could use a quanternion and save less numbers, but I like my intuitive vectors...)
		//  also, I ran a test using double values on a 64bit linux machine to test the error
		//  that develops between the orthogonality between the up and dir vectors.  Running
		//  10 million random rotations on dir and up, I get only a ~.005 error in the degree (meaning
		//  the degree between the two vectors was ~90.005 instead of 90 exactly).  I believe this to
		//  be a reasonable error, as if we apply a new rotation every .01 seconds in the simulation,
		//  we only get an error that is greater than about 0.005 degrees after 27.7 simulation hours,
		//  which is ok for me...  (note, in the same test, the magnitude of the vectors changes by
		//  less than 0.00001 (they were 0.999991 instead of exactly 1)
		double pos[3];
		double dir[3];
		double up[3];
		double rotMat[3][3];
		
		double lastPos[3];
		double lastDir[3];
		double lastUp[3];
		
		
		//Cell properties
		
		double currentTime;
		int currentMovement;
		int lastMovement;
		double currentLigConc;
		
		//Cell Movement Stats
		bool shouldOutputDirection;
		
		
		
		// With shape = 4 scale = 18.32045567939674 location = -4.606176605609249 it
		static const double DefaultGammaAlpha = 4.;
		static const double DefaultGammaBeta = 1./18.32045567939674;
		static const double DefaultGammaOffset = -4.606176605609249;
		
		static const double DefaultSpeed = 20; //uM per second
		static const double DefaultRotDifConst = 0.06205; // rad^2/s
		
		//protected update functions

		
		
		
		//Keep a copy of the chemotaxis system and environment
		System * chemotaxisSystem;
		Environment *env;
		
		//Also keep track of the output stream for the cell properties
		ofstream outputFileStream;
		void outputCellHeader();
		void outputCellValues();
		
		double runTimeSpent;
		
	private:
		static unsigned int cellIdCounter;
		ChemotaxisUtil::GammaSampler *gs;
};

}


#endif /*CELL_HH_*/
