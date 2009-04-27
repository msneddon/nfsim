#ifndef NG_CELL_HH_
#define NG_CELL_HH_

#include <math.h>
#include <fstream>

#include "../../../NFcore/NFcore.hh"
#include "../../../NFutil/NFutil.hh"
#include "util.hh"

#include "environment.hh"



using namespace std;
using namespace NFcore;


class Environment;

/** The basic unit of the NG chemotaxis simulation.
 * This is the core class needed to run simulations of swimming cells in a
 * 3D environment.  It is basically a wrapper for a NG network free simulator
 * object (@see ng.hh) that has additional functions that allow it to swim
 * according to the correct E.Coli statistics.  It allows interaction with
 * the given environment, but does not yet have the capability to interact
 * with other cells yet.
 */
class AgentCell
{
	public:

		AgentCell(
				System *s,
				Environment *e,
				double speed,
				double rotDiffusionConstant,
				double cheYpThreshold,
				double startTime,
				string outputDirectoryPath
				);

		~AgentCell();

		double stepTo(double endTime, double dt);
		void equilibriate(double duration, double dt);

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

		static const int APART = 0;
		static const int BUNDLED = 1;

		void swimToNewPosition(double elapsedTime);


		void outputCellValues();

	protected:


		void changeDirRandom();
		void changeDirDistribution();
		void changeDirRotDiffusion(double elapsedTime);
		void moveToNewPosition(double elapsedTime);

		double speed;
		double rotDiffusionConstant;




		double * cheYhistory;
		int cheYhisPos;

		int cheYpThreshold;// = 1795;
		int motorState;// = this->CCW;
		int flagellaState;// = this->APART;
		int lastFlagellaState;
		double apartDuration;// = 0;
		bool droppingTumble;// = false;

		double boxcarTimeWidth;// = 0.3;
		int cheYhistorySize;// = (int)round(boxcarTimeWidth/dt);

		double meanCheYp;// = 0;
		double cheYhistorySum;// = 0;




		//Identify the cell with a number
		int cellId;

		//Keep track of current 3d position, direction, and local orientation
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


		//Distribution of new swim direction angles, matched to experimental observations
		// With shape = 4, scale = 18.32045567939674, location = -4.606176605609249
		// These variables are now initialized in cell.cpp
		static const double DefaultGammaAlpha;// = 4.;
		static const double DefaultGammaBeta;// = 1./18.32045567939674;
		static const double DefaultGammaOffset;// = -4.606176605609249;

		static const double DefaultSpeed;// = 20; //uM per second
		static const double DefaultRotDifConst;// = 0.06205; // rad^2/s



		//Keep a copy of the chemotaxis system and environment
		System * system;
		Environment *env;

		//Also keep track of the output stream for the cell properties
		ofstream motorFileStream;
		ofstream outputFileStream;
		void outputCellHeader();

		double runTimeSpent;


		string fileName;
		unsigned int binaryFileOutputCounter;

		static const int CW = 1;
		static const int CCW = 0;



	private:
		static unsigned int cellIdCounter;
		ChemotaxisUtil::GammaSampler *gs;

};


#endif /*NG_CELL_HH_*/

