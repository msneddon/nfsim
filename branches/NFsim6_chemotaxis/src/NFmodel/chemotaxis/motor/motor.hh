#ifndef MOTOR_HH_
#define MOTOR_HH_


#include "../ng/ng.hh"
#include "../chemotaxis.hh"
#include "../ng_param/ng_param.hh"
#include "../../../NFreactions/rate_change/NFrateChangeReactions.hh"


using namespace NFcore;


namespace Motor {


MoleculeType * makeMotor(System *s);


ReactionClass * createMotor_cw2ccw(MoleculeType * Motor, Observable * cheYp, ofstream *o);
ReactionClass * createMotor_ccw2cw(MoleculeType * Motor, Observable * cheYp, ofstream *o);

void ANrxn_addMotorSwitching_Full(System *s, MoleculeType *motor, MoleculeType *cheY, const char * motorStatsOutputFilename);

void createRxn_MotorSwitching(System *s, MoleculeType *motor, MoleculeType *cheY, NGparam &p);
ReactionClass * createMotorRxn_cw2ccw(MoleculeType * Motor, Observable * cheYp, ofstream *o, NGparam &p);
ReactionClass * createMotorRxn_ccw2cw(MoleculeType * Motor, Observable * cheYp, ofstream *o, NGparam &p);


	///////////////////////////////////////////////////////////////////////
	// Reaction for the Motors depend on observable CheYp
	class MotorSwitch : public ObsDrxn
	{
		public:
			MotorSwitch( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int n_observables,
								Observable ** obs,
								const char * motorRotationStateName,
								int newMotorRotation,
								double Kd,
								double g0,
								double g1,
								double omega,
								double cellVolume,
								ofstream *o);
			virtual ~MotorSwitch();
			
			virtual double update_a();
			virtual void transformReactants(Molecule ** reactants, int nReactants);
			
			static const int CW = 1;
			static const int CCW = 0;
			
		protected:
			int motorRotationIndex;
			int newMotorRotation;
			
			//Saved motor parameters
			double Kd;
			double g1;
			double wFactor;
			double concentrationFactor;
			
			ofstream *o;
	};
	
}


#endif /*MOTOR_HH_*/
