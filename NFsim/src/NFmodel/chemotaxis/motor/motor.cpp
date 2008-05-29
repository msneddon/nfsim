#include "motor.hh"

#include <math.h>


using namespace Motor;
//
//
//MoleculeType * Motor::makeMotor(System *s)
//{
//	int numOfBsites = 0;
//	const char ** bSiteNames = new const char * [numOfBsites];
//	
//	int numOfStates = 1;
//	const char ** stateNames = new const char * [numOfStates];
//	stateNames[0] = "r";  //The rotation state of the motor (either CW=1 or CCW=0)
//	
//	int * stateValues = new int [numOfStates];
//	stateValues[0] = Param::motorStartState;  //start not swimming 
//	
//	MoleculeType *motor = new MoleculeType("Motor",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
//	motor->populateWithDefaultMolecules(Param::motorCount);
//	return motor;
//}
//
//
//
//
//void Motor::ANrxn_addMotorSwitching_Full(System *s, MoleculeType *motor, MoleculeType *cheY, const char * motorStatsOutputFilename)
//{
//	//For now, we create another observable - but it may prove faster to pass
//	//a cheYp observable that we used before
//	TemplateMolecule *Yp = new TemplateMolecule(cheY);
//	Yp->addStateValue("p",PHOS);
//	Observable * obs_cheYp = new Observable("Yp(motor)",Yp);
//	cheY->addObservable(obs_cheYp);
//	
//	//Here, open a file to write
//	if(NG::Param::outputMotorStats)
//	{
//		ofstream *o = new ofstream();
//		o->open(motorStatsOutputFilename);
//		o->setf(ios::scientific);
//		(*o)<<"#\tcdtime\tmotor\tevent"<<endl;
//		s->addReaction(createMotor_cw2ccw(motor, obs_cheYp, o));
//		s->addReaction(createMotor_ccw2cw(motor, obs_cheYp, o));
//	}
//	else
//	{
//		s->addReaction(createMotor_cw2ccw(motor, obs_cheYp, NULL));
//		s->addReaction(createMotor_ccw2cw(motor, obs_cheYp, NULL));
//		
//	}
//	
//	
//}
//
//
//void Motor::createRxn_MotorSwitching(System *s, MoleculeType *motor, MoleculeType *cheY, NGparam &p)
//{
//	//For now, we create another observable - but it may prove faster to pass
//	//a cheYp observable that we used before
//	TemplateMolecule *Yp = new TemplateMolecule(cheY);
//	Yp->addStateValue(p.get_nameCheYphosState(),PHOS);
//	Observable * obs_cheYp = new Observable("Yp(motor)",Yp);
//	cheY->addObservable(obs_cheYp);
//		
//	//Here, open a file to write
//	if(p.get_outputMotorState())
//	{
//		ofstream *o = new ofstream();
//		o->open(p.get_fileNameMotorOutput());
//		o->setf(ios::scientific);
//		(*o)<<"#\tcdtime\tmotor\tevent"<<endl;
//		s->addReaction(createMotorRxn_cw2ccw(motor, obs_cheYp, o, p));
//		s->addReaction(createMotorRxn_ccw2cw(motor, obs_cheYp, o, p));
//	}
//	else
//	{
//		s->addReaction(createMotorRxn_cw2ccw(motor, obs_cheYp, NULL, p));
//		s->addReaction(createMotorRxn_ccw2cw(motor, obs_cheYp, NULL, p));
//	}
//}
//
//
//
//
//ReactionClass * Motor::createMotorRxn_cw2ccw(MoleculeType * Motor, Observable * cheYp, ofstream *o, NGparam &p)
//{
//	int n_reactants = 1;
//	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
//	
//	reactantTemplates[0] = new TemplateMolecule(Motor);
//	reactantTemplates[0]->addStateValue(p.get_nameMotorRotationState(),MotorSwitch::CW);  //Set constraint that we are in CW state
//	
//	int n_observables = 1;
//	Observable ** obs = new Observable *[n_observables];
//	obs[0] = cheYp;
//	
//	return new MotorSwitch("Motor(CW -> CCW)",n_reactants, reactantTemplates, n_observables, obs, p.get_nameMotorRotationState(), 
//		MotorSwitch::CCW, p.get_motKd(), p.get_motg0(), p.get_motg1(), p.get_motOmega(), p.get_cellVolume(),o);
//}
//
//ReactionClass * Motor::createMotorRxn_ccw2cw(MoleculeType * Motor, Observable * cheYp, ofstream *o, NGparam &p)
//{
//	int n_reactants = 1;
//	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
//	
//	reactantTemplates[0] = new TemplateMolecule(Motor);
//	reactantTemplates[0]->addStateValue(p.get_nameMotorRotationState(),MotorSwitch::CCW);  //Set constraint that we are in CW state
//	
//	int n_observables = 1;
//	Observable ** obs = new Observable *[n_observables];
//	obs[0] = cheYp;
//	
//	return new MotorSwitch("Motor(CCW -> CW)",n_reactants, reactantTemplates, n_observables, obs, p.get_nameMotorRotationState(), 
//		MotorSwitch::CW, p.get_motKd(), p.get_motg0(), p.get_motg1(), p.get_motOmega(), p.get_cellVolume(),o);
//}
//
//
//
//
//ReactionClass * Motor::createMotor_cw2ccw(MoleculeType * Motor, Observable * cheYp, ofstream *o)
//{
//	int n_reactants = 1;
//	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
//	
//	reactantTemplates[0] = new TemplateMolecule(Motor);
//	reactantTemplates[0]->addStateValue("r",MotorSwitch::CW);  //Set constraint that we are in CW state
//	
//	int n_observables = 1;
//	Observable ** obs = new Observable *[n_observables];
//	obs[0] = cheYp;
//	
//	return new MotorSwitch("Motor(CW -> CCW)",n_reactants, reactantTemplates, n_observables, obs, "r", 
//		MotorSwitch::CCW, Param::Mot_Kd, Param::Mot_g0, Param::Mot_g1, Param::Mot_omega, Param::CELL_VOLUME,o);
//}
//
//ReactionClass * Motor::createMotor_ccw2cw(MoleculeType * Motor, Observable * cheYp, ofstream *o)
//{
//	int n_reactants = 1;
//	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
//	
//	reactantTemplates[0] = new TemplateMolecule(Motor);
//	reactantTemplates[0]->addStateValue("r",MotorSwitch::CCW);  //Set constraint that we are in CW state
//	
//	int n_observables = 1;
//	Observable ** obs = new Observable *[n_observables];
//	obs[0] = cheYp;
//	
//	return new MotorSwitch("Motor(CCW -> CW)",n_reactants, reactantTemplates, n_observables, obs, "r", 
//		MotorSwitch::CW, Param::Mot_Kd, Param::Mot_g0, Param::Mot_g1, Param::Mot_omega, Param::CELL_VOLUME,o);
//}
//
//
//
////////////////////////////  Motor Switch reactions
//MotorSwitch::MotorSwitch( char * name, 
//							int n_reactants, 
//							TemplateMolecule ** reactantTemplates,
//							int n_observables,
//							Observable ** obs,
//							const char * motorRotationStateName,
//							int newMotorRotation,
//							double Kd,
//							double g0,
//							double g1,
//							double omega,
//							double cellVolume,
//							ofstream *o) : ObsDrxn(name, n_reactants, reactantTemplates, n_observables, obs)
//{
//	//A reference for parameters
//	//this->Kd = 3.06;
//	//double g0 = 6.7;
//	//this->g1 = 80;
//	//double omega = 400;
//	//	~eColiVolume = 1.41e-15; Liters
//	
//	this->Kd = Kd;
//	this->g1 = g1;
//	//w = ( (omega*g0*sqrt(32)) / (2*pi) ) * exp(-g0);
//	//this->wFactor = ( omega*g0*sqrt(32.0) / (2*M_PI) ) * exp(-g0);
//	this->wFactor = 1.02;
//	
//	double avogadro = 6.02214179e23; // molecules per mole
//	// 1uM = 10^-6 mole / L
//	this->concentrationFactor = 1.0 / (cellVolume*avogadro*1e-6); // now in units uM / mole
//	// then we just multiply the number of molecules by this factor to get the concentration of molecules
//	
//	this->motorRotationIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(motorRotationStateName);
//	this->newMotorRotation = newMotorRotation;
//	this->o = o;
//}							
//
//MotorSwitch::~MotorSwitch() 
//{
//	///////////////////////////////
//	///  WARNING!!! THIS IS CLEARLY A HACK!!
//	///  there should be a better way to output this info to the file
//	///  but this way is the easiest to code for now.  Update this!!
//	///  (we want a better way to close out this stream!)
//	if(newMotorRotation==1)
//	{
//		//cout<<"here"<<endl;
//		//cout<<o<<endl;
//		if(o!=NULL) {
//		o->flush(); 
//		o->close(); 
//		delete o;
//		}
//	}
//};
//
//
//double MotorSwitch::update_a()
//{
//	//Return zero if we don't have anything
//	//if(getReactantCount(0)==0) { a=0; rate=0; return 0; }
//	
//	//Get CheY value
//	int YpNum = obs[0]->getCount();
//	double Yp = concentrationFactor*(double)YpNum;
//	//cout<<Yp<<endl;
//	
//	//Calculate free energy and rate
//	double deltaG = (g1/2.0) * ((1.0/2.0) - (Yp / (Kd + Yp)) );
//	if(newMotorRotation==CW)
//	{
//		//This is the rate for CCW->CW
//		rate = wFactor * exp(-deltaG);	
//	}
//	else if(newMotorRotation==CCW)
//	{
//		//This is the rate for CW->CCW
//		rate = wFactor * exp(deltaG);
//	}
//	else
//	{
//		cerr<<"!!! Error - invalid rotation in Motor reaction definition (called from MotorSwitch)"<<endl;	
//	}
//	
//	
//	//Only one reactant (the motor) so we can just multiply by the reactant count
//	a = rate*getReactantCount(0);
//	return a;
//}
//void MotorSwitch::transformReactants(Molecule ** reactants, int nReactants)
//{
//	if(o!=NULL)
//	{
//	(*o)<<"\t"<<reactants[0]->getMoleculeType()->getSystem()->getCurrentTime()<<"\t";
//	(*o)<<reactants[0]->getMoleculeID()<<"\t"<<newMotorRotation<<endl;
//	}
//	reactants[0]->setState(motorRotationIndex, newMotorRotation);
//}
