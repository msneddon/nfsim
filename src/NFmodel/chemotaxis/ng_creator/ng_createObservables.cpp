#include <sstream>

#include "ng_creator.hh"




void createObservable_ReceptorMethLevel(MoleculeType *recDimer, NGparam &p)
{
	int numberOfSites = p.get_receptorDimerNumberOfMethSites();
	
	
	for (int m=0; m<=numberOfSites; m++)
	{
		TemplateMolecule *rec = new TemplateMolecule(recDimer);
		rec->addStateValue("m",m);
		
//		std::stringstream out;
//		out << "T(m="<<m<<")";
//		std::string s = out.str();
		Observable * obsRec = new Observable("T",rec);
		recDimer->addObservable(obsRec);
	}
}



void createObservable_MotorState(MoleculeType * motor, NGparam &p)
{
	TemplateMolecule *motor_cw = new TemplateMolecule(motor);
	motor_cw->addStateValue(p.get_nameMotorRotationState(),Motor::MotorSwitch::CW);
	Observable * obsMotor_cw = new Observable("Mot(CW)",motor_cw);
	motor->addObservable(obsMotor_cw);
	
	TemplateMolecule *motor_ccw = new TemplateMolecule(motor);
	motor_ccw->addStateValue(p.get_nameMotorRotationState(),Motor::MotorSwitch::CCW);
	Observable * obsMotor_ccw = new Observable("Mot(CCW)",motor_ccw);
	motor->addObservable(obsMotor_ccw);
}