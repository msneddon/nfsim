//#include <sstream>
//
//#include "ng_creator.hh"
//
//
//
//
//void createObservable_ReceptorMethLevel(MoleculeType *recDimer, NGparam &p)
//{
//	int numberOfSites = p.get_receptorDimerNumberOfMethSites();
//	
//	
//	for (int m=0; m<=numberOfSites; m++)
//	{
//		TemplateMolecule *rec = new TemplateMolecule(recDimer);
//		rec->addStateValue("m",m);
//		
////		std::stringstream out;
////		out << "T(m="<<m<<")";
////		std::string s = out.str();
//		Observable * obsRec = new Observable("T",rec);
//		recDimer->addObservable(obsRec);
//	}
//}
//
//
//
//void createObservable_MotorState(MoleculeType * motor, NGparam &p)
//{
//	TemplateMolecule *motor_cw = new TemplateMolecule(motor);
//	motor_cw->addStateValue(p.get_nameMotorRotationState(),Motor::MotorSwitch::CW);
//	Observable * obsMotor_cw = new Observable("Mot(CW)",motor_cw);
//	motor->addObservable(obsMotor_cw);
//	
//	TemplateMolecule *motor_ccw = new TemplateMolecule(motor);
//	motor_ccw->addStateValue(p.get_nameMotorRotationState(),Motor::MotorSwitch::CCW);
//	Observable * obsMotor_ccw = new Observable("Mot(CCW)",motor_ccw);
//	motor->addObservable(obsMotor_ccw);
//}
//
//
//
//
//
//
//void createObservable_CheRB_boundToTether(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer, NGparam &p)
//{
//	TemplateMolecule *tempRecDimer = new TemplateMolecule(recDimer);
//	TemplateMolecule *tempCheR = new TemplateMolecule(cheR);
//	TemplateMolecule::bind(tempRecDimer, "tether", tempCheR, "te");
//	Observable * obs = new Observable("CheR-t-Rec",tempCheR);
//	cheR->addObservable(obs);
//	
//	tempRecDimer = new TemplateMolecule(recDimer);
//	TemplateMolecule *tempCheB = new TemplateMolecule(cheB);
//	TemplateMolecule::bind(tempRecDimer, "tether", tempCheB, "te");
//	obs = new Observable("CheB-t-Rec",tempCheB);
//	cheB->addObservable(obs);
//}
//
//
//void createObservable_CheRB_boundToTetherAndActiveSite(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer, NGparam &p)
//{
//	TemplateMolecule *tempRecDimer = new TemplateMolecule(recDimer);
//	TemplateMolecule *tempCheR = new TemplateMolecule(cheR);
//	//TemplateMolecule::bind(tempRecDimer, "tether", tempCheR, "te");
//	TemplateMolecule::bind(tempRecDimer, "active", tempCheR, "av");
//	Observable * obs = new Observable("CheR-t:a-Rec",tempCheR);
//	cheR->addObservable(obs);
//	
//	tempRecDimer = new TemplateMolecule(recDimer);
//	TemplateMolecule *tempCheB = new TemplateMolecule(cheB);
//	//TemplateMolecule::bind(tempRecDimer, "tether", tempCheB, "te");
//	TemplateMolecule::bind(tempRecDimer, "active", tempCheB, "av");
//	obs = new Observable("CheB-t:a-Rec",tempCheB);
//	cheB->addObservable(obs);
//	
//	
//	tempRecDimer = new TemplateMolecule(recDimer);
//	tempCheB = new TemplateMolecule(cheB);
//	TemplateMolecule::bind(tempRecDimer, "active", tempCheB, "av");
//	obs = new Observable("CheB-a-Rec",tempCheB);
//	cheB->addObservable(obs);
//}
//
//
//
//
//
//
//void createObservable_CheAphos(MoleculeType * cheA, NGparam &p)
//{
//	TemplateMolecule *tempMol = new TemplateMolecule(cheA);
//	tempMol->addStateValue("p",PHOS);
//	Observable * obs = new Observable("CheA(p)",tempMol);
//	cheA->addObservable(obs);
//}
//
//void createObservable_CheBphos(MoleculeType * cheB, NGparam &p)
//{
//	TemplateMolecule *tempMol = new TemplateMolecule(cheB);
//	tempMol->addStateValue("p",PHOS);
//	Observable * obs = new Observable("CheB(p)",tempMol);
//	cheB->addObservable(obs);
//}
//
//void createObservable_CheYphos(MoleculeType * cheY, NGparam &p)
//{
//	TemplateMolecule *tempMol = new TemplateMolecule(cheY);
//	tempMol->addStateValue("p",PHOS);
//	Observable * obs = new Observable("CheY(p)",tempMol);
//	cheY->addObservable(obs);
//}
//
//void createObservable_MotorObs(MoleculeType * motor, NGparam &p)
//{
//	TemplateMolecule *motor_cw = new TemplateMolecule(motor);
//	motor_cw->addStateValue("r",1);
//	Observable * obsMotor_cw = new Observable("Mot(CW)",motor_cw);
//	motor->addObservable(obsMotor_cw);
//	
//	TemplateMolecule *motor_ccw = new TemplateMolecule(motor);
//	motor_ccw->addStateValue("r",0);
//	Observable * obsMotor_ccw = new Observable("Mot(CCW)",motor_ccw);
//	motor->addObservable(obsMotor_ccw);
//}