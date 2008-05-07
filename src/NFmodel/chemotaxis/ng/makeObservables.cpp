#include "ng.hh"


using namespace NG;

void NG::addReceptorMethLevelObs(MoleculeType *recDimer)
{
	TemplateMolecule *rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",0);
	Observable * obsRec = new Observable("T(m=0)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",1);
	obsRec = new Observable("T(m=1)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",2);
	obsRec = new Observable("T(m=2)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",3);
	obsRec = new Observable("T(m=3)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",4);
	obsRec = new Observable("T(m=4)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",5);
	obsRec = new Observable("T(m=5)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",6);
	obsRec = new Observable("T(m=6)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",7);
	obsRec = new Observable("T(m=7)",rec);
	recDimer->addObservable(obsRec);
	
	rec = new TemplateMolecule(recDimer);
	rec->addStateValue("m",8);
	obsRec = new Observable("T(m=8)",rec);
	recDimer->addObservable(obsRec);
}



void NG::addCheRB_boundToTether(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer)
{
	TemplateMolecule *tempRecDimer = new TemplateMolecule(recDimer);
	TemplateMolecule *tempCheR = new TemplateMolecule(cheR);
	TemplateMolecule::bind(tempRecDimer, "tether", tempCheR, "te");
	Observable * obs = new Observable("CheR-t-Rec",tempCheR);
	cheR->addObservable(obs);
	
	tempRecDimer = new TemplateMolecule(recDimer);
	TemplateMolecule *tempCheB = new TemplateMolecule(cheB);
	TemplateMolecule::bind(tempRecDimer, "tether", tempCheB, "te");
	obs = new Observable("CheB-t-Rec",tempCheB);
	cheB->addObservable(obs);
}


void NG::addCheRB_boundToTetherAndActiveSite(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer)
{
	TemplateMolecule *tempRecDimer = new TemplateMolecule(recDimer);
	TemplateMolecule *tempCheR = new TemplateMolecule(cheR);
	//TemplateMolecule::bind(tempRecDimer, "tether", tempCheR, "te");
	TemplateMolecule::bind(tempRecDimer, "active", tempCheR, "av");
	Observable * obs = new Observable("CheR-t:a-Rec",tempCheR);
	cheR->addObservable(obs);
	
	tempRecDimer = new TemplateMolecule(recDimer);
	TemplateMolecule *tempCheB = new TemplateMolecule(cheB);
	//TemplateMolecule::bind(tempRecDimer, "tether", tempCheB, "te");
	TemplateMolecule::bind(tempRecDimer, "active", tempCheB, "av");
	obs = new Observable("CheB-t:a-Rec",tempCheB);
	cheB->addObservable(obs);
	
	
	tempRecDimer = new TemplateMolecule(recDimer);
	tempCheB = new TemplateMolecule(cheB);
	TemplateMolecule::bind(tempRecDimer, "active", tempCheB, "av");
	obs = new Observable("CheB-a-Rec",tempCheB);
	cheB->addObservable(obs);
}



void NG::addCheAphos(MoleculeType * cheA)
{
	TemplateMolecule *tempMol = new TemplateMolecule(cheA);
	tempMol->addStateValue("p",PHOS);
	Observable * obs = new Observable("CheA(p)",tempMol);
	cheA->addObservable(obs);
}

void NG::addCheBphos(MoleculeType * cheB)
{
	TemplateMolecule *tempMol = new TemplateMolecule(cheB);
	tempMol->addStateValue("p",PHOS);
	Observable * obs = new Observable("CheB(p)",tempMol);
	cheB->addObservable(obs);
}

void NG::addCheYphos(MoleculeType * cheY)
{
	TemplateMolecule *tempMol = new TemplateMolecule(cheY);
	tempMol->addStateValue("p",PHOS);
	Observable * obs = new Observable("CheY(p)",tempMol);
	cheY->addObservable(obs);
}

void NG::addMotorObs(MoleculeType * motor)
{
	TemplateMolecule *motor_cw = new TemplateMolecule(motor);
	motor_cw->addStateValue("r",1);
	Observable * obsMotor_cw = new Observable("Mot(CW)",motor_cw);
	motor->addObservable(obsMotor_cw);
	
	TemplateMolecule *motor_ccw = new TemplateMolecule(motor);
	motor_ccw->addStateValue("r",0);
	Observable * obsMotor_ccw = new Observable("Mot(CCW)",motor_ccw);
	motor->addObservable(obsMotor_ccw);
}


