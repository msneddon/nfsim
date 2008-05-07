#include "ng.hh"





using namespace NG;








MoleculeType * NG::makeRecDimer(System * s)
{
	int numOfBsites = 8;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "t0"; //Bonds to link this receptor to others in the hex neighborhood
	bSiteNames[1] = "t1";
	bSiteNames[2] = "t2";
	bSiteNames[3] = "t3";
	bSiteNames[4] = "t4";
	bSiteNames[5] = "t5";
	bSiteNames[6] = "tether";  //Tethering site for CheB or CheR
	bSiteNames[7] = "asite";  //Can bind CheB or CheR at the same site that catalyzes meth / demeth
	
	int numOfStates = 2;
	const char ** stateNames = new const char * [numOfStates];
	stateNames[0] = "m";  //Can be methylated
	stateNames[1] = "type"; //Has a type (ie Tar, Tsr, etc)
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = 2;  //start with methylation level 2 (out of 8 possible)
	stateValues[1] = TAR;  //Assume Tar unless we change it
	
	MoleculeType *receptorDimer = new MoleculeType("ReceptorDimer",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	return receptorDimer;
}



MoleculeType * NG::makeCheA(System * s)
{
	int numOfBsites = 0;
	const char ** bSiteNames = new const char * [numOfBsites];
	
	int  numOfStates = 1;
	const char ** stateNames = new const char * [numOfStates];
	stateNames[0] = "p";  //Can be phosphorylated
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = Param::cheA_phosState;  //Begins unphosphorylated
	
	MoleculeType *cheA = new MoleculeType("CheA",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheA->populateWithDefaultMolecules(Param::cheAcount);
	return cheA;
}


MoleculeType * NG::makeCheR(System *s)
{
	int numOfBsites = 2;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "te";  //Can bind receptor at tethering site
	bSiteNames[1] = "av";  //Can bind receptor at active site
	
	int numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	
	int * stateValues = new int [numOfStates];
	
	MoleculeType *cheR = new MoleculeType("CheR",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheR->populateWithDefaultMolecules(Param::cheRcount);
	cout<<"Created CheR with "<<Param::cheRcount<<" molecules \n";
	return cheR;
}

MoleculeType * NG::makeCheB(System *s)
{
	int numOfBsites = 2;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "te";  //Can bind receptor at tethering site
	bSiteNames[1] = "av";  //Can bind receptor at active site
	
	int numOfStates = 1;
	const char ** stateNames = new const char * [numOfStates];
	stateNames[0] = "p";  //Can be phosphorylated
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = Param::cheB_phosState;
	
	MoleculeType *cheB = new MoleculeType("CheB",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheB->populateWithDefaultMolecules(Param::cheBcount);
	cout<<"Created CheB with "<<Param::cheBcount<<" molecules \n";
	return cheB;
}





MoleculeType * NG::makeCheY(System *s)
{
	int numOfBsites = 0;
	const char ** bSiteNames = new const char * [numOfBsites];
	
	int numOfStates = 1;
	const char ** stateNames = new const char * [numOfStates];
	stateNames[0] = "p";  //Can be phosphorylated
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = Param::cheY_phosState;
	
	MoleculeType *cheY = new MoleculeType("CheY",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheY->populateWithDefaultMolecules(Param::cheYcount);
	cout<<"Created CheY with "<<Param::cheYcount<<" molecules \n";
	return cheY;
}






void NG::makeTarMolecules(MoleculeType * receptorDimer)
{
	for(int i=0; i<Param::tarDimerCount; i++)
	{
		//create a tar, note that we don't have to add it to anything
		//because all moleculeType init is done in the constructor
		Molecule * tar = new Molecule(receptorDimer);
		tar->setState("type",TAR);
		tar->setState("m", Param::tarMethLevel);
	}
}

void NG::makeTsrMolecules(MoleculeType * receptorDimer)
{
	for(int i=0; i<Param::tsrDimerCount; i++)
	{
		//create a tsr, note that we don't have to add it to anything
		//because all moleculeType init is done in the constructor
		Molecule * tsr = new Molecule(receptorDimer);
		tsr->setState("type",TSR);
		tsr->setState("m", Param::tsrMethLevel);
	}
}








////////////////////////////// New Functions ...














