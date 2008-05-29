

#include "ng_creator.hh"






/** Creates Receptor Dimer Molecule Types and Molecule Instances.
 * 
 * 
 */
MoleculeType * createReceptorDimer(System * s, NGparam &p)
{
	//First set up the binding sites we need for the model
	int numOfBsites; string* bSiteNames;
	if(p.get_useNeighborRxns() && p.get_useTether()) {
		numOfBsites = 8; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = "t0";
		bSiteNames[1] = "t1";
		bSiteNames[2] = "t2";
		bSiteNames[3] = "t3";
		bSiteNames[4] = "t4";
		bSiteNames[5] = "t5";	
		bSiteNames[6] = p.get_nameReceptorTetherSite();
		bSiteNames[7] = p.get_nameReceptorActiveSite();
	} else if(!p.get_useNeighborRxns() && p.get_useTether()) {
		numOfBsites = 2; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameReceptorTetherSite();
		bSiteNames[1] = p.get_nameReceptorActiveSite();
	} else if(!p.get_useNeighborRxns() && !p.get_useTether()) {
		numOfBsites = 1; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameReceptorActiveSite();
	} else {
		cerr<<"Receptor Dimers cannot use neighborhood reactions if it is not using tethering reactions!"<<endl;
		cerr<<"Check your parameters (in class NGparam).  Quitting."<<endl;
		exit(1);
	}
	
	//Here we set the states that the receptors can take
	int numOfStates = 2;
	string* stateNames = new string [numOfStates];
	stateNames[0] = p.get_nameReceptorMethState();  //Can be methylated
	stateNames[1] = p.get_nameReceptorTypeState(); //Has a type (ie Tar, Tsr, etc)
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = 2;  //start with methylation level 2 (out of 8 possible)
	stateValues[1] = TAR;  //Assume Tar unless we change it
	
	//Create the actual moleculeType
	MoleculeType *receptorDimer = new MoleculeType("ReceptorDimer",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	
	//Create the actual instances of the Tar dimer
	int level = 0; int levelCount = 0; int loopCount=0;
	for(int r=0; r<p.get_tarCount(); r++, levelCount++, loopCount=0)
	{	
		
		Molecule * tar = receptorDimer->genDefaultMolecule();
		tar->setState(p.get_nameReceptorTypeState(),TAR);
		
		while( levelCount>=p.get_tarMethLevelCount(level) ) {
			levelCount = 0; level++;
			if(level>8) { level=0; loopCount++; }
			if(loopCount>2) {
				cerr<<"You seemed to have set all Tsr Level counts to zero!  I cannot work with that!"<<endl;
				cerr<<"I think I'll quit now."<<endl;
				exit(1);
			}
		}
		tar->setState(p.get_nameReceptorMethState(),level);
	}
	
	//And here create the instances of the Tsr dimer
	level = 0; levelCount = 0; loopCount=0;
	for(int r=0; r<p.get_tsrCount(); r++, levelCount++, loopCount=0)
	{
		Molecule * tsr = receptorDimer->genDefaultMolecule();
		tsr->setState(p.get_nameReceptorTypeState(),TSR);
		
		while( levelCount>=p.get_tsrMethLevelCount(level) ) {
			levelCount = 0; level++;
			if(level>8) { level=0; loopCount++; }
			if(loopCount>2) {
				cerr<<"You seemed to have set all Tsr Level counts to zero!  I cannot work with that!"<<endl;
				cerr<<"I think I'll quit now."<<endl;
				exit(1);
			}
		}
		tsr->setState(p.get_nameReceptorMethState(),level);
	}
	
	return receptorDimer;
}


/** Creates CheA Molecule Types and Molecule Instances.
 * 
 * 
 */
MoleculeType * createCheA(System *s, NGparam &p)
{
	int numOfBsites = 0;
	string * bSiteNames = new string [numOfBsites];
	
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = p.get_nameCheAphosState();  //Can be phosphorylated
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = p.get_cheA_phosState();
	
	MoleculeType *cheA = new MoleculeType("CheA",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheA->populateWithDefaultMolecules(p.get_cheAcount());
	return cheA;
}


/** Creates CheR Molecule Types and Molecule Instances.
 * 
 * 
 */
MoleculeType * createCheR(System *s, NGparam &p)
{
	int numOfBsites; string* bSiteNames;
	if(p.get_useTether()) {
		numOfBsites = 2; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameCheRtetherSite();
		bSiteNames[1] = p.get_nameCheRactiveSite();
	} else {
		numOfBsites = 1; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameCheRactiveSite();
	} 
	
	int numOfStates = 0;
	string * stateNames = new string [numOfStates];
	int * stateValues = new int [numOfStates];
	
	MoleculeType *cheR = new MoleculeType("CheR",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheR->populateWithDefaultMolecules(p.get_cheRcount());
	return cheR;
}


/**  Creates CheB Molecule Types and Molecule Instances.
 * 
 * 
 */
MoleculeType * createCheB(System *s, NGparam &p)
{
	//Set binding sites depending on whether or not we can be tethered
	int numOfBsites; string * bSiteNames;
	if(p.get_useTether()) {
		numOfBsites = 2; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameCheBtetherSite();
		bSiteNames[1] = p.get_nameCheBactiveSite();
	} else {
		numOfBsites = 1; bSiteNames = new string [numOfBsites];
		bSiteNames[0] = p.get_nameCheBactiveSite();
	}
	
	//Set the state that we can be phosphorylated
	int numOfStates; string * stateNames; int * stateValues;
	if(p.get_useCheBFeedback())
	{
		numOfStates = 1;
		stateNames = new string [numOfStates];
		stateNames[0] = p.get_nameCheBphosState();  //Can be phosphorylated
		stateValues = new int [numOfStates];
		stateValues[0] = p.get_cheB_phosState();
	} else {
		numOfStates = 0;
		stateNames = new string [numOfStates];
		stateValues = new int [numOfStates];
	}
	
	MoleculeType *cheB = new MoleculeType("CheB",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheB->populateWithDefaultMolecules(p.get_cheBcount());
	return cheB;
}



/**  Creates CheY Molecule Types and Molecule Instances.
 * 
 * 
 * 
 */
MoleculeType * createCheY(System *s, NGparam &p)
{
	int numOfBsites = 0;
	string * bSiteNames = new string [numOfBsites];
	
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = p.get_nameCheYphosState();  //Can be phosphorylated
	int * stateValues = new int [numOfStates];
	stateValues[0] = p.get_cheY_phosState();
	
	MoleculeType *cheY = new MoleculeType("CheY",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	cheY->populateWithDefaultMolecules(p.get_cheYcount());
	return cheY;
}



/**  Creates Motor Molecule Types and Molecule Instances.
 * 
 * 
 * 
 */
MoleculeType * createMotor(System *s, NGparam &p)
{
	int numOfBsites = 0;
	string * bSiteNames = new string [numOfBsites];
	
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = p.get_nameMotorRotationState();  //The rotation state of the motor (either CW=1 or CCW=0)
	
	int * stateValues = new int [numOfStates];
	stateValues[0] = p.get_motorStartState();
	
	MoleculeType *motor = new MoleculeType("Motor",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	motor->populateWithDefaultMolecules(p.get_motorCount());
	return motor;
}





