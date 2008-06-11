 


#include "ng_param.hh"
#include <iostream>

using namespace std;


NGparam::NGparam()
{
	initEverything();
	initNames();
	setDefaultOutputFileNames();
}

NGparam::~NGparam() {}

/*  A system with all the basic components needed for swimming, including
 * receptor clusters, methylation/demethylation, feedback, and a motor, but
 * excluding tethering and assistance neighborhoods */
void NGparam::setLiteSystem()
{
	initNames();
	
	cellVolume = 1.41e-15; //Liters
	
	//Global Model Parameters
	useTether = false;
	useRetether = false;
	useNeighborRxns = false;
	useCheBFeedback = true;
	useCheA = true;
	useCheY = true;
	useCheB = true;
	useCheR = true;
	useMotor = true;
	useRB_aSiteBindingDependsOnAvail=true;
	
	//Set the growth conditions to rich media
	setRichMediaConcentrations();
	
	//Set energy values and parameters for the receptors / clusters (From Will)
	tarPerCluster = 6;
	tsrPerCluster = 12;
	
	/* will's model
	useWillCluster = true;	
	e0 = 3.36;//
	e1 = -(0.063*18)/(tarPerCluster+tsrPerCluster);
	asp_Koff_TAR = 0.0094e-3;
	asp_Kon_TAR = 0.017e-3;
	asp_Koff_TSR = 150e-3;
	asp_Kon_TSR = 245e-3; */
	
	
	/* wingreen */
	useWillCluster = false;
	asp_Koff_TAR = 0.02e-3; // This is 0.02 mM or 0.00002 M
	asp_Kon_TAR = 0.5e-3; // mM
	asp_Koff_TSR = 100e-3; // mM
	asp_Kon_TSR = 1000; //Yes, 1000 M, or 10^6 mM
	
	
	
	numberOfDimerMethylationSites = 8; // 8 sites in a dimer
	
	
	
	//Starting conditions of the molecules / environment
	aspartateConcentration = 0;
	cheA_phosState = 0;
	cheB_phosState = 1;
	cheY_phosState = 0;
	motorStartState = 0;
	
	tarMethLevel[0] = 3;   tsrMethLevel[0] = 3;
	tarMethLevel[1] = 2;   tsrMethLevel[1] = 2; 
	tarMethLevel[2] = 1;   tsrMethLevel[2] = 1;
	tarMethLevel[3] = 1;   tsrMethLevel[3] = 1;
	tarMethLevel[4] = 0;   tsrMethLevel[4] = 0;
	tarMethLevel[5] = 0;   tsrMethLevel[5] = 0;
	tarMethLevel[6] = 0;   tsrMethLevel[6] = 0;
	tarMethLevel[7] = 0;   tsrMethLevel[7] = 0;
	tarMethLevel[8] = 0;   tsrMethLevel[8] = 0;
	
	
	
	//Methylation Rates from Will's model (scaled for receptor dimers)
	double KM_R_michaelisConstant = 0.2; //uM
	double KM_B_michaelisConstant = 1.8; //uM
	
	//Catylitic rates from will's model
	this->CHER_meth_RECEPTOR = 1.6; //1.35; //13.5; // s^-1
	this->CHEB_demeth_RECEPTOR = 2.16; //2.16; // s^-1
	
	if(!useWillCluster) {
		CHER_meth_RECEPTOR = 1.2; //1.35; //13.5; // s^-1
		CHEB_demeth_RECEPTOR = 3; //2.16; // s^-1
	}
	
	//All we have to do is fix the unbinding rates, and we get the binding rates
	//from will's parameters.
	CHER_unbind_ACTIVE = 30; // per second
	CHEB_unbind_ACTIVE = 30; // per second
	
	FREE_CHER_bind_ACTIVE = (CHER_meth_RECEPTOR + CHER_unbind_ACTIVE) / (KM_R_michaelisConstant);
	FREE_CHER_bind_ACTIVE = FREE_CHER_bind_ACTIVE*1e6 / (cellVolume*NFutil::NA);
	
	FREE_CHEB_bind_ACTIVE = (CHEB_demeth_RECEPTOR + CHEB_unbind_ACTIVE) / (KM_B_michaelisConstant);
	FREE_CHEB_bind_ACTIVE = FREE_CHEB_bind_ACTIVE*1e6 / (cellVolume*NFutil::NA);
	
	
	//Phosphorylation Rates
	AUTO_phos_CHEA = 20;  // s^-1
	
	CHEA_phos_CHEB = 10e6 / (cellVolume*NFutil::NA); // 10 / uM s
	AUTO_dephos_CHEB = 1; // s^-1
	
	CHEA_phos_CHEY = 90e6 / (cellVolume*NFutil::NA); // 90 / uM s;;
	AUTO_dephos_CHEY = 18; // use 11 for the case where binding depends on available sites;  // s^-1
	
	if(!useWillCluster) {
			CHER_meth_RECEPTOR = 1; //1.35; //13.5; // s^-1
			AUTO_dephos_CHEY = 12;
	}
	
	
	//Motor Parameters
	Mot_Kd = 3.06;
	Mot_g1 = 35;

	//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
	Mot_g0 = 6.7;
	Mot_omega = 400;
	
	
	
	
	//Outputters
	outputActivity = false;
	outputTarMethState = false;
	outputTsrMethState = false;
	outputAllMethState = true;
	outputCheAphos = true;
	outputCheBphos = true;
	outputCheYphos = true;
	outputMotorState = true;
	
	
	//Rates that apply to tethering should all be set to zero
	//as they are not used
	FREE_CHER_bind_TETHER = 0;
	FREE_CHEB_bind_TETHER = 0;	
	CHER_unbind_TETHER = 0;
	CHEB_unbind_TETHER = 0;
	TETHERED_CHER_bind_ACTIVE = 0;
	TETHERED_CHER_bind_NEIGHBOR_ACTIVE = 0;
	TETHERED_CHEB_bind_ACTIVE = 0;
	TETHERED_CHEB_bind_NEIGHBOR_ACTIVE = 0;
	ACTIVE_BOUND_CHER_bind_TETHER = 0;
	ACTIVE_BOUND_CHEB_bind_TETHER = 0;
}



void NGparam::setFullLite()
{
	
	cellVolume = 1.41e-15; //Liters
		
		
	//Global Model Parameters
	useTether = false;
	useRetether = false;
	useNeighborRxns = false;
	useCheBFeedback = true;
	useCheA = true;
	useCheY = true;
	useCheB = true;
	useCheR = true;
	useMotor = false;
		
		
	//Molecule Counts
	tarCount = 100;
	tsrCount = 100;
	cheAcount = 50;
	cheRcount = 140;
	cheBcount = 240;
	cheYcount = 8000;
	motorCount = 5;
		
		
	//Outputters
	outputActivity = true;
	outputTarMethState = false;
	outputTsrMethState = false;
	outputAllMethState = true;
	outputCheAphos = true;
	outputCheBphos = false;
	outputCheYphos = true;
	outputMotorState = false;
		
		
	//Starting conditions of the molecules / environment
	aspartateConcentration = 0;
	cheA_phosState = 0;
	cheB_phosState = 1;
	cheY_phosState = 1;
	motorStartState = 0;
	
	tarMethLevel[0] = 3;   tsrMethLevel[0] = 3;
	tarMethLevel[1] = 2;   tsrMethLevel[1] = 2; 
	tarMethLevel[2] = 1;   tsrMethLevel[2] = 1;
	tarMethLevel[3] = 1;   tsrMethLevel[3] = 1;
	tarMethLevel[4] = 0;   tsrMethLevel[4] = 0;
	tarMethLevel[5] = 0;   tsrMethLevel[5] = 0;
	tarMethLevel[6] = 0;   tsrMethLevel[6] = 0;
	tarMethLevel[7] = 0;   tsrMethLevel[7] = 0;
	tarMethLevel[8] = 0;   tsrMethLevel[8] = 0;
		
	numberOfDimerMethylationSites = 10;
	asp_Koff_TAR = 0.0094e-3;
	asp_Kon_TAR = 0.17e-3;
	asp_Koff_TSR = 150e-3;
	asp_Kon_TSR = 245e-3;
	
		

	//Receptor Cluster Parameters
	tarPerCluster=1;
	tsrPerCluster=1;
	
	e0=0;
	e1=0;
	
	TAR_freeEnergyOffset[0] = 1;		TSR_freeEnergyOffset[0] = 1;
	TAR_freeEnergyOffset[1] = 2;		TSR_freeEnergyOffset[2] = 3;
	TAR_freeEnergyOffset[3] = 2;		TSR_freeEnergyOffset[3] = 2;
	TAR_freeEnergyOffset[4] = 1;		TSR_freeEnergyOffset[4] = 1;
	TAR_freeEnergyOffset[5] = 0;		TSR_freeEnergyOffset[5] = 0;
	TAR_freeEnergyOffset[6] = 0;		TSR_freeEnergyOffset[6] = 0;
	TAR_freeEnergyOffset[7] = 0;		TSR_freeEnergyOffset[7] = 0;
	TAR_freeEnergyOffset[8] = 0;		TSR_freeEnergyOffset[8] = 0;
	
	
/*	double asp_Koff_TAR;
	double asp_Kon_TAR;
	double asp_Koff_TSR;
	double asp_Kon_TSR;
		
		
		
		//Rate Constants
		
		//Tethering Rates
		double FREE_CHER_bind_TETHER;
		double FREE_CHEB_bind_TETHER;	
		double CHER_unbind_TETHER;
		double CHEB_unbind_TETHER;
		
		double TETHERED_CHER_bind_ACTIVE;
		double TETHERED_CHER_bind_NEIGHBOR_ACTIVE;
		double TETHERED_CHEB_bind_ACTIVE;
		double TETHERED_CHEB_bind_NEIGHBOR_ACTIVE;

		//Active Site Rates
		double FREE_CHER_bind_ACTIVE;
		double FREE_CHEB_bind_ACTIVE;
		double CHER_unbind_ACTIVE;
		double CHEB_unbind_ACTIVE;
		
		double ACTIVE_BOUND_CHER_bind_TETHER;
		double ACTIVE_BOUND_CHEB_bind_TETHER;
		
		//Meth / Demeth Rates
		static double CHER_meth_RECEPTOR;   // per second
		static double CHEB_demeth_RECEPTOR; // per second
		
		//Phosphorylation Rates
		static double AUTO_phos_CHEA;
		static double CHEA_phos_CHEY;
		static double CHEA_phos_CHEB;
		static double AUTO_dephos_CHEY;
		static double AUTO_dephos_CHEB;
		
		
		
		//Motor Parameters
		double Mot_Kd;
		double Mot_g0;
		
		//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
		static double Mot_g1;
		static double Mot_omega;*/
		
		
	
		setWillParameters();
}


void NGparam::setCheR_2x() { cheRcount*=2; }
void NGparam::setCheR_4x() { cheRcount*=4; }
void NGparam::setCheR_8x() { cheRcount*=8; }


void NGparam::setRichMediaConcentrations() {
	tarCount = 4500/2; //Divide by two because we use dimers
	tsrCount = 9000/2; 
	cheAcount = 4500;
	cheRcount = 140;
	cheBcount = 240;
	cheYcount = 8200;
	//AUTO_dephos_CHEY = proportional to cheZ (3200)
	motorCount = 1;
}

void NGparam::setMinimalMediaConcentrations() {
	tarCount = 8000/2; //Divide by two because we use dimers
	tsrCount = 16000/2;
	cheAcount = 5100;
	cheRcount = 160;
	cheBcount = 270;
	cheYcount = 6300;
	//AUTO_dephos_CHEY = proportional to cheZ (2700)
	motorCount = 1;
}

void NGparam::setInitReceptorMethToZero() {
	tarMethLevel[0] = 1;   tsrMethLevel[0] = 1;
	tarMethLevel[1] = 0;   tsrMethLevel[1] = 0; 
	tarMethLevel[2] = 0;   tsrMethLevel[2] = 0;
	tarMethLevel[3] = 0;   tsrMethLevel[3] = 0;
	tarMethLevel[4] = 0;   tsrMethLevel[4] = 0;
	tarMethLevel[5] = 0;   tsrMethLevel[5] = 0;
	tarMethLevel[6] = 0;   tsrMethLevel[6] = 0;
	tarMethLevel[7] = 0;   tsrMethLevel[7] = 0;
	tarMethLevel[8] = 0;   tsrMethLevel[8] = 0;
}
void NGparam::setInitReceptorMethToEight() {
	tarMethLevel[0] = 0;   tsrMethLevel[0] = 0;
	tarMethLevel[1] = 0;   tsrMethLevel[1] = 0; 
	tarMethLevel[2] = 0;   tsrMethLevel[2] = 0;
	tarMethLevel[3] = 0;   tsrMethLevel[3] = 0;
	tarMethLevel[4] = 0;   tsrMethLevel[4] = 0;
	tarMethLevel[5] = 0;   tsrMethLevel[5] = 0;
	tarMethLevel[6] = 0;   tsrMethLevel[6] = 0;
	tarMethLevel[7] = 0;   tsrMethLevel[7] = 0;
	tarMethLevel[8] = 1;   tsrMethLevel[8] = 1;
}


void NGparam::setToNoiseOutput() {
	moleculeOutputFileName = "/home/msneddon/Desktop/will_test/noise/NG_basicOutput.txt";
	keyFileName = "/home/msneddon/Desktop/will_test/noise/NG_groupKey.txt";
	activityFileName = "/home/msneddon/Desktop/will_test/noise/NG_pOnVals.txt";
	motorOutputFileName = "/home/msneddon/Desktop/will_test/noise/NG_motorOutput.txt";
	
	outputActivity = false;
	outputTarMethState = false;
	outputTsrMethState = false;
	outputAllMethState = false;
	outputCheAphos = false;
	outputCheBphos = false;
	outputCheYphos = true;
	outputMotorState = false;
}

void NGparam::setToActivityOutput() {
	outputActivity = true;
	outputTarMethState = false;
	outputTsrMethState = false;
	outputAllMethState = true;
	outputCheAphos = true;
	outputCheBphos = true;
	outputCheYphos = true;
	outputMotorState = false;
}


void NGparam::setWillParameters() {
	
	//Set concentrations (roughly from Li and Hazelbauer)
	tarCount = 4500; //4500/10;  
	tsrCount = 9000; //9000/10;
	cheAcount = 4500; //4500/10;
	cheRcount =  255; //85;
	cheBcount = 255; //509;
	cheYcount = 8236; //8200/10;
	motorCount = 0;
	
	//15.9
	//this->useCheB=false;
	//this->setInitReceptorMethToZero();
	
	//setInitReceptorMethToZero(); 
	//number of clusters = 750
	this->tarPerCluster = 6; //9; //(6*75)/10;
	this->tsrPerCluster = 12; //18; //(12*75)/10;
	numberOfDimerMethylationSites = 4;
	tarMethLevel[0] = 0;   tsrMethLevel[0] = 0;
	tarMethLevel[1] = 1;   tsrMethLevel[1] = 1; 
	tarMethLevel[2] = 2;   tsrMethLevel[2] = 2;
	tarMethLevel[3] = 1;   tsrMethLevel[3] = 1;
	tarMethLevel[4] = 0;   tsrMethLevel[4] = 0;
	tarMethLevel[5] = 0;   tsrMethLevel[5] = 0;
	tarMethLevel[6] = 0;   tsrMethLevel[6] = 0;
	tarMethLevel[7] = 0;   tsrMethLevel[7] = 0;
	tarMethLevel[8] = 0;   tsrMethLevel[8] = 0;
	
	
	//Set offset energy for receptor clusters
	e0 = 3.36*1;// / 18;  //Divide by 18 because will has 18 dimers per cluster
	e1 = -(0.063); //*36)/(tarPerCluster+tsrPerCluster);
/*	for(int m=0; m<=8; m++) {
		TAR_freeEnergyOffset[m] = e0+e1*m;
		TSR_freeEnergyOffset[m] = e0+e1*m;
	}*/
	
	
	//Ligand binding and unbinding to tar and tsr
	asp_Koff_TAR = 0.0094e-3;
	asp_Kon_TAR = 0.017e-3;
	asp_Koff_TSR = 150e-3;
	asp_Kon_TSR = 245e-3;
	
	
	//Methylation Rates from Will's model (scaled for receptor dimers)
	double KM_R_michaelisConstant = 0.2; //uM
	double KM_B_michaelisConstant = 1; //uM
	
	//Catylitic rates from will's model
	this->CHER_meth_RECEPTOR = 0.6; //1.35; //13.5; // s^-1
	this->CHEB_demeth_RECEPTOR = 0.8; //2.16; // s^-1
	
	//All we have to do is fix the unbinding rates, and we get the binding rates
	//from will's parameters.
	CHER_unbind_ACTIVE = 50; // per second
	CHEB_unbind_ACTIVE = 50; // per second
	
	FREE_CHER_bind_ACTIVE = (CHER_unbind_ACTIVE) / (KM_R_michaelisConstant);
	FREE_CHER_bind_ACTIVE = FREE_CHER_bind_ACTIVE*1e6 / (cellVolume*NFutil::NA);
	
	FREE_CHEB_bind_ACTIVE = (CHEB_unbind_ACTIVE) / (KM_B_michaelisConstant);
	FREE_CHEB_bind_ACTIVE = FREE_CHEB_bind_ACTIVE*1e6 / (cellVolume*NFutil::NA);
	
	
	//Phosphorylation Rates
	AUTO_phos_CHEA = 20;  // s^-1
	
	CHEA_phos_CHEB = 10e6 / (cellVolume*NFutil::NA); // 10 / uM s
	AUTO_dephos_CHEB = 1; // s^-1
	
	CHEA_phos_CHEY = 80e6 / (cellVolume*NFutil::NA); // 90 / uM s;;
	AUTO_dephos_CHEY = 15;  // s^-1
	
	
	
	//All other rates (ie tethering or whatever) should be set to zero
	//As they are not a part of will's model
	FREE_CHER_bind_TETHER = 0;
	FREE_CHEB_bind_TETHER = 0;	
	CHER_unbind_TETHER = 0;
	CHEB_unbind_TETHER = 0;
	TETHERED_CHER_bind_ACTIVE = 0;
	TETHERED_CHER_bind_NEIGHBOR_ACTIVE = 0;
	TETHERED_CHEB_bind_ACTIVE = 0;
	TETHERED_CHEB_bind_NEIGHBOR_ACTIVE = 0;
	ACTIVE_BOUND_CHER_bind_TETHER = 0;
	ACTIVE_BOUND_CHEB_bind_TETHER = 0;
}

















void NGparam::setFullSystem()
{
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}







/**
 * 
 * 
 */
void NGparam::initNames() 
{
	receptorTetherSiteName = "tether";
	receptorActiveSiteName = "active";
	receptorMethStateName = "m";
	receptorTypeStateName = "type";
	
	cheRtetherSiteName = "te";
	cheBtetherSiteName = "te";
	cheRactiveSiteName = "av";
	cheBactiveSiteName = "av";
	
	cheAphosSiteName = "p";
	cheBphosSiteName = "p";
	cheYphosSiteName = "p";
	motorRotationStateName = "ro";
}


/**
 * 
 * 
 */
void NGparam::setDefaultOutputFileNames() 
{
	moleculeOutputFileName = "/home/msneddon/Desktop/will_test/NG_basicOutput.txt";
	keyFileName = "/home/msneddon/Desktop/will_test/NG_groupKey.txt";
	activityFileName = "/home/msneddon/Desktop/will_test/NG_pOnVals.txt";
	motorOutputFileName = "/home/msneddon/Desktop/will_test/NG_motorOutput.txt";
}

void NGparam::setOutputFileNames(const char * moleculeOutput, const char * key, const char * activity, const char * motor) 
{
	moleculeOutputFileName = moleculeOutput;
	keyFileName = key;
	activityFileName = activity;
	motorOutputFileName = motor;
}



/**
 * Blindly initialize everything to a zero or useless value.  This is always called in the constructor to
 * make sure every paramater starts with some value.
 * 
 */
void NGparam::initEverything()
{
	moleculeOutputFileName="";
	keyFileName="";
	activityFileName="";
	motorOutputFileName="";
	outputActivity=false;
	outputTarMethState=false;
	outputTsrMethState=false;
	outputAllMethState=false;
	outputCheAphos=false;
	outputCheBphos=false;
	outputCheYphos=false;
	outputMotorState=false;
	receptorTetherSiteName="";
	receptorActiveSiteName="";
	receptorMethStateName="";
	receptorTypeStateName="";
	cheRtetherSiteName="";
	cheBtetherSiteName="";
	cheRactiveSiteName="";
	cheBactiveSiteName="";
	cheAphosSiteName="";
	cheBphosSiteName="";
	cheYphosSiteName="";
	motorRotationStateName="";
	cellVolume=1;
	useTether=false;
	useRetether=false;
	useNeighborRxns=false;
	useCheBFeedback=false;
	useCheA=false;
	useCheY=false;
	useCheB=false;
	useCheR=false;
	useMotor=false;
	useWillCluster = true;
	useRB_aSiteBindingDependsOnAvail=false;
	tarCount=0;
	tsrCount=0;
	cheAcount=0;
	cheRcount=0;
	cheBcount=0;
	cheYcount=0;
	motorCount=0;
	aspartateConcentration=0;
	cheA_phosState=0;
	cheB_phosState=0;
	cheY_phosState=0;
	motorStartState=0;
	for(int m=0; m<=8; m++) {
		tarMethLevel[m]=0;
		tsrMethLevel[m]=0;
		TAR_freeEnergyOffset[m] = 0;
		TSR_freeEnergyOffset[m] = 0;
	}
	tarPerCluster=0;
	tsrPerCluster=0;
	e0=0;
	e1=0;
	asp_Koff_TAR=1;
	asp_Kon_TAR=1;
	asp_Koff_TSR=1;
	asp_Kon_TSR=1;
	numberOfDimerMethylationSites=1;
	
	//Rate Constants
	FREE_CHER_bind_TETHER=0;
	FREE_CHEB_bind_TETHER=0;	
	CHER_unbind_TETHER=0;
	CHEB_unbind_TETHER=0;
	TETHERED_CHER_bind_ACTIVE=0;
	TETHERED_CHER_bind_NEIGHBOR_ACTIVE=0;
	TETHERED_CHEB_bind_ACTIVE=0;
	TETHERED_CHEB_bind_NEIGHBOR_ACTIVE=0;
	FREE_CHER_bind_ACTIVE=0;
	FREE_CHEB_bind_ACTIVE=0;
	CHER_unbind_ACTIVE=0;
	CHEB_unbind_ACTIVE=0;
	ACTIVE_BOUND_CHER_bind_TETHER=0;
	ACTIVE_BOUND_CHEB_bind_TETHER=0;
	CHER_meth_RECEPTOR=0;
	CHEB_demeth_RECEPTOR=0;
	AUTO_phos_CHEA=0;
	CHEA_phos_CHEY=0;
	CHEA_phos_CHEB=0;
	AUTO_dephos_CHEY=0;
	AUTO_dephos_CHEB=0;
	Mot_Kd=1;
	Mot_g0=1;
	
	//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
	Mot_g1=1;
	Mot_omega=1;
	
}





