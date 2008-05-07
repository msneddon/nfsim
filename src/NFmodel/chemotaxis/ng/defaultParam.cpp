#include "ng.hh"


using namespace NG;


//////////////////////////////
/////Default Parameters

//Main simulation parameters
double Param::aspConc = 0;
bool Param::outputGroupStats = false;
bool Param::outputMotorStats = false;
bool Param::useFullSystem = true;
bool Param::useTetherRxns = true;
bool Param::useCheBFeedbackLoop = true;

bool Param::useCheAPhosObs = false;
bool Param::useCheBPhosObs = false;
bool Param::useMotorObs = true;
bool Param::useMethLevelObs = false;
bool Param::useRB_boundTetherObs = false;
bool Param::useRB_boundActiveSiteObs = false;



////////////////////////////////////////////////////////////////
//Initial counts of molecules and size of clusters


// /*Rich media (from li and hazelbauer, 2004)
int cfactor = 1;

int Param::receptorClusterCount = 395/cfactor;
int Param::tarDimerPerClusterCount = 6*cfactor;
int Param::tsrDimerPerClusterCount = 13*cfactor;
int Param::cheAPerClusterCount = 17*cfactor;

int Param::cheRcount = 140;//140;
int Param::cheBcount = 240;
int Param::cheYcount = 8200;
int Param::motorCount = 1; 
// */


/* Minimal media
int Param::receptorClusterCount = 700;
int Param::tarDimerPerClusterCount = 6;
int Param::tsrDimerPerClusterCount = 13;
int Param::cheAPerClusterCount = 11;

int Param::cheRcount = 160;
int Param::cheBcount = 270;
int Param::cheYcount = 6300;
int Param::motorCount = 10;
*/



//Starting conditions of the molecules
int Param::cheB_phosState = PHOS;
int Param::cheA_phosState = PHOS;
int Param::cheY_phosState = UNPHOS;
int Param::motorStartState = Motor::MotorSwitch::CCW;
int Param::tarMethLevel = 2; 
int Param::tsrMethLevel = 2;


////////////////////////////////////////////////////////////////
//Rate constants

//CheR & CheB tethering rxns
//double Param::FREE_CHER_bind_TETHER = 0.94e6/(Param::NA*Param::CELL_VOLUME);//0.01e6/(NA*CELL_VOLUME);
//double Param::FREE_CHEB_bind_TETHER = 7.4e6/(Param::NA*Param::CELL_VOLUME); //
//double Param::TETHERED_CHER_bind_ACTIVE = 10; // s^-1
//double Param::TETHERED_CHEB_bind_ACTIVE = 10; // s^-1
//double Param::CHER_unbind_TETHER = 0.2;//0.2; // s^-1
//double Param::CHEB_unbind_TETHER = 0.2; // s^-1
//
//
//
//double Param::FREE_CHER_bind_ACTIVE = 1.0e6/(Param::NA*Param::CELL_VOLUME);//2.85e6/(Param::NA*Param::CELL_VOLUME);//0.0001e6/(NA*CELL_VOLUME);
//double Param::FREE_CHEB_bind_ACTIVE = 1.0e6/(Param::NA*Param::CELL_VOLUME);//1.98e6/(Param::NA*Param::CELL_VOLUME);//0.000006e6/(NA*CELL_VOLUME);
//double Param::CHER_unbind_ACTIVE = 0.4;
//double Param::CHEB_unbind_ACTIVE = 0.4;
//double Param::ACTIVE_BOUND_CHER_bind_TETHER = 3; //3;//4; //20; // s^-1
//double Param::ACTIVE_BOUND_CHEB_bind_TETHER = 3;//4; //20; // s^-1
//
//
//double Param::CHER_meth_RECEPTOR = 4.5;  //0.4 // per second
//double Param::CHEB_demeth_RECEPTOR = 7.2; //0.6; // per second
//
//double Param::AUTO_phos_CHEA = 11;    // per second
//double Param::CHEA_phos_CHEY = 90e6/(Param::NA*Param::CELL_VOLUME);  // 100e6/(Na*V)
//double Param::CHEA_phos_CHEB = 10e6/(Param::NA*Param::CELL_VOLUME);   // 10e6/(Na*V)
//double Param::AUTO_dephos_CHEY = 35;  // per second
//double Param::AUTO_dephos_CHEB = 1;   // per second

double Param::FREE_CHER_bind_TETHER = 4e6/(Param::NA*Param::CELL_VOLUME);//0.01e6/(NA*CELL_VOLUME);
double Param::FREE_CHEB_bind_TETHER = 3e6/(Param::NA*Param::CELL_VOLUME); //
double Param::TETHERED_CHER_bind_ACTIVE = 14;//10; // s^-1
double Param::TETHERED_CHEB_bind_ACTIVE = 14;//10; // s^-1
double Param::CHER_unbind_TETHER = 0.5; //0.2;//0.2; // s^-1
double Param::CHEB_unbind_TETHER = 0.5; //0.2; // s^-1

double Param::FREE_CHER_bind_ACTIVE = 0.002e6/(Param::NA*Param::CELL_VOLUME);//2.85e6/(Param::NA*Param::CELL_VOLUME);//0.0001e6/(NA*CELL_VOLUME);
double Param::FREE_CHEB_bind_ACTIVE = 0.002e6/(Param::NA*Param::CELL_VOLUME);//1.98e6/(Param::NA*Param::CELL_VOLUME);//0.000006e6/(NA*CELL_VOLUME);
double Param::CHER_unbind_ACTIVE = 50;
double Param::CHEB_unbind_ACTIVE = 50;
double Param::ACTIVE_BOUND_CHER_bind_TETHER = 2; //3;//4; //20; // s^-1
double Param::ACTIVE_BOUND_CHEB_bind_TETHER = 2;//4; //20; // s^-1

double Param::CHER_meth_RECEPTOR = 6;  //0.4 // per second
double Param::CHEB_demeth_RECEPTOR = 7; //0.6; // per second

double Param::AUTO_phos_CHEA = 15;    // per second
double Param::CHEA_phos_CHEY = 90e6/(Param::NA*Param::CELL_VOLUME);  // 100e6/(Na*V)
double Param::CHEA_phos_CHEB = 10e6/(Param::NA*Param::CELL_VOLUME);   // 10e6/(Na*V)
double Param::AUTO_dephos_CHEY = 24;  // per second
double Param::AUTO_dephos_CHEB = 1;   // per second


////////////////////////////////////////////////////////////////
//Receptor ligand binding constants
//double Param::asp_Koff_TAR = 0.0094e-3; //0.02e-3; // This is 0.02 mM or 0.00002 M
//double Param::asp_Kon_TAR = 0.017e-3; //0.5e-3; // mM
//double Param::asp_Koff_TSR = 150e-3; //100e-3; // mM
//double Param::asp_Kon_TSR = 245e-3; //1000; //Yes, 1000 M, or 10^6 mM

double Param::asp_Koff_TAR = 0.02e-3; // This is 0.02 mM or 0.00002 M
double Param::asp_Kon_TAR = 0.5e-3; // mM
double Param::asp_Koff_TSR = 100e-3; // mM
double Param::asp_Kon_TSR = 1000; //Yes, 1000 M, or 10^6 mM



//Motor Parameters
double Param::Mot_Kd = 3.06;
double Param::Mot_g0 = 6.7;

//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
double Param::Mot_g1 = 35;
double Param::Mot_omega = 100;

//Calculated values from parameters
int Param::sizeOfCluster = Param::tarDimerPerClusterCount+Param::tsrDimerPerClusterCount;
int Param::tarDimerCount = Param::tarDimerPerClusterCount*Param::receptorClusterCount;
int Param::tsrDimerCount = Param::tsrDimerPerClusterCount*Param::receptorClusterCount;
int Param::cheAcount = Param::cheAPerClusterCount*Param::receptorClusterCount;





/// End parameters
//////////////////////////////



void Param::refactor()
{
	Param::sizeOfCluster = Param::tarDimerPerClusterCount+Param::tsrDimerPerClusterCount;
	Param::tarDimerCount = Param::tarDimerPerClusterCount*Param::receptorClusterCount;
	Param::tsrDimerCount = Param::tsrDimerPerClusterCount*Param::receptorClusterCount;
	Param::cheAcount = Param::cheAPerClusterCount*Param::receptorClusterCount;
}

