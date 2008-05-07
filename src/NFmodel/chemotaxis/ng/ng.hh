#ifndef NG_HH_
#define NG_HH_


#include "../chemotaxis.hh"
#include "../../../NFreactions/rate_change/NFrateChangeReactions.hh"

//Define some global constants for use in remembering index
//values in the system






using namespace NFcore;

namespace NG {
	
	
	void run(int argc, char *argv[]);
	
	
	// A class that nicely wraps all the parameters of the system
	class Param {
		
		public:
			//Universal Values
			static const double CELL_VOLUME = 1.41e-15;
			static const double NA = 6.0221415e23;
		
			//Main simulation parameters
			static double aspConc;
			static bool outputGroupStats;
			static bool outputMotorStats;
			static bool useFullSystem;
			static bool useTetherRxns;
			static bool useCheBFeedbackLoop;
			
			
			//Select the observables to use
			static bool useCheAPhosObs;
			static bool useCheBPhosObs;
			static bool useMotorObs;
			static bool useMethLevelObs;
			static bool useRB_boundTetherObs;
			static bool useRB_boundActiveSiteObs;
			
			
	
			////////////////////////////////////////////////////////////////
			//Initial counts of molecules and size of clusters
			static int receptorClusterCount;
			static int tarDimerPerClusterCount;
			static int tsrDimerPerClusterCount;
			static int cheAPerClusterCount;
			
			static int cheRcount;
			static int cheBcount;
			static int cheYcount;
			static int motorCount; 
			
			
			
			
			
			//Starting conditions of the molecules
			static int cheB_phosState;
			static int cheA_phosState;
			static int cheY_phosState;
			static int motorStartState;
			static int tarMethLevel; 
			static int tsrMethLevel;
			
			
			////////////////////////////////////////////////////////////////
			//Rate constants
			
			
			
			//Tethering rate
			// Binding of pentapeptide NWETF (C-terminus of major receptors) (Yi & Weis, 2002)
			//  Barnakov, A. N., Barnakova, L. A., & Hazelbauer, G. L. (2002) J. Biol. Chem. 277, 42151-42156
			//  Kd = 11uM for CheR
			//  kd = 150uM for CheB
			//  I assume that the unbinding rate is 0.1/s, as in Hansen, Endres, Wingreen, Plos Comp Bio, 2007,
			//  therefore, rate of cheR binding is ~0.01 / uM s
			//  and cheB binding is ~0.0006/uM s
		
	///	double RB_BindingFactor = 5;
		
	//	double R_BIND_TETHER_RATE = RB_BindingFactor* //0.01e6/(Na*CellVolume);
	//	double B_BIND_TETHER_RATE = RB_BindingFactor*0.0006e6/(Na*CellVolume); //0.0006e6/(Na*CellVolume);
	//	double R_UNBIND_TETHER_RATE = RB_BindingFactor*0.1;
	//	double B_UNBIND_TETHER_RATE = RB_BindingFactor*0.1;
	
	//	double autoPhosARate = 23.5;  //kpauto
	//	double BpBindReceptorRate = 1.982140e6/(Na*CellVolume);  //kpTB
	//	double RBindReceptorRate = 2.857143e6/(Na*CellVolume);  //kpTR
	//	double BUnbindReceptorRate = 1.25; //kmTB
		//double RUnbindReceptorRate = 1.25; //kmTR
		
	//	double BpDemethReceptorRate = 0.6; //1.4;//0.6; //kcTB
	//	double RmethReceptorRate = 0.75; //0.3; //0.75; //kcTR
	
			
			//CheR & CheB tethering rxns
			static double FREE_CHER_bind_TETHER; // = 0.01e6/(NA*CELL_VOLUME);
			static double FREE_CHEB_bind_TETHER; //
			static double TETHERED_CHER_bind_ACTIVE; // s^-1
			static double TETHERED_CHEB_bind_ACTIVE; // s^-1
			static double CHER_unbind_TETHER; // s^-1
			static double CHEB_unbind_TETHER; // s^-1
			
			
			
			static double FREE_CHER_bind_ACTIVE;//0.0001e6/(NA*CELL_VOLUME);
			static double FREE_CHEB_bind_ACTIVE;//0.000006e6/(NA*CELL_VOLUME);
			static double CHER_unbind_ACTIVE;
			static double CHEB_unbind_ACTIVE;
			static double ACTIVE_BOUND_CHER_bind_TETHER; // s^-1
			static double ACTIVE_BOUND_CHEB_bind_TETHER; // s^-1
			
			
			static double CHER_meth_RECEPTOR;   // per second
			static double CHEB_demeth_RECEPTOR; // per second
			
			static double AUTO_phos_CHEA;    // per second
			static double CHEA_phos_CHEY;  // 100e6/(Na*V)
			static double CHEA_phos_CHEB;   // 10e6/(Na*V)
			static double AUTO_dephos_CHEY;  // per second
			static double AUTO_dephos_CHEB ;   // per second
			
			
			
			////////////////////////////////////////////////////////////////
			//Receptor ligand binding constants
			static double asp_Koff_TAR; // This is 0.02 mM or 0.00002 M
			static double asp_Kon_TAR; // mM
			static double asp_Koff_TSR; // mM
			static double asp_Kon_TSR; //Yes, 1000 M, or 10^6 mM
			
			//Motor Parameters
			static double Mot_Kd;
			static double Mot_g0;
			
			//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
			static double Mot_g1;
			static double Mot_omega;
		
			//Calculated values from parameters
			static int sizeOfCluster;
			static int tarDimerCount;
			static int tsrDimerCount;
			static int cheAcount;
			
			
			static void refactor();
	};
	
	
	
	
	System * init_NG_system(
		const char * outputFilename,
		const char * keyFileName,
		const char * pOnFileName,
		const char * motorFileName);
	
	
	
	
	MoleculeType * makeRecDimer(System * s);
	MoleculeType * makeCheA(System * s);
	MoleculeType * makeCheR(System *s);
	MoleculeType * makeCheB(System *s);
	MoleculeType * makeCheY(System *s);
	
	void makeTarMolecules(MoleculeType * receptorDimer);
	void makeTsrMolecules(MoleculeType * receptorDimer);
	
	
	ReactionClass * createMotor_cw2ccw(MoleculeType * Motor, Observable * cheYp, ofstream *o);
	ReactionClass * createMotor_ccw2cw(MoleculeType * Motor, Observable * cheYp, ofstream *o);
	
	
	void rxn_freeRB_bind_unbind_tether(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB);
	void rxn_freeRB_bind_unbind_active(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB);
	void rxn_tetheredRB_bind_active(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB);
	
	void rxn_RB_meth_demeth(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB);
	
	void rxn_retether(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB);
	
	
	void rxn_CheA_auto_phos(System * s, MoleculeType * cheA);
	
	
	void rxn_CheA_phos_CheY(System * s, MoleculeType * cheA, MoleculeType *cheY);
	void ANrxn_add_CheY_autodephosFull(System *s, MoleculeType * cheY);
	void rxn_CheA_phos_CheBFull(System * s, MoleculeType * cheA, MoleculeType *cheB);
	void rxn_add_CheB_autodephosFull(System *s, MoleculeType * cheB);
	
	
	
	
	void addReceptorMethLevelObs(MoleculeType *recDimer);
	void addCheRB_boundToTether(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer);
	void addCheRB_boundToTetherAndActiveSite(MoleculeType * cheR, MoleculeType * cheB, MoleculeType * recDimer);
	void addCheAphos(MoleculeType * cheA);
	void addCheBphos(MoleculeType * cheB);
	void addCheYphos(MoleculeType * cheY);
	void addMotorObs(MoleculeType * motor);
	
	
	void rxn_tetheredRB_bind_activeNeighbor(System *s, MoleculeType *recDimer, MoleculeType *che, bool isCheR,char *neighborSiteName, char *reciprocalNeighborSiteName);
	
	void run_NG_stepper(System *s);
	void run_NG_stepAndRemove(System *s);
	void run_NG_noiseTest(System *s, int time, int samples);
	
	
	
	
	///////////////////////////////////////////////////////////////////////
	// DimerGroup definition - to keep track of available active sites per dimer
	class DimerGroup:public Group
	{
		public:
		
			DimerGroup(char *name, System *s, int methStateIndex, unsigned int numOfMethSites);
			virtual ~DimerGroup();
			
			virtual void addToGroup(Molecule * m);
			virtual void notify(Molecule *changedMolecule, int oldStateValue, int newStateValue);
			
			/* function that lets us update ligand concentration for this group */
			virtual void updateGroupProperty(double * values, int n_values);
			virtual void printDetails();
			
			static const int METH_SITE_LEVEL = 0;
			static const int FREE_SITE_LEVEL = 1;
			
		protected:
			bool isFull;
			unsigned int numOfMethSites;
	};
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Allows binding of cheR or cheB to the active site of methylation
	class NG_FreeActiveSiteBinding : public DORrxn 
	{
		public:
			NG_FreeActiveSiteBinding( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int DORreactantIndex,
								char * DORgroupName,
								int DORgroupValueIndex,
								double baseRate,
								const char *cheBindingSite,
								const char *recActiveSite);
			virtual ~NG_FreeActiveSiteBinding() {};
			virtual void transformReactants(Molecule ** reactants, int nReactants);
			
		protected:
			int cheBindingSiteIndex;
			int recActiveSiteIndex;
	};
	
	
	// Allows a tethered cheR or cheB to bind the active site of methylation
	class NG_TetheredActiveSiteBinding : public DORrxn 
	{
		public:
			NG_TetheredActiveSiteBinding( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int DORreactantIndex,
								char * DORgroupName,
								int DORgroupValueIndex,
								double baseRate,
								int cheBindingSiteIndex,
								const char *recActiveSite,
								const char *tetherSite);
			virtual ~NG_TetheredActiveSiteBinding() {};
			virtual void transformReactants(Molecule ** reactants, int nReactants);
			
		protected:
			int cheBindingSiteIndex;
			int recActiveSiteIndex;
			int tetherSiteIndex;
	};
	
	// Allows a tethered cheR or cheB to bind the active site of a neighbor
	class NG_TetheredActiveNeighborSiteBinding : public DORrxn 
	{
		public:
			NG_TetheredActiveNeighborSiteBinding( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int DORreactantIndex,
								char * DORgroupName,
								int DORgroupValueIndex,
								double baseRate,
								const char *neighborSite,
								const char *recTetherSite,
								const char *recActiveSite,
								int cheBindingSiteIndex);
			virtual ~NG_TetheredActiveNeighborSiteBinding() {};
			virtual void transformReactants(Molecule ** reactants, int nReactants);
			
		protected:
			int neighborSiteIndex;
			int recTetherSiteIndex;
			int cheActiveSiteIndex;
			int recActiveSiteIndex;
	};
	
	
	
	// Allows for the methylation or demethylation or receptors
	class NG_BR_MethDemeth : public DORrxn
	{
		public:
			NG_BR_MethDemeth( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int DORreactantIndex,
								char * DORgroupName,
								int DORgroupValueIndex,
								double baseRate,
								const char *receptorMethStateName,
								const char *receptorActiveSiteName,
								int deltaMethState);
			virtual ~NG_BR_MethDemeth() { };
			virtual void transformReactants(Molecule ** reactants, int nReactants);
		
		protected:
		
			int receptorMethStateIndex;  //this is the state that maintains receptor meth level
			int deltaMethState;          //this will be +1 or -1 depending if we are meth or demeth
			int receptorActiveSiteName;  //this is the index values of the bonds that connect receptors
			                             //in the cluster (will be of length 6 for a hex grid)
	};
	
	
	// Allows for the methylation or demethylation or receptors
	class NG_BR_MethDemeth2 : public ReactionClass
	{
		public:
			NG_BR_MethDemeth2( char * name, 
											int n_reactants, 
											TemplateMolecule ** reactantTemplates,
											double rate,
											const char *receptorMethStateName,
											const char *receptorActiveSiteName,
											int deltaMethState);
			virtual ~NG_BR_MethDemeth2() { };
			virtual void transformReactants(Molecule ** reactants, int nReactants);
		
		protected:
		
			int receptorMethStateIndex;  //this is the state that maintains receptor meth level
			int deltaMethState;          //this will be +1 or -1 depending if we are meth or demeth
			int receptorActiveSiteName;  //this is the index values of the bonds that connect receptors
			                             //in the cluster (will be of length 6 for a hex grid)
	};
	
	

	
	//Allows a moleucle to retether to a receptor if it is bound to the active site
	class NG_retether : public ReactionClass {
	
		public:
			NG_retether(	char * name, 
									int n_reactants, 
									TemplateMolecule ** reactantTemplates, 
									double rate,
									int recTetherSiteIndex,
									int cheActiveSiteIndex,
									int cheTetherSiteIndex);
			virtual ~NG_retether() {};
			virtual void transformReactants(Molecule ** reactants, int nReactants);
		private:
			int recTetherSiteIndex;
			int cheActiveSiteIndex;
			int cheTetherSiteIndex;
	};
	
	
	
	class NG_AutoPhosA : public DORrxn 
	{
		public:
			NG_AutoPhosA( char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates,
								int DORreactantIndex,
								char * DORgroupName,
								int DORgroupValueIndex,
								double baseRate,
								const char * stateName,
								int newStateValue);
			virtual ~NG_AutoPhosA() {};
			virtual void transformReactants(Molecule ** reactants, int nReactants);
			
		protected:
			int stateIndex;
			int newStateValue;
	};
}







#endif /*NG_HH_*/
