#ifndef NG_PARAM_HH_
#define NG_PARAM_HH_


#include "../../../NFutil/NFutil.hh"



/**
 * 
 */
class NGparam {
	
	
	public:
	
		NGparam();
		~NGparam();
	
		void setLiteSystem();
		void setFullSystem();
		
		
		void setDefaultOutputFileNames();
		void setOutputFileNames(const char * moleculeOutput, const char * key, const char * activity, const char * motor);
		
		void setFullLite();
		void setRichMediaConcentrations();
		void setMinimalMediaConcentrations();
		
		
		void setInitReceptorMethToZero();
		void setInitReceptorMethToEight();
		
		void setToNoiseOutput();
		void setToActivityOutput();
		void setWillParameters();
	
		
		
		
		void setCheR_2x();
		void setCheR_4x();
		void setCheR_8x();
		
		//Global Parameters
		bool get_useTether() const { return useTether; };
		bool get_useRetether() const { return useRetether; };
		bool get_useNeighborRxns() const { return useNeighborRxns; };
		bool get_useCheBFeedback() const { return useCheBFeedback; };
		bool get_useCheA() const { return useCheA; };
		bool get_useCheY() const { return useCheY; };
		bool get_useCheB() const { return useCheB; };
		bool get_useCheR() const { return useCheR; };
		bool get_useMotor() const { return useMotor; };
		
		//Molecule Counts
		int get_tarCount() const { return tarCount; };
		int get_tsrCount() const { return tsrCount; };
		int get_cheAcount() const { return cheAcount; };
		int get_cheRcount() const { return cheRcount; };
		int get_cheBcount() const { return cheBcount; };
		int get_cheYcount() const { return cheYcount; };
		int get_motorCount() const { return motorCount; };
		
		//Starting conditions of the molecules / environment
		double get_aspartateConcentration() const { return aspartateConcentration; };
		int get_cheA_phosState() const { return cheA_phosState; };
		int get_cheB_phosState() const { return cheB_phosState; };
		int get_cheY_phosState() const { return cheY_phosState; };
		int get_motorStartState() const { return motorStartState; };
		int get_tarMethLevelCount(int level) const { return tarMethLevel[level]; }; 
		int get_tsrMethLevelCount(int level) const { return tsrMethLevel[level]; };
		
		//Names of binding sites and states of molecules
		const char * get_nameReceptorTetherSite() const { return receptorTetherSiteName; };
		const char * get_nameReceptorActiveSite() const { return receptorActiveSiteName; };
		const char * get_nameReceptorMethState() const { return receptorMethStateName; };
		const char * get_nameReceptorTypeState() const { return receptorTypeStateName; };
		
		const char * get_nameCheRtetherSite() const { return cheRtetherSiteName; };
		const char * get_nameCheBtetherSite() const { return cheBtetherSiteName; };
		const char * get_nameCheRactiveSite() const { return cheRactiveSiteName; };
		const char * get_nameCheBactiveSite() const { return cheBactiveSiteName; };
		
		const char * get_nameCheAphosState() const { return cheAphosSiteName; };
		const char * get_nameCheBphosState() const { return cheBphosSiteName; };
		const char * get_nameCheYphosState() const { return cheYphosSiteName; };
		const char * get_nameMotorRotationState() const { return motorRotationStateName; };
		

		
		
		
		//Cluster parameters
		
		unsigned int get_clusterTarCountPerCluster() const { return tarPerCluster; };
		unsigned int get_clusterTsrCountPerCluster() const { return tsrPerCluster; };
		
		double get_clusterEnergy_e0() const { return e0; };
		double get_clusterEnergy_e1() const { return e1; };
		 
		double get_clusterTarFreeEnergyOffset(int level) const { return TAR_freeEnergyOffset[level]; }; 
		double get_clusterTsrFreeEnergyOffset(int level) const { return TSR_freeEnergyOffset[level]; };
		
		double get_clusterTarAspKoff() const { return asp_Koff_TAR; };
		double get_clusterTarAspKon() const { return asp_Kon_TAR; };
		double get_clusterTsrAspKoff() const { return asp_Koff_TSR; };
		double get_clusterTsrAspKon() const { return asp_Kon_TSR; };
		
		
		
		
		
		
		
		
		
		
		//Filenames
		const char * get_fileNameMoleculeOutput() const { return moleculeOutputFileName; };
		const char * get_fileNameGroupKey() const { return keyFileName; };
		const char * get_fileNameGroupActivity() const { return activityFileName; };
		const char * get_fileNameMotorOutput() const { return motorOutputFileName; };
		
		
		//Get observables to output
		bool get_outputActivity() const { return outputActivity; };
		bool get_outputTarMethState() const { return outputTarMethState; };
		bool get_outputTsrMethState() const { return outputTsrMethState; };
		bool get_outputAllMethState() const { return outputAllMethState; };
		bool get_outputCheAphos() const { return outputCheAphos; };
		bool get_outputCheBphos() const { return outputCheBphos; };
		bool get_outputCheYphos() const { return outputCheYphos; };
		bool get_outputMotorState() const { return outputMotorState; };
		
		
		
		
		
		
		//Rate constants
		double get_cellVolume() const { return cellVolume; };
		double get_rateFREE_CHER_bind_TETHER() const { return FREE_CHER_bind_TETHER; };
		double get_rateFREE_CHEB_bind_TETHER() const { return FREE_CHEB_bind_TETHER; };
		double get_rateCHER_unbind_TETHER() const { return CHER_unbind_TETHER; };
		double get_rateCHEB_unbind_TETHER() const { return CHEB_unbind_TETHER; };
		
		double get_rateTETHERED_CHER_bind_ACTIVE() const { return TETHERED_CHER_bind_ACTIVE; };
		double get_rateTETHERED_CHER_bind_NEIGHBOR_ACTIVE() const { return TETHERED_CHER_bind_NEIGHBOR_ACTIVE; };
		double get_rateTETHERED_CHEB_bind_ACTIVE() const { return TETHERED_CHEB_bind_ACTIVE; };
		double get_rateTETHERED_CHEB_bind_NEIGHBOR_ACTIVE() const { return TETHERED_CHEB_bind_NEIGHBOR_ACTIVE; };
		
		double get_rateFREE_CHER_bind_ACTIVE() const { return FREE_CHER_bind_ACTIVE; };
		double get_rateFREE_CHEB_bind_ACTIVE() const { return FREE_CHEB_bind_ACTIVE; };
		double get_rateCHER_unbind_ACTIVE() const { return CHER_unbind_ACTIVE; };
		double get_rateCHEB_unbind_ACTIVE() const { return CHEB_unbind_ACTIVE; };
		
		double get_rateACTIVE_BOUND_CHER_bind_TETHER() const { return ACTIVE_BOUND_CHER_bind_TETHER; };
		double get_rateACTIVE_BOUND_CHEB_bind_TETHER() const { return ACTIVE_BOUND_CHEB_bind_TETHER; };
		
		double get_rateCHER_meth_RECEPTOR() const { return CHER_meth_RECEPTOR; };
		double get_rateCHEB_demeth_RECEPTOR() const { return CHEB_demeth_RECEPTOR; };
		
		double get_rateAUTO_phos_CHEA() const { return AUTO_phos_CHEA; };
		double get_rateCHEA_phos_CHEY() const { return CHEA_phos_CHEY; };
		double get_rateCHEA_phos_CHEB() const { return CHEA_phos_CHEB; };
		double get_rateAUTO_dephos_CHEY() const { return AUTO_dephos_CHEY; };
		double get_rateAUTO_dephos_CHEB() const { return AUTO_dephos_CHEB; };			
		
		
		unsigned int get_receptorDimerNumberOfMethSites() const { return numberOfDimerMethylationSites; };
		
		
		double get_motKd() const { return Mot_Kd; };
		double get_motg0() const { return Mot_g0; };
		double get_motg1() const { return Mot_g1; };
		double get_motOmega() const { return Mot_omega; };
		
	
	protected:
		
		void initNames();
		void initEverything();
		
		//Output FileNames
		const char * moleculeOutputFileName;
		const char * keyFileName;
		const char * activityFileName;
		const char * motorOutputFileName;
		
		//Outputters
		bool outputActivity;
		bool outputTarMethState;
		bool outputTsrMethState;
		bool outputAllMethState;
		bool outputCheAphos;
		bool outputCheBphos;
		bool outputCheYphos;
		bool outputMotorState;
		
		
		
		
		//Binding site and State names (so that mistakes are minimized)
		const char * receptorTetherSiteName;
		const char * receptorActiveSiteName;
		const char * receptorMethStateName;
		const char * receptorTypeStateName;
		
		const char * cheRtetherSiteName;
		const char * cheBtetherSiteName;
		const char * cheRactiveSiteName;
		const char * cheBactiveSiteName;
		
		const char * cheAphosSiteName;
		const char * cheBphosSiteName;
		const char * cheYphosSiteName;
		const char * motorRotationStateName;
		
		
		
		double cellVolume;
		
		
		//Global Model Parameters
		bool useTether;
		bool useRetether;
		bool useNeighborRxns;
		bool useCheBFeedback;
		bool useCheA;
		bool useCheY;
		bool useCheB;
		bool useCheR;
		bool useMotor;
		
		bool useWillCluster;
		
		/* This parameter controls whether or not binding to the active
		 * site of a receptor depends on the methylation / demethylation level
		 * so that CheR binds faster to the active site if there are more
		 * sites where CheR can methylate */
		bool useRB_aSiteBindingDependsOnAvail;
		
		
		//Molecule Counts
		int tarCount;
		int tsrCount;
		int cheAcount;
		int cheRcount;
		int cheBcount;
		int cheYcount;
		int motorCount;
		
		

		
		
		//Starting conditions of the molecules / environment
		double aspartateConcentration;
		int cheA_phosState;
		int cheB_phosState;
		int cheY_phosState;
		int motorStartState;
		int tarMethLevel [9]; 
		int tsrMethLevel [9];
		
		
		

		//Receptor Cluster Parameters
		unsigned int tarPerCluster;
		unsigned int tsrPerCluster;
		
		//For will's parameters;
		double e0;
		double e1;
		
		double TAR_freeEnergyOffset [9];
		double TSR_freeEnergyOffset [9];
		double asp_Koff_TAR;
		double asp_Kon_TAR;
		double asp_Koff_TSR;
		double asp_Kon_TSR;
		unsigned int numberOfDimerMethylationSites;
		
		
		
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
		double CHER_meth_RECEPTOR;   // per second
		double CHEB_demeth_RECEPTOR; // per second
		
		//Phosphorylation Rates
		double AUTO_phos_CHEA;
		double CHEA_phos_CHEY;
		double CHEA_phos_CHEB;
		double AUTO_dephos_CHEY;
		double AUTO_dephos_CHEB;
		
		
		
		//Motor Parameters
		double Mot_Kd;
		double Mot_g0;
		
		//Warning!  these change nothing!  (now, wFactor = 1.02, see AN_system_full.cpp)
		double Mot_g1;
		double Mot_omega;
		
		
	
	
	
		
		
		
		//Tethering rate
		// Binding of pentapeptide NWETF (C-terminus of major receptors) (Yi & Weis, 2002)
		//  Barnakov, A. N., Barnakova, L. A., & Hazelbauer, G. L. (2002) J. Biol. Chem. 277, 42151-42156
		//  Kd = 11uM for CheR
		//  kd = 150uM for CheB
		//  I assume that the unbinding rate is 0.1/s, as in Hansen, Endres, Wingreen, Plos Comp Bio, 2007,
		//  therefore, rate of cheR binding is ~0.01 / uM s
		//  and cheB binding is ~0.0006/uM s
	


};
		
		

#endif /*NG_PARAM_HH_*/

		
		
		
		
		