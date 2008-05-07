#include "ng.hh"
#include <stdio.h>
#include <stdlib.h>

using namespace NG;
using namespace ReceptorCluster;
using namespace Motor;





void NG::run(int argc, char *argv[])
{
	cout<<"Running The NG Chemotaxis Simulation..."<<endl;
	
	char * outputFilename = "/home/msneddon/Desktop/NG_sim/defaultOut/NG_basicOutput.txt";
	char * keyFileName = "/home/msneddon/Desktop/NG_sim/defaultOut/NG_groupKey.txt";
	char * pOnFileName = "/home/msneddon/Desktop/NG_sim/defaultOut/NG_pOnVals.txt";
	char * motorFileName = "/home/msneddon/Desktop/NG_sim/defaultOut/NG_motorOutput.txt";
	
	if(argc>4)
	{
		outputFilename = argv[4];
		keyFileName = argv[5];
		pOnFileName = argv[6];
		motorFileName = argv[7];
	}
	
	if(argc>3)
	{
		if(strncmp(argv[3],"noise",5)==0)
		{
			Param::aspConc = 0;
			Param::outputGroupStats = false;
			Param::outputMotorStats = true;
			Param::useFullSystem = true;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = true;

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = false;
			Param::useMethLevelObs = false;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			
			
//			Param::FREE_CHER_bind_TETHER = 3e6/(Param::NA*Param::CELL_VOLUME);//0.01e6/(NA*CELL_VOLUME);
//			Param::FREE_CHEB_bind_TETHER = 2e6/(Param::NA*Param::CELL_VOLUME); //
//			Param::TETHERED_CHER_bind_ACTIVE = 14;//10; // s^-1
//			Param::TETHERED_CHEB_bind_ACTIVE = 14;//10; // s^-1
//			Param::CHER_unbind_TETHER = 0.01; //0.2;//0.2; // s^-1
//			Param::CHEB_unbind_TETHER = 0.01; //0.2; // s^-1
//
//			Param::FREE_CHER_bind_ACTIVE = 0.002e6/(Param::NA*Param::CELL_VOLUME);//2.85e6/(Param::NA*Param::CELL_VOLUME);//0.0001e6/(NA*CELL_VOLUME);
//			Param::FREE_CHEB_bind_ACTIVE = 0.002e6/(Param::NA*Param::CELL_VOLUME);//1.98e6/(Param::NA*Param::CELL_VOLUME);//0.000006e6/(NA*CELL_VOLUME);
//			Param::CHER_unbind_ACTIVE = 3;
//			Param::CHEB_unbind_ACTIVE = 3;
//			Param::ACTIVE_BOUND_CHER_bind_TETHER = 2; //3;//4; //20; // s^-1
//			Param::ACTIVE_BOUND_CHEB_bind_TETHER = 2;//4; //20; // s^-1
//
//			Param::CHER_meth_RECEPTOR = 7;  //0.4 // per second
//			Param::CHEB_demeth_RECEPTOR = 6; //0.6; // per second
//
//			Param::AUTO_phos_CHEA = 15;    // per second
//			Param::CHEA_phos_CHEY = 90e6/(Param::NA*Param::CELL_VOLUME);  // 100e6/(Na*V)
//			Param::CHEA_phos_CHEB = 10e6/(Param::NA*Param::CELL_VOLUME);   // 10e6/(Na*V)
//			Param::AUTO_dephos_CHEY = 20;  // per second
//			Param::AUTO_dephos_CHEB = 1;   // per second
			
			
			
			
			
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			run_NG_noiseTest(s, 10000, 1000000);
			s->printAllReactions();
			delete s;		
			return;
				
		}
		else if(strncmp(argv[3],"stepper",7)==0)
		{
			Param::aspConc = 0;
			Param::outputGroupStats = true;
			Param::outputMotorStats = true;
			Param::useFullSystem = true;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = false;

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = true;
			Param::useMethLevelObs = true;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			
			
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			run_NG_stepper(s);
			s->printAllReactions();
			delete s;	
			return;
		}
		else if(strncmp(argv[3],"singlestep",10)==0)
		{
			Param::aspConc = 0;
			Param::outputGroupStats = true;
			Param::outputMotorStats = true;
			Param::useFullSystem = true;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = true;

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = true;
			Param::useMethLevelObs = true;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			
			
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			run_NG_stepAndRemove(s);
			s->printAllReactions();
			delete s;	
			
			return;
		}
		else if(strcmp(argv[3],"cheR")==0)
		{
			cout<<"Trying to read what cheR level you want..."<<endl;
			char *cheR = argv[8];
			int cheRlevel = atoi(cheR);
			Param::cheRcount = Param::cheRcount*cheRlevel;
			
			Param::aspConc = 0;
			Param::outputGroupStats = false;
			Param::outputMotorStats = true;
			Param::useFullSystem = true;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = true;

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = false;
			Param::useMethLevelObs = false;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			run_NG_noiseTest(s, 10000, 1000000);
			s->printAllReactions();
			delete s;	
			
			return;
		}
		
		else if(strcmp(argv[3],"cheR_km_test")==0)
		{
			//cout<<"Trying to read what cheR level you want..."<<endl;
			char *charClusterCount = argv[8];
			int clusterCount= atoi(charClusterCount);
			Param::receptorClusterCount = clusterCount;
			//Param::cheRcount = Param::cheRcount*cheRlevel;
			
			//cout<<"here"<<endl;
			Param::aspConc = 0;
			Param::outputGroupStats = false;
			Param::outputMotorStats = false;
			Param::useFullSystem = false;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = true;
			

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = false;
			Param::useMethLevelObs = true;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			
			Param::cheBcount = 0;
			Param::cheRcount = 140;
			
			Param::tarMethLevel = 0; 
			Param::tsrMethLevel = 0;
			
			Param::refactor();
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			s->sim(10,100);
			//s->equilibriate()
			//run_NG_noiseTest(s, 1000, 1000);
			s->printAllReactions();
			delete s;	
			
			return;
		}
		else if(strcmp(argv[3],"quickTest")==0)
		{
			
			//cout<<"here"<<endl;
			Param::aspConc = 0;
			Param::outputGroupStats = false;
			Param::outputMotorStats = false;
			Param::useFullSystem = true;
			Param::useTetherRxns = true;
			Param::useCheBFeedbackLoop = true;
			

			Param::useCheAPhosObs = false;
			Param::useCheBPhosObs = false;
			Param::useMotorObs = false;
			Param::useMethLevelObs = true;
			Param::useRB_boundTetherObs = false;
			Param::useRB_boundActiveSiteObs = false;
			

			
//			
			Param::cheRcount = 1000;//140;
			Param::cheBcount = 1000;
			Param::cheYcount = 8000;
			Param::motorCount = 0; 
			Param::receptorClusterCount = 1;
			Param::cheAPerClusterCount = 100;
//			
//			
//			Param::tarMethLevel = 2; 
//			Param::tsrMethLevel = 2;
			
			Param::refactor();
			
			
			
			
			System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
			s->outputAllObservableNames();
			s->sim(5000,5000);
			//s->equilibriate()
			//run_NG_noiseTest(s, 1000, 1000);
			//s->stepTo(100);
			s->printAllReactions();
			delete s;	
			
			return;
		}
		
		else
		{
			cout<<"Could not identify test.  exiting."<<endl;
			return;	
		}
	}
	else
	{
		System *s = init_NG_system(outputFilename, keyFileName, pOnFileName, motorFileName);
		s->outputAllObservableNames();
		s->equilibriate(500,50);
		s->sim(1000,1000);
		s->printAllReactions();
		delete s;	
		return;	
	}
}


void NG::run_NG_stepAndRemove(System *s)
{
	
	
	int n_values = 1;
	double * ligandConcentation = new double[n_values];
	
	
	int time = 1000; int samples = 100;
	
	s->equilibriate(1500,150);
	//s->sim(1000,1000);
	s->sim(300,30);
	ligandConcentation[0] = 0.00003;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	    
	ligandConcentation[0] = 0;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	    
}

void NG::run_NG_noiseTest(System *s, int time, int samples)
{
	int n_values = 1;
	double * ligandConcentation = new double[n_values];
	
	s->equilibriate(300,100);
	ligandConcentation[0] = 0;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
}




void NG::run_NG_stepper(System *s)
{
	int n_values = 1;
	double * ligandConcentation = new double[n_values];
	ligandConcentation[0] = 0.00001;
	
	int time = 600; int samples = 60;
	
	s->equilibriate(1000,100);
	//s->sim(1000,1000);
	s->sim(300,30);
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	    
	ligandConcentation[0] = 0.0001;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples); 
	    
	ligandConcentation[0] = 0.001;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	    
	ligandConcentation[0] = 0.01;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	    
	ligandConcentation[0] = 0.1;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	
	ligandConcentation[0] = 1;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
	
	ligandConcentation[0] = 10;
	s->updateAllGroupProperty(ligandConcentation, n_values);
	s->sim(time,samples);
}




//The Next Generation Chemotaxis System
System * NG::init_NG_system(
	const char * outputFilename,
	const char * keyFileName,
	const char * pOnFileName,
	const char * motorFileName)
{
	
	//////////////////////////////
	// Check some of the receptor cluster parameters for consistency
	if((Param::tarDimerCount + Param::tsrDimerCount)%Param::sizeOfCluster!=0)
	{ cerr<<"Counts of receptors and cluster sizes don't match!"<<endl; exit(1);}
	
	
	/////////////////////////////
	//Make the system
	System *s = new System("The NG Chemotaxis System");
	s->registerOutputFileLocation(outputFilename);
	
	
	/////////////////////////////
	//Create all the molecule types and molecules
	MoleculeType *recDimer = makeRecDimer(s);
	makeTarMolecules(recDimer);
	makeTsrMolecules(recDimer);
	MoleculeType *cheR = makeCheR(s);
	MoleculeType *cheB = makeCheB(s);
	
	
	//Create and populate the clusters
	vector <Molecule *> recDimerMolecules;
	for(int r=0; r<recDimer->getMoleculeCount(); r++)
		recDimerMolecules.push_back(recDimer->getMolecule(r));
	createFixedSizeClustersWithConstantRatio(recDimerMolecules,
		Param::sizeOfCluster, true, s, 
		Param::asp_Koff_TAR, 
		Param::asp_Kon_TAR, 
		Param::asp_Koff_TSR, 
		Param::asp_Kon_TSR, 
		Param::aspConc );
	
	if(Param::useFullSystem)
	{
		MoleculeType *cheA = makeCheA(s);
		MoleculeType *cheY = makeCheY(s);
		MoleculeType *motor = Motor::makeMotor(s);
		
		vector <Molecule *> cheAMolecules;
		for(int a=0; a<cheA->getMoleculeCount(); a++)
			cheAMolecules.push_back(cheA->getMolecule(a));
		distributeRandomCheAtoClusters_Evenly(s,cheAMolecules, Param::receptorClusterCount);
		
		
		
		rxn_CheA_auto_phos(s, cheA);
		rxn_CheA_phos_CheY(s, cheA, cheY);
		ANrxn_add_CheY_autodephosFull(s, cheY);
		
		rxn_CheA_phos_CheBFull(s, cheA, cheB);
		rxn_add_CheB_autodephosFull(s, cheB);
		ANrxn_addMotorSwitching_Full(s, motor, cheY, motorFileName);
		
		
		if(Param::useCheAPhosObs) addCheAphos(cheA);
		if(Param::useCheBPhosObs) addCheBphos(cheB);
		if(Param::useMotorObs) addMotorObs(motor);
		
		
	}
	
	for(int r=0; r<recDimer->getMoleculeCount(); r++)
	{
		Group *g = new DimerGroup(DIMER_GROUP_NAME,s,recDimer->getStateIndex("m"),8);
		g->addToGroup(recDimer->getMolecule(r));	
	}
	
	recDimer->printDetails();
	
	
	if(Param::useTetherRxns)
	{
		rxn_freeRB_bind_unbind_tether(s, recDimer, cheR, cheB);
		rxn_tetheredRB_bind_active(s, recDimer, cheR, cheB);
	 	rxn_retether(s, recDimer, cheR, cheB);
	}

	NG::rxn_freeRB_bind_unbind_active(s, recDimer, cheR, cheB);
	NG::rxn_RB_meth_demeth(s, recDimer, cheR, cheB);
	
	
	
	
	// Observables to output
	if(Param::useMethLevelObs) NG::addReceptorMethLevelObs(recDimer);
	if(Param::useRB_boundTetherObs) NG::addCheRB_boundToTether(cheR, cheB, recDimer);
	if(Param::useRB_boundActiveSiteObs) NG::addCheRB_boundToTetherAndActiveSite(cheR, cheB, recDimer);
	
	
	
	if(Param::outputGroupStats)
	{
		vector <TemplateMolecule *> tm;
		vector <const char *> templateNames;
		vector <const char *> filenames;
		vector <unsigned int> values;
		
		TemplateMolecule *rec = new TemplateMolecule(recDimer);
		rec->addStateValue("type",TAR);
		tm.push_back(rec);
		templateNames.push_back("TarDimer");
		
		rec = new TemplateMolecule(recDimer);
		rec->addStateValue("type",TSR);
		tm.push_back(rec);
		templateNames.push_back("TsrDimer");
		
		filenames.push_back(pOnFileName);
		values.push_back(P_ON_INDEX);
		
		GroupOutputter * g = new GroupOutputter(s, CLUSTER_NAME, keyFileName, tm, templateNames, filenames, values);
		s->addGroupOutputter(g);
	}
	
	//Return the simulation system
	s->prepareForSimulation();
	return s;	
}
