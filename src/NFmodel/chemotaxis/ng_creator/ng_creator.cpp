

#include "ng_creator.hh"



using namespace NFcore;
using namespace NG;
using namespace ReceptorCluster;




System * create(
		const char *systemName, 
		NGparam &p)
{

	/////////////////////////////////////////////////////////////////////////////////////
	//Make the actual system
	System *s = new System(systemName);
	s->registerOutputFileLocation(p.get_fileNameMoleculeOutput());
	
	
	/////////////////////////////////////////////////////////////////////////////////////
	//Create all the molecule types and initial molecule instances
	MoleculeType *receptorDimer = createReceptorDimer(s,p);
	MoleculeType *cheA = 0;
	MoleculeType *cheR = 0;
	MoleculeType *cheB = 0;
	MoleculeType *cheY = 0;
	MoleculeType *motor = 0;
	
	if(p.get_useCheA())
		cheA = createCheA(s,p);
	if(p.get_useCheR())
		cheR = createCheR(s,p);
	if(p.get_useCheB())
		cheB = createCheB(s,p);
	if(p.get_useCheY())
		cheY = createCheY(s,p);
	if(p.get_useMotor()) {
		motor = createMotor(s,p);
	}
		
	
	
	
	//Create and populate the clusters of receptors
	vector <Molecule *> recDimerMolecules;
	for(int r=0; r<receptorDimer->getMoleculeCount(); r++)
		recDimerMolecules.push_back(receptorDimer->getMolecule(r));
	int clusterCount = createFixedSizeClustersWithConstantRatio(s,recDimerMolecules,p);
	
	
	
	//Add each dimer to its own group so it can bind cheR with a rate proportional to its available active sites
	for(int r=0; r<receptorDimer->getMoleculeCount(); r++) {
		Group *g = new DimerGroup(DIMER_GROUP_NAME,s,receptorDimer->getStateIndex(p.get_nameReceptorMethState()), p.get_receptorDimerNumberOfMethSites());
		g->addToGroup(receptorDimer->getMolecule(r));	
	}
	
	
	
	
	
	
	if(p.get_useCheA()) {
		//add cheA to clusters and create the autoPhos cheA reaction
		vector <Molecule *> cheAMolecules;
		for(int a=0; a<cheA->getMoleculeCount(); a++)
			cheAMolecules.push_back(cheA->getMolecule(a));
		distributeRandomCheAtoClusters_Evenly(s,cheAMolecules, clusterCount);
		
		createRxn_CheA_auto_phos(s,cheA,p);
	}
	
	if(p.get_useCheY() && p.get_useCheA()) 
	{
		//add CheA phos CheY
		createRxn_CheA_phos_CheY(s,cheA,cheY,p);
		
		//add cheY autodephos
		createRxn_add_CheY_autodephos(s,cheY,p);
	}
	
	if(p.get_useCheR()) {
		if(p.get_useTether()) {
			//add bind / unbind tether
			//add tethered binds active site
		}
		
		if(p.get_useNeighborRxns()) {
			//add bind neighbor active site
		}
		
		if(p.get_useRetether()) {
			//add active site rebinds tether
		}
		
		//Active site binding and methylation reactions
		createRxn_freeR_bind_unbind_active(s,receptorDimer, cheR, p);
		createRxn_R_meth(s, receptorDimer, cheR, p);
	}
	
	if(p.get_useCheB()) {
		
		if(p.get_useCheBFeedback() && p.get_useCheA())
		{
			//add CheA phos CheB
			//add CheB autodephos
			createRxn_CheA_phos_CheB(s,cheA,cheB,p);
			createRxn_add_CheB_autodephos(s,cheB,p);
			
		}
		
		if(p.get_useTether())
		{
			//add bind unbind tether
		}
		
		
		if(p.get_useNeighborRxns()) {
			//add bind neighbor active site
		}
				
		if(p.get_useRetether()) {
			//add active site rebinds tether
		}
				
		//Active site binding and demethylation reactions
		createRxn_freeB_bind_unbind_active(s,receptorDimer, cheB, p);
		createRxn_B_demeth(s, receptorDimer, cheB, p);
	}
	
	

	
	
	if(p.get_useCheY() && p.get_useMotor())
	{
		//add motor reactions
		Motor::createRxn_MotorSwitching(s,motor,cheY,p);
	}
	
	
	
	
	
	//Be sure to output everything we need...
	

		//if(Param::useRB_boundTetherObs) NG::addCheRB_boundToTether(cheR, cheB, recDimer);
		//if(Param::useRB_boundActiveSiteObs) NG::addCheRB_boundToTetherAndActiveSite(cheR, cheB, recDimer);
	
	
	
	if(p.get_outputTarMethState()) {}
	if(p.get_outputTsrMethState()) {}
	if(p.get_outputAllMethState()) { createObservable_ReceptorMethLevel(receptorDimer, p); }
	
	if(p.get_outputCheAphos() && p.get_useCheA()) { addCheAphos(cheA); } 
	
	if(p.get_outputCheBphos()) {addCheBphos(cheB);}
	if(p.get_outputCheYphos() && !p.get_useMotor()) { addCheYphos(cheY); }
	if(p.get_outputMotorState()) { createObservable_MotorState(motor, p); }
	
//	addCheRB_boundToTetherAndActiveSite(cheR, cheB, receptorDimer);
	
	
	
	//This code sets up the outputter for the activity of the clusters
	if(p.get_outputActivity())
	{
		vector <TemplateMolecule *> tm;
		vector <const char *> templateNames;
		vector <const char *> filenames;
		vector <unsigned int> values;
		
		TemplateMolecule *rec = new TemplateMolecule(receptorDimer);
		rec->addStateValue(p.get_nameReceptorTypeState(),TAR);
		tm.push_back(rec);
		templateNames.push_back("TarDimer");
		
		rec = new TemplateMolecule(receptorDimer);
		rec->addStateValue(p.get_nameReceptorTypeState(),TSR);
		tm.push_back(rec);
		templateNames.push_back("TsrDimer");
		
		filenames.push_back(p.get_fileNameGroupActivity());
		values.push_back(P_ON_INDEX);
		
		GroupOutputter * g = new GroupOutputter(s, CLUSTER_NAME, p.get_fileNameGroupKey(), tm, templateNames, filenames, values);
		s->addGroupOutputter(g);
	}
	
	
	
	
	
	//////////////////////////////
	// Check some of the receptor cluster parameters for consistency
//	if((Param::tarDimerCount + Param::tsrDimerCount)%Param::sizeOfCluster!=0)
//	{ cerr<<"Counts of receptors and cluster sizes don't match!"<<endl; exit(1);}
	
	
	
	//Create and populate the clusters
/*	vector <Molecule *> recDimerMolecules;
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
		Group *g = new DimerGroup(DIMER_GROUP_NAME,s,recDimer->getStateIndex("m"));
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
	
	
	
	
*/
	
	
	
	
	
	
	s->prepareForSimulation();
	s->outputAllObservableNames();
	
	s->getMoleculeType(0)->getMolecule(0)->printDetails();
	return s;
}























