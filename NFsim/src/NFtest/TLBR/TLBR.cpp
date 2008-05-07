#include "TLBR.hh"



#include <math.h>


using namespace NFtest_TLBR;

System * runTLBRSystem(int N_l, int N_r, double cTot, double beta, double koff, char * outputFileName, double simTime, bool outputObs, bool outputAvg, bool outputFinalDist);



void NFtest_TLBR::run(int argc, char *argv[])
{
	cout<<"Running TLBR tests..."<<endl;
	
	if(argc<=4)
	{
		cout<<"Not enough arguements!  "<<endl;
		cout<<"You need to specify first which TLBR test to run and an output filename"<<endl;
		return;
	}
	
	char * filename = argv[4];
	
	
	if( argc>3 )
	{
		if(strncmp(argv[3],"1",1)==0)
		{
//			cout<<" Running test 1: sample trajectory 1 (slope)"<<endl;
//			int N_l = 2000;
//			int N_r = 3000;
//			double ctot = 0.11;
//			double beta = 16.8;
//			double koff = 0.01;
//			double endTime = 10000;
//			cout<<"simming for "<< endTime <<" seconds "<<endl;
//			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false, true, false);
//			delete s;
//			return;
			
			for(int k=0; k<9; k++) {
			cout<<" Running test 1: sample trajectory"<<endl;
			int N_l = 42000;
			int N_r = 3000;
			double koff = 0.01;
			double ctot = 0.378;
			double beta = 0.3;
			double endTime = 3000;
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false, true, false);
			delete s;
			
			}
			return;
			
		}
		else if(strncmp(argv[3],"2",1)==0)
		{
			cout<<" Running test 2: sample trajectory 2 (peak)"<<endl;
			int N_l = 50000;
			int N_r = 3000;
			double ctot = 2.7;
			double beta = 16.8;
			double koff = 0.01;
			double endTime = 10000;
			cout<<"simming for "<< endTime <<" seconds "<<endl;
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false,true,false);
			delete s;
			return;
		}
		else if(strncmp(argv[3],"3",1)==0)
		{
			cout<<" Running test 3: final distribution before PT"<<endl;
			int N_l = 42000;
			int N_r = 3000;
			double koff = 0.01;
			double ctot = 0.378;
			double beta = 0.3;
			double endTime = 1200;
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false, false, true);
			delete s;
			return;
		}
		else if(strncmp(argv[3],"4",1)==0)
		{
			cout<<" Running test 4: Observables below PT"<<endl;
			
			//Below PT
			int N_l = 42000;
			int N_r = 3000;
			double ctot = 0.378;
			double beta = 0.3;
			double koff = 0.01;
			double endTime = 2000;
			
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, true, false, false);
			delete s;
			return;
		}
		else if(strncmp(argv[3],"5",1)==0)
		{
			cout<<" Running test 5: final distribution magic"<<endl;
			int N_l = 1000;
			int N_r = 3000;
			double koff = 0.01;
			double ctot = 0.054;
			double beta = 16.8;
			double endTime = 2000;
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false, false, true);
			delete s;
			return;
		}
		else if(strncmp(argv[3],"6",1)==0)
		{
			cout<<" Running test 6: final distribution above PT"<<endl;
			int N_l = 42000;
			int N_r = 3000;
			double koff = 0.01;
			double ctot = 0.378;
			double beta = 90;
			double endTime = 1000;
			System *s = runTLBRSystem(N_l, N_r, ctot, beta, koff, filename, endTime, false, false, true);
			delete s;
			return;
		}
	}
	else
	{
		cout<<"You didn't specify which TLBR setup to run."<<endl;
		cout<<"  Your options:  1) Sample Trajectory 1 (medium)"<<endl;
		cout<<"                 2) Sample Trajectory 2 (long)"<<endl;
		cout<<"                 3) Final distribution before PT (fast)"<<endl;
		cout<<"                 4) Final distribution after PT (very long)"<<endl;
		return;
	}
}





MoleculeType * TLBR_createL(System * s, int count);
MoleculeType * TLBR_createR(System * s, int count);

void TLBR_createFreeBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate);
void TLBR_createUnBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate);
void TLBR_createCrossLinkingRxns(System * s, MoleculeType * L, MoleculeType *R, double rate);

void addObservables(System * s, MoleculeType * L, MoleculeType *R);

System * runTLBRSystem(int N_l, int N_r, double cTot, double beta, double koff, char * outputFileName,double endTime, bool outputObs, bool outputAvg, bool outputFinalDist)
{
	
	//cout.setf(ios::scientific);
	
	double freeBindRate = (cTot*koff)/(3.0*(double)N_l);
	double crossLinkRate = (beta*koff)/(double)N_r;
	
//	
//	cout<<"cTot: " <<cTot<<endl;
//	cout<<"koff: " << koff << endl;
//	cout<<"N_L: "  << N_l << endl;
//	cout<<"free Rate: "<<freeBindRate<<endl;
//	cout<<beta<<endl;
//	cout<<crossLinkRate<<endl;
//	exit(0);
	
	System * s = new System("TLBR",true);
	MoleculeType * L = TLBR_createL(s,N_l);
	MoleculeType * R = TLBR_createR(s,N_r);
	
	
	TLBR_createFreeBindingRxns(s,L,R,freeBindRate);
	TLBR_createUnBindingRxns(s,L,R,koff);
	TLBR_createCrossLinkingRxns(s,L,R,crossLinkRate);
	
	if(outputObs) addObservables(s,L,R);
	
	s->prepareForSimulation();
	
	
	/////
	s->registerOutputFileLocation(outputFileName);
	if(outputObs || outputAvg) s->outputAllObservableNames();
	
	double simTime = s->getCurrentTime();
	double dt = 100;
	double nextStoppingTime = dt;
	double avg=-1;
	if(outputAvg) avg = s->outputMeanCount(R);
	
	if(outputObs) s->outputAllObservableCounts();
	while(simTime<endTime)
	{
		simTime = s->stepTo(nextStoppingTime);
		if(outputAvg) avg = s->outputMeanCount(R);
		cout<<"at time: "<<simTime<<"  avg: "<<avg<<endl;
		
		if(outputObs) s->outputAllObservableCounts();
		nextStoppingTime += dt;
	}
	//double avg = s->outputMeanCount(R);
	s->printAllReactions();
	
	avg = s->calculateMeanCount(R);
	cout<<"Finished...  avg:"<<avg<<endl;
	//s->printAllComplexes();
	
	//If we aren't outputting the trajectory, make sure that we output the
	//final aggregate size distribution
	if(outputFinalDist) s->outputMoleculeTypeCountPerComplex(R);
	if(outputFinalDist) s->outputMoleculeTypeCountPerComplex(R);
	
	return s;
}

MoleculeType * TLBR_createL(System * s, int count)
{
	int numOfBsites = 3;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "r0";
	bSiteNames[1] = "r1";
	bSiteNames[2] = "r2";
	int numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *L = new MoleculeType("L",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}

MoleculeType * TLBR_createR(System * s, int count)
{
	int numOfBsites = 2;
	const char ** bSiteNames = new const char * [numOfBsites];
	bSiteNames[0] = "l0";
	bSiteNames[1] = "l1";
	int numOfStates = 0;
	const char ** stateNames = new const char * [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *R = new MoleculeType("R",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	R->populateWithDefaultMolecules(count);
	return R;	
}



void TLBR_createFreeBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate)
{
	// Reaction r0 binds l0
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	ReactionClass *r = new ReactionSimpleBinding("FreeBinding(r0-l0)",n_reactants,reactantTemplates,rate,"r0","l0");
	s->addReaction(r);
	
	// Reaction r1 binds l0
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("FreeBinding(r1-l0)",n_reactants,reactantTemplates,rate,"r1","l0");
	s->addReaction(r);
	
	// Reaction r2 binds l0
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("FreeBinding(r2-l0)",n_reactants,reactantTemplates,rate,"r2","l0");
	s->addReaction(r);
	
	// Reaction r0 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("FreeBinding(r0-l1)",n_reactants,reactantTemplates,rate,"r0","l1");
	s->addReaction(r);
	
	// Reaction r1 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("FreeBinding(r1-l1)",n_reactants,reactantTemplates,rate,"r1","l1");
	s->addReaction(r);
	
	// Reaction r2 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("FreeBinding(r2-l1)",n_reactants,reactantTemplates,rate,"r2","l1");
	s->addReaction(r);
}


void TLBR_createUnBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate)
{
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(R);
	reactantTemplates[0]->addOccupiedBindingSite("l0");
	ReactionClass *r = new ReactionUnbinding("Unbind(l0)",n_reactants,reactantTemplates,rate,"l0");
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	n_reactants = 1;
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(R);
	reactantTemplates[0]->addOccupiedBindingSite("l1");
	r = new ReactionUnbinding("Unbind(l1)",n_reactants,reactantTemplates,rate,"l1");
	r->setTraversalLimit(2);
	s->addReaction(r);
}




void TLBR_createCrossLinkingRxns(System * s, MoleculeType * L, MoleculeType *R, double rate)
{
	int traversalLimit = 2;
	bool doNotAllowSameComplex = true;


	///////// r0 binds l0
	
	/////////// 1
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	ReactionClass *r = new ReactionSimpleBinding("CrossLink(r0-l0)",n_reactants,reactantTemplates,rate,"r0","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	/////////// 2
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r0-l0)",n_reactants,reactantTemplates,rate,"r0","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	/////////// 3
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r0-l0)",n_reactants,reactantTemplates,rate,"r0","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	
	///////// r1 binds l0
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r1-l0)",n_reactants,reactantTemplates,rate,"r1","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r1-l0)",n_reactants,reactantTemplates,rate,"r1","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r1-l0)",n_reactants,reactantTemplates,rate,"r1","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	///////// r2 binds l0
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r2-l0)",n_reactants,reactantTemplates,rate,"r2","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r2-l0)",n_reactants,reactantTemplates,rate,"r2","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l0");
	r = new ReactionSimpleBinding("CrossLink(r2-l0)",n_reactants,reactantTemplates,rate,"r2","l0");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	///////////
	//////////
	////////
	
	///////// r0 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r0-l1)",n_reactants,reactantTemplates,rate,"r0","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r0-l1)",n_reactants,reactantTemplates,rate,"r0","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r0-l1)",n_reactants,reactantTemplates,rate,"r0","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	///////// r1 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r1-l1)",n_reactants,reactantTemplates,rate,"r1","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r1-l1)",n_reactants,reactantTemplates,rate,"r1","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r1-l1)",n_reactants,reactantTemplates,rate,"r1","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	///////// r2 binds l1
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addEmptyBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r2-l1)",n_reactants,reactantTemplates,rate,"r2","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addEmptyBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r2-l1)",n_reactants,reactantTemplates,rate,"r2","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(L);
	reactantTemplates[0]->addOccupiedBindingSite("r0");
	reactantTemplates[0]->addOccupiedBindingSite("r1");
	reactantTemplates[0]->addEmptyBindingSite("r2");
	reactantTemplates[1] = new TemplateMolecule(R);
	reactantTemplates[1]->addEmptyBindingSite("l1");
	r = new ReactionSimpleBinding("CrossLink(r2-l1)",n_reactants,reactantTemplates,rate,"r2","l1");
	if(doNotAllowSameComplex) r->setDoNotAllowSameComplex();
	r->setTraversalLimit(traversalLimit);
	s->addReaction(r);
}



void addObservables(System * s, MoleculeType * L, MoleculeType *R)
{
//	TemplateMolecule *tR = new TemplateMolecule(R);
//	tR->addEmptyBindingSite("l0");
//	tR->addEmptyBindingSite("l1");
//	Observable * obsR = new Observable("R_free",tR);
//	R->addObservable(obsR);
//	
//	tR = new TemplateMolecule(R);
//	obsR = new Observable("R_tot",tR);
//	R->addObservable(obsR);
	
	
	//////////
	TemplateMolecule *tL = new TemplateMolecule(L);
	tL->addEmptyBindingSite("r0");
	tL->addEmptyBindingSite("r1");
	tL->addEmptyBindingSite("r2");
	Observable * obsL = new Observable("L(-,-,-)",tL);
	L->addObservable(obsL);
	
	/////////////
	tL = new TemplateMolecule(L);
	tL->addOccupiedBindingSite("r0");
	tL->addEmptyBindingSite("r1");
	tL->addEmptyBindingSite("r2");
	obsL = new Observable("L(+,-,-)",tL);
	L->addObservable(obsL);
	
	tL = new TemplateMolecule(L);
	tL->addEmptyBindingSite("r0");
	tL->addOccupiedBindingSite("r1");
	tL->addEmptyBindingSite("r2");
	obsL = new Observable("L(-,+,-)",tL);
	L->addObservable(obsL);
	
	tL = new TemplateMolecule(L);
	tL->addEmptyBindingSite("r0");
	tL->addEmptyBindingSite("r1");
	tL->addOccupiedBindingSite("r2");
	obsL = new Observable("L(-,-,+)",tL);
	L->addObservable(obsL);
	
	
	/////////////
	tL = new TemplateMolecule(L);
	tL->addOccupiedBindingSite("r0");
	tL->addOccupiedBindingSite("r1");
	tL->addEmptyBindingSite("r2");
	obsL = new Observable("L(+,+,-)",tL);
	L->addObservable(obsL);
	
	tL = new TemplateMolecule(L);
	tL->addOccupiedBindingSite("r0");
	tL->addEmptyBindingSite("r1");
	tL->addOccupiedBindingSite("r2");
	obsL = new Observable("L(+,-,+)",tL);
	L->addObservable(obsL);
	
	tL = new TemplateMolecule(L);
	tL->addEmptyBindingSite("r0");
	tL->addOccupiedBindingSite("r1");
	tL->addOccupiedBindingSite("r2");
	obsL = new Observable("L(-,+,+)",tL);
	L->addObservable(obsL);
	
	
	
	//////////////
	tL = new TemplateMolecule(L);
	tL->addOccupiedBindingSite("r0");
	tL->addOccupiedBindingSite("r1");
	tL->addOccupiedBindingSite("r2");
	obsL = new Observable("L(+,+,+)",tL);
	L->addObservable(obsL);
	
	
	
	
	
	
}



