#include "tlbr.hh"



#include <math.h>


using namespace NFcore;
using namespace NFtest_tlbr;


void NFtest_tlbr::run(map<string,string> &argMap)
{
	/*
	 * To run TLBR
	 * 
	 * ./NFsim -test tlbr -preset 1
	 * 
	 * -n_L = starting number of ligand molecules
	 * -n_R = starting number of receptor molecules
	 * -cTot 
	 * -beta
	 * -koff
	 * -sim = length of run
	 * 
	 * 
	 * 
	 */
	
	int preset=-1;
	preset=NFinput::parseAsInt(argMap,"preset",preset);
	
	//Set up the default Parameters
	int n_L = 42000;
	int n_R = 3000;
	double koff = 0.01;
	double cTot = 0.378;
	double beta = 0.3;
	double simTime = 3000;
	double dt = 100;
	
	bool outputObservables = false;
	bool outputAvg = true;
	bool outputFinalDist = false;
	
	//If a valid preset was entered, then set the values
	//to the requested values
	if(preset==-1) {}  //This means we are not using any preset values...
	else if(preset==1) {
		cout<<"Loading preset parameters 1: Ramping Trajectory"<<endl;
		n_L = 2000;
		n_R = 3000;
		cTot = 0.11;
		beta = 16.8;
		koff = 0.01;
		simTime = 10000;
		outputObservables = false;
		outputAvg = true;
		outputFinalDist = false;
	}
	else if(preset==2) {
		cout<<"Loading preset parameters 2: Peaking Trajectory"<<endl;
		n_L = 50000;
		n_R = 3000;
		cTot = 2.7;
		beta = 16.8;
		koff = 0.01;
		simTime = 10000;
		outputObservables = false;
		outputAvg = true;
		outputFinalDist = false;
	} else if(preset==3) {
		cout<<"Loading preset parameters 3: Final Distribution before PT"<<endl;
		n_L = 42000;
		n_R = 3000;
		koff = 0.01;
		cTot = 0.378;
		beta = 0.3;
		simTime = 1200;
		outputObservables = false;
		outputAvg = false;
		outputFinalDist = true;
	} else if(preset==4) {
		cout<<"Loading preset parameters 4: Final Distribution above PT"<<endl;
		n_L = 42000;
		n_R = 3000;
		koff = 0.01;
		cTot = 0.378;
		beta = 90;
		simTime = 1000;
		outputObservables = false;
		outputAvg = false;
		outputFinalDist = true;
	} else if(preset==5) {
		cout<<"Loading preset parameters 5: Final Distribution magic"<<endl;
		n_L = 1000;
		n_R = 3000;
		koff = 0.01;
		cTot = 0.054;
		beta = 16.8;
		simTime = 2000;
		outputObservables = false;
		outputAvg = false;
		outputFinalDist = true;
	}
	else {
		cout<<"!! Warning: The preset value of "<<preset<<" that you gave me does nothing!"<<endl;
	}
	
	
	n_L=NFinput::parseAsInt(argMap,"n_L",n_L);
	n_R=NFinput::parseAsInt(argMap,"n_R",n_R);
	koff=NFinput::parseAsDouble(argMap,"koff",koff);
	cTot=NFinput::parseAsDouble(argMap,"cTot",cTot);
	beta=NFinput::parseAsDouble(argMap,"beta",beta);
	simTime=NFinput::parseAsDouble(argMap,"sim",simTime);
	dt=NFinput::parseAsDouble(argMap,"dt",dt);
	
	
	string filename="testTlbrOut";
	if(argMap.find("out")!=argMap.end()) {
		filename = argMap.find("out")->second;
	}
	if(filename.empty()) {
		cout<<"No filename given: using testTlbrOut[#]_nf.out"<<endl;
	}
	
	
	
	int runs=NFinput::parseAsInt(argMap,"runs",1);
	string sRuns;
	std::stringstream out;
	out << runs;
	sRuns = out.str();
	
	
	
	for(int r=0; r<runs; r++) {
		
		string sRuns;
		std::stringstream out;
		out << r;
		sRuns = out.str();
	
		cout<<"Executing Run Number: "<<r<<endl<<endl;
		runSystem(n_L, n_R, cTot, beta, koff,
				filename+sRuns+string("_nf.out"), simTime, dt,
				outputObservables, outputAvg, outputFinalDist);
		
	}
	
	
/*	
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
		cout<<"  Your options:  -r 1: Sample Trajectory 1 (medium)"<<endl;
		cout<<"                 -r 2: Sample Trajectory 2 (long)"<<endl;
		cout<<"                 -r 3: Final distribution before PT (fast)"<<endl;
		cout<<"                 -r 4: Final distribution after PT (very long)"<<endl;
		return;
	}
	*/
}


void NFtest_tlbr::runSystem(int n_L, int n_R, double cTot, double beta, double koff,
		string outputFileName, double simTime, double dt,
		bool outputObservables, bool outputAvg, bool outputFinalDist)
{
	double freeBindRate = (cTot*koff)/(3.0*(double)n_L);
	double crossLinkRate = (beta*koff)/(double)n_R;
	
	//Create the system with ligands, receptors, and reactions
	System * s = new System("TLBR",true);
	MoleculeType * L = createL(s,n_L);
	MoleculeType * R = createR(s,n_R);
	
	createFreeBindingRxns(s,L,R,freeBindRate);
	createUnbindingRxns(s,R,koff);
	createCrossLinkingRxns(s,L,R,crossLinkRate);
	
	//add the observables if we so choose
	//if(outputObservables) addObservables(s,L,R);
		
	//Prepare to run

	s->printAllMoleculeTypes();
	s->prepareForSimulation();
	s->registerOutputFileLocation(outputFileName.c_str());
	if(outputObservables || outputAvg) s->outputAllObservableNames();
	
	
	double currentTime = 0;
	double nextStoppingTime = dt;
	double avg=-1;
	if(outputAvg) avg = s->outputMeanCount(R);
		
	if(outputObservables) s->outputAllObservableCounts();
	int step=0;
	while(currentTime<simTime)
	{
		currentTime = s->stepTo(nextStoppingTime);
		if(outputAvg) avg = s->outputMeanCount(R);
		if(step%((int)(100/dt))==0)cout<<"at time: "<<currentTime<<"  avg: "<<avg<<endl;
			
		if(outputObservables) s->outputAllObservableCounts();
		nextStoppingTime += dt;
		step++;
	}
	
	s->printAllReactions();
		
	avg = s->calculateMeanCount(R);
	cout<<"Finished...  avg complex size: "<<avg<<" receptors"<<endl;
	//s->printAllComplexes();
		
	//If we aren't outputting the trajectory, make sure that we output the
	//final aggregate size distribution
	if(outputFinalDist) s->outputMoleculeTypeCountPerComplex(R);
	if(outputFinalDist) s->outputMoleculeTypeCountPerComplex(R);
	
	delete s;
}





MoleculeType * NFtest_tlbr::createL(System * s, int count)
{
	int numOfBsites = 3;
	string* bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "r0";
	bSiteNames[1] = "r1";
	bSiteNames[2] = "r2";
	int numOfStates = 0;
	string* stateNames = new string [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *L = new MoleculeType("L",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	L->populateWithDefaultMolecules(count);
	return L;	
}

MoleculeType * NFtest_tlbr::createR(System * s, int count)
{
	int numOfBsites = 2;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "l0";
	bSiteNames[1] = "l1";
	int numOfStates = 0;
	string * stateNames = new string [numOfStates];
	int * stateValues = new int [numOfStates];
	MoleculeType *R = new MoleculeType("R",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	R->populateWithDefaultMolecules(count);
	return R;	
}









//Create reactions where a free ligand binds a receptor
void NFtest_tlbr::createFreeBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate)
{
	
	{ // Reaction r0 binds l0
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r0", rTemp, "l0");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r0-l0)",rate,ts));
	}
	
	{ // Reaction r1 binds l0
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r1", rTemp, "l0");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r1-l0)",rate,ts));
	}
	
	{ // Reaction r2 binds l0
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r2", rTemp, "l0");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r2-l0)",rate,ts));
	}
	
	{ // Reaction r0 binds l1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l1");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r0", rTemp, "l1");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r0-l1)",rate,ts));
	}
	
	{ // Reaction r1 binds l1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l1");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r1", rTemp, "l1");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r1-l1)",rate,ts));
	}
	
	{ // Reaction r2 binds l1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l1");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addBindingTransform(lTemp,"r2", rTemp, "l1");
		ts->finalize();
		
		s->addReaction(new BasicRxnClass("FreeBinding(r2-l1)",rate,ts));
	}
}


void NFtest_tlbr::createUnbindingRxns(System * s, MoleculeType * R, double rate)
{
	{ // Unbind l0
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addOccupiedBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addUnbindingTransform(rTemp,"l0");
		ts->finalize();
		
		ReactionClass *r = new BasicRxnClass("Unbind(l0)",rate,ts);
		r->setTraversalLimit(2);
		s->addReaction(r);
	}
	{ // Unbind l1
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addOccupiedBindingSite("l1");
		
		vector <TemplateMolecule *> templates; templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		ts->addUnbindingTransform(rTemp,"l1");
		ts->finalize();
		
		ReactionClass *r = new BasicRxnClass("Unbind(l1)",rate,ts);
		r->setTraversalLimit(2);
		s->addReaction(r);
	}
}




void NFtest_tlbr::createCrossLinkingRxns(System * s, MoleculeType * L, MoleculeType *R, double rate)
{
	int traversalLimit = 2;
	bool doNotAllowSameComplex = true;


	///////// r0 binds l0 ////////////////////////////////////////
	{/////////// Variant 1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addOccupiedBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 1)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 2
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addOccupiedBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 2)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 3
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addOccupiedBindingSite("r1");
		lTemp->addOccupiedBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 3)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	
	
	///////// r1 binds l0 ////////////////////////////////////////
	{/////////// Variant 1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addOccupiedBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 1)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 2
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addOccupiedBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 2)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 3
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addOccupiedBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addOccupiedBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 3)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	
	
	///////// r2 binds l0 ////////////////////////////////////////
	{/////////// Variant 1
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addOccupiedBindingSite("r0");
		lTemp->addEmptyBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 1)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 2
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addEmptyBindingSite("r0");
		lTemp->addOccupiedBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 2)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	{/////////// Variant 3
		TemplateMolecule *lTemp = new TemplateMolecule(L);
		lTemp->addOccupiedBindingSite("r0");
		lTemp->addOccupiedBindingSite("r1");
		lTemp->addEmptyBindingSite("r2");
		
		TemplateMolecule *rTemp = new TemplateMolecule(R);
		rTemp->addEmptyBindingSite("l0");
		
		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
		
		TransformationSet *ts = new TransformationSet(templates);
		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
		ts->finalize();
		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 3)",rate,ts);
		r->setTraversalLimit(traversalLimit);
		s->addReaction(r);
	}
	
	////////////////////////////////////////////////////////////////////////
	
	
	///////// r0 binds l1 ////////////////////////////////////////
		{/////////// Variant 1
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addEmptyBindingSite("r0");
			lTemp->addOccupiedBindingSite("r1");
			lTemp->addEmptyBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 1)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 2
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addEmptyBindingSite("r0");
			lTemp->addEmptyBindingSite("r1");
			lTemp->addOccupiedBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 2)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 3
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addEmptyBindingSite("r0");
			lTemp->addOccupiedBindingSite("r1");
			lTemp->addOccupiedBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 3)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		
		
		///////// r1 binds l1 ////////////////////////////////////////
		{/////////// Variant 1
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addEmptyBindingSite("r0");
			lTemp->addEmptyBindingSite("r1");
			lTemp->addOccupiedBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 1)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 2
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addOccupiedBindingSite("r0");
			lTemp->addEmptyBindingSite("r1");
			lTemp->addEmptyBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 2)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 3
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addOccupiedBindingSite("r0");
			lTemp->addEmptyBindingSite("r1");
			lTemp->addOccupiedBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 3)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		
		
		///////// r2 binds l1 ////////////////////////////////////////
		{/////////// Variant 1
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addOccupiedBindingSite("r0");
			lTemp->addEmptyBindingSite("r1");
			lTemp->addEmptyBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 1)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 2
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addEmptyBindingSite("r0");
			lTemp->addOccupiedBindingSite("r1");
			lTemp->addEmptyBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 2)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
		{/////////// Variant 3
			TemplateMolecule *lTemp = new TemplateMolecule(L);
			lTemp->addOccupiedBindingSite("r0");
			lTemp->addOccupiedBindingSite("r1");
			lTemp->addEmptyBindingSite("r2");
			
			TemplateMolecule *rTemp = new TemplateMolecule(R);
			rTemp->addEmptyBindingSite("l1");
			
			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
			
			TransformationSet *ts = new TransformationSet(templates);
			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
			ts->finalize();
			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 3)",rate,ts);
			r->setTraversalLimit(traversalLimit);
			s->addReaction(r);
		}
}










/*



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

*/

