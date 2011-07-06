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
	double cTot = 0.84; //0.378;
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
	}else if(preset==6) {
		cout<<"Loading preset parameters 5: Speed"<<endl;
		n_L = 42000;
		n_R = 3000;
		koff = 0.01;
		cTot = 0.378;
		beta = 90;
		simTime = 1000;
		outputObservables = false;
		outputAvg = false;
		outputFinalDist = false;
	}else if(preset==7) {
		cout<<"Loading preset parameters 5: Memory"<<endl;
		n_L = 42000*23;
		n_R = 3000*23;
		koff = 0.01;
		cTot = 0.378;
		beta = 90;
		simTime = 0.01;
		outputObservables = false;
		outputAvg = false;
		outputFinalDist = false;
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

}


void NFtest_tlbr::runSystem(int n_L, int n_R, double cTot, double beta, double koff,
		string outputFileName, double simTime, double dt,
		bool outputObservables, bool outputAvg, bool outputFinalDist)
{
	double freeBindRate = (cTot*koff)/(3.0*(double)n_L);
	double crossLinkRate = (beta*koff)/(double)n_R;

	cout<<"beta: " <<beta<<endl;
	cout<<"cTot: " <<cTot<<endl;
	cout<<"kOff: " <<koff<<endl;
	cout<<"N_l:  " <<n_L<<endl;
	cout<<"N_r:  " <<n_R<<endl;


	//Create the system with ligands, receptors, and reactions
	System * s = new System("TLBR",true);
	vector<vector<string> > v;
	MoleculeType * L = createL(s,n_L);
	L->addEquivalentComponents(v);
	MoleculeType * R = createR(s,n_R);
	R->addEquivalentComponents(v);


	createUnbindingRxns(s,R,koff);
	createFreeBindingRxns(s,L,R,freeBindRate);
	createCrossLinkingRxns(s,L,R,crossLinkRate);




	//add the observables if we so choose
	//if(outputObservables) addObservables(s,L,R);

	//Prepare to run

	s->printAllMoleculeTypes();
	s->prepareForSimulation();
	s->registerOutputFileLocation(outputFileName.c_str());

	/////////////////  added here for calculating scalling...
	//s->equilibriate(200,20);
	s->sim(simTime,3);

	string x = "";
	cout<<"waiting for you to enter something."<<endl;
	cin>>x;

	s->printAllReactions();
	// NETGEN -- redirected call to the ComplexList object at s->allComplexes
	cout << "Finished...  avg complex size: " << (s->getAllComplexes()).calculateMeanCount(R) << " receptors" << endl;

	delete s;
	return;
	/////////////////////////////////////

	if(outputObservables || outputAvg) s->outputAllObservableNames();


	double currentTime = 0;
	double nextStoppingTime = dt;
	double avg=-1;
	// NETGEN -- redirected call to the ComplexList object at s->allComplexes
	if(outputAvg) avg = (s->getAllComplexes()).outputMeanCount(R);

	if(outputObservables) s->outputAllObservableCounts();
	int step=0;
	while(currentTime<simTime)
	{
		currentTime = s->stepTo(nextStoppingTime);
		if(outputAvg) avg = (s->getAllComplexes()).outputMeanCount(R);
		if(step%((int)(100/dt))==0)cout<<"at time: "<<currentTime<<"  avg: "<<avg<<endl;

		if(outputObservables) s->outputAllObservableCounts();
		nextStoppingTime += dt;
		step++;
	}

	s->printAllReactions();

	// NETGEN -- redirected call to the ComplexList object at s->allComplexes
	avg = (s->getAllComplexes()).calculateMeanCount(R);
	cout<<"Finished...  avg complex size: "<<avg<<" receptors"<<endl;
	//(s->allComplexes).printAllComplexes();

	//If we aren't outputting the trajectory, make sure that we output the
	//final aggregate size distribution
	// NETGEN -- redirected call to the ComplexList object at s->allComplexes
	if(outputFinalDist) (s->getAllComplexes()).outputMoleculeTypeCountPerComplex(R);
	if(outputFinalDist) (s->getAllComplexes()).outputMoleculeTypeCountPerComplex(R);

	delete s;
}





MoleculeType * NFtest_tlbr::createL(System * s, int count)
{
	vector <string> compName;
	vector <string> defaultCompState;
	vector < vector <string> > possibleCompStates;

	compName.push_back("r0");
	defaultCompState.push_back("No State");
	vector <string> possibleR0states;
	possibleCompStates.push_back(possibleR0states);

	compName.push_back("r1");
	defaultCompState.push_back("No State");
	vector <string> possibleR1states;
	possibleCompStates.push_back(possibleR1states);

	compName.push_back("r2");
	defaultCompState.push_back("No State");
	vector <string> possibleR2states;
	possibleCompStates.push_back(possibleR2states);

	MoleculeType *L = new MoleculeType("L", compName, defaultCompState, possibleCompStates, s);
	L->populateWithDefaultMolecules(count);
	return L;
}

MoleculeType * NFtest_tlbr::createR(System * s, int count)
{
	vector <string> compName;
	vector <string> defaultCompState;
	vector < vector <string> > possibleCompStates;

	compName.push_back("l0");
	defaultCompState.push_back("No State");
	vector <string> possibleL0states;
	possibleCompStates.push_back(possibleL0states);

	compName.push_back("l1");
	defaultCompState.push_back("No State");
	vector <string> possibleL1states;
	possibleCompStates.push_back(possibleL1states);

	MoleculeType *R = new MoleculeType("R", compName, defaultCompState, possibleCompStates, s);
	R->populateWithDefaultMolecules(count);
	return R;
}









//Create reactions where a free ligand binds a receptor
void NFtest_tlbr::createFreeBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate)
{
//	{ // Reaction r0 binds l0
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r0", rTemp, "l0");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r0-l0)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//
//	{ // Reaction r1 binds l0
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r1", rTemp, "l0");
//
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r1-l0)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//
//	{ // Reaction r2 binds l0
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r2", rTemp, "l0");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r2-l0)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//
//	{ // Reaction r0 binds l1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l1");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r0", rTemp, "l1");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r0-l1)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//
//	{ // Reaction r1 binds l1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l1");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r1", rTemp, "l1");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r1-l1)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//
//	{ // Reaction r2 binds l1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l1");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp ); templates.push_back( lTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addBindingTransform(lTemp,"r2", rTemp, "l1");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("FreeBinding(r2-l1)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
}


void NFtest_tlbr::createUnbindingRxns(System * s, MoleculeType * R, double rate)
{

//	{ // Unbind l0
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addOccupiedBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addUnbindingTransform(rTemp,"l0",NULL, "");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("Unbind(l0)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
//	{ // Unbind l1
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addOccupiedBindingSite("l1");
//
//		vector <TemplateMolecule *> templates; templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		ts->addUnbindingTransform(rTemp,"l1", NULL, "");
//		ts->finalize();
//
//		ReactionClass *r = new BasicRxnClass("Unbind(l1)",rate,ts);
//		r->setTraversalLimit(2);
//		s->addReaction(r);
//	}
}



void NFtest_tlbr::createCrossLinkingRxns(System * s, MoleculeType * L, MoleculeType *R, double rate)
{
//	int traversalLimit = 2;
//	bool doNotAllowSameComplex = true;
//
//	///////// r0 binds l0 ////////////////////////////////////////
//	{/////////// Variant 1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addOccupiedBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 1)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 2
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addOccupiedBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 2)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 3
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addOccupiedBindingSite("r1");
//		lTemp->addOccupiedBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r0",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r0-l0, 3)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//
//
//	///////// r1 binds l0 ////////////////////////////////////////
//	{/////////// Variant 1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addOccupiedBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 1)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 2
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addOccupiedBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 2)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 3
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addOccupiedBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addOccupiedBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r1",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r1-l0, 3)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//
//
//	///////// r2 binds l0 ////////////////////////////////////////
//	{/////////// Variant 1
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addOccupiedBindingSite("r0");
//		lTemp->addEmptyBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 1)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 2
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addEmptyBindingSite("r0");
//		lTemp->addOccupiedBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 2)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//	{/////////// Variant 3
//		TemplateMolecule *lTemp = new TemplateMolecule(L);
//		lTemp->addOccupiedBindingSite("r0");
//		lTemp->addOccupiedBindingSite("r1");
//		lTemp->addEmptyBindingSite("r2");
//
//		TemplateMolecule *rTemp = new TemplateMolecule(R);
//		rTemp->addEmptyBindingSite("l0");
//
//		vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//		TransformationSet *ts = new TransformationSet(templates);
//		if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l0"); }
//		else { ts->addBindingTransform(lTemp,"r2",rTemp,"l0"); }
//		ts->finalize();
//		ReactionClass *r=new BasicRxnClass("CrossLink(r2-l0, 3)",rate,ts);
//		r->setTraversalLimit(traversalLimit);
//		s->addReaction(r);
//	}
//
//	////////////////////////////////////////////////////////////////////////
//
//
//	///////// r0 binds l1 ////////////////////////////////////////
//		{/////////// Variant 1
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addEmptyBindingSite("r0");
//			lTemp->addOccupiedBindingSite("r1");
//			lTemp->addEmptyBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 1)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 2
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addEmptyBindingSite("r0");
//			lTemp->addEmptyBindingSite("r1");
//			lTemp->addOccupiedBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 2)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 3
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addEmptyBindingSite("r0");
//			lTemp->addOccupiedBindingSite("r1");
//			lTemp->addOccupiedBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r0",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r0",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r0-l1, 3)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//
//
//		///////// r1 binds l1 ////////////////////////////////////////
//		{/////////// Variant 1
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addEmptyBindingSite("r0");
//			lTemp->addEmptyBindingSite("r1");
//			lTemp->addOccupiedBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 1)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 2
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addOccupiedBindingSite("r0");
//			lTemp->addEmptyBindingSite("r1");
//			lTemp->addEmptyBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 2)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 3
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addOccupiedBindingSite("r0");
//			lTemp->addEmptyBindingSite("r1");
//			lTemp->addOccupiedBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r1",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r1",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r1-l1, 3)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//
//
//		///////// r2 binds l1 ////////////////////////////////////////
//		{/////////// Variant 1
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addOccupiedBindingSite("r0");
//			lTemp->addEmptyBindingSite("r1");
//			lTemp->addEmptyBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 1)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 2
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addEmptyBindingSite("r0");
//			lTemp->addOccupiedBindingSite("r1");
//			lTemp->addEmptyBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 2)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
//		{/////////// Variant 3
//			TemplateMolecule *lTemp = new TemplateMolecule(L);
//			lTemp->addOccupiedBindingSite("r0");
//			lTemp->addOccupiedBindingSite("r1");
//			lTemp->addEmptyBindingSite("r2");
//
//			TemplateMolecule *rTemp = new TemplateMolecule(R);
//			rTemp->addEmptyBindingSite("l1");
//
//			vector <TemplateMolecule *> templates; templates.push_back( lTemp ); templates.push_back( rTemp );
//
//			TransformationSet *ts = new TransformationSet(templates);
//			if(doNotAllowSameComplex) { ts->addBindingSeparateComplexTransform(lTemp,"r2",rTemp,"l1"); }
//			else { ts->addBindingTransform(lTemp,"r2",rTemp,"l1"); }
//			ts->finalize();
//			ReactionClass *r=new BasicRxnClass("CrossLink(r2-l1, 3)",rate,ts);
//			r->setTraversalLimit(traversalLimit);
//			s->addReaction(r);
//		}
}











