#ifdef __CDT_PARSER__
#define CLOCKS_PER_SEC
#endif


#include "NFcore.hh"

#include <math.h>
#include <fstream>
#include "../NFscheduler/NFstream.h"
#include "../NFscheduler/Scheduler.h"

#define ATOT_TOLERANCE 1e-9

using namespace std;
using namespace NFcore;

int System::NULL_EVENT_COUNTER = 0;


System::System(string name)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = false;     // NETGEN -- is this needed?
	// NETGEN
	allComplexes.setSystem( this );
	allComplexes.setUseComplex( false );

	this->outputGlobalFunctionValues=false;
	this->globalMoleculeLimit = 100000;
	rxnIndexMap=0;
	useBinaryOutput=false;
	outputEventCounter=false;
	globalEventCounter=0;
	onTheFlyObservables=true;
	universalTraversalLimit=-1;
	ds=0;
	selector = 0;
	csvFormat = false;
}


System::System(string name, bool useComplex)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;

	this->useComplex = useComplex;    // NETGEN -- is this needed?
	// NETGEN
	allComplexes.setSystem( this );
	allComplexes.setUseComplex( useComplex );

	this->outputGlobalFunctionValues=false;
	this->globalMoleculeLimit = 100000;

	rxnIndexMap=0;
	useBinaryOutput=false;
	onTheFlyObservables=true;
	outputEventCounter=false;
	globalEventCounter=0;
	universalTraversalLimit=-1;
	ds=0;
	selector = 0;
	csvFormat = false;
}

System::System(string name, bool useComplex, int globalMoleculeLimit)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = useComplex;  // NETGEN -- is this needed?
	// NETGEN
	allComplexes.setSystem( this );
	allComplexes.setUseComplex( useComplex );

	this->globalMoleculeLimit=globalMoleculeLimit;
	this->outputGlobalFunctionValues=false;

	rxnIndexMap=0;
	useBinaryOutput=false;
	outputEventCounter=false;
	globalEventCounter=0;
	onTheFlyObservables=true;
	universalTraversalLimit=-1;
	ds=0;
	selector = 0;
	csvFormat = false;
}



System::~System()
{
	if(ds!=0) delete ds;

	if(selector!=0) delete selector;

	//Delete the rxnIndexMap array
	if(rxnIndexMap!=NULL) {
		for(unsigned int r=0; r<allReactions.size(); r++)
			if(rxnIndexMap[r]!=NULL) { delete [] rxnIndexMap[r]; }
		delete [] rxnIndexMap;
	}

	//Need to delete reactions
	ReactionClass *r;
	while(allReactions.size()>0)
	{
		r = allReactions.back();
		allReactions.pop_back();
		delete r;
	}

	//Delete all observables of this type that exist
	Observable *o;
	while(obsToOutput.size()>0)
	{
		o = obsToOutput.back();
		obsToOutput.pop_back();
		delete o;
	}


	//Delete all MoleculeTypes (which deletes all molecules and templates)
	MoleculeType *s;
	while(allMoleculeTypes.size()>0)
	{
		s = allMoleculeTypes.back();
		allMoleculeTypes.pop_back();
		delete s;
	}

	// NETGEN -- not needed, complexList managed in its own class
	/*
	//Delete all the complexes
	Complex *c;
	while(allComplexes.size()>0)
	{
		c = allComplexes.back();
		allComplexes.pop_back();
		delete c;
	}
    */

	GlobalFunction *gf;
	while(this->globalFunctions.size()>0)
	{
		gf = globalFunctions.back();
		globalFunctions.pop_back();
		delete gf;
	}

	LocalFunction *lf;
	while(this->localFunctions.size()>0)
	{
		lf = localFunctions.back();
		localFunctions.pop_back();
		delete lf;
	}

	CompositeFunction *cf;
	while(this->compositeFunctions.size()>0)
	{
		cf = compositeFunctions.back();
		compositeFunctions.pop_back();
		delete cf;
	}


	nextReaction = 0;


	//Need to delete reactions
	Outputter *op;
	while(allOutputters.size()>0)
	{
		op = allOutputters.back();
		allOutputters.pop_back();
		delete op;
	}

	//Close our connections to output files
	outputFileStream.close();

	propensityDumpStream.close();
}


void System::setOutputToBinary()
{
	this->useBinaryOutput = true;
	if(outputFileStream.is_open()) {
		outputFileStream.close();
		cerr<<"Error!! You are trying to switch the output of this system to Binary, but\n";
		cerr<<"you already have an open file stream that is not binary!  The results are\n";
		cerr<<"therefore unpredictable.  It would be better if you fix this problem first.\n";
		cerr<<"This problem is caused when you call 'setOutputToBinary()' after you call\n";
		cerr<<"registerOutputFileLocation().\n";
		cerr<<"So I'm just going to stop now."<<endl;
		exit(1);
	}
}

void System::turnOff_OnTheFlyObs() {
	this->onTheFlyObservables=false;
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
		(*rxnIter)->turnOff_OnTheFlyObs();
}

int System::getNumOfSpeciesObs() const {
	return (int)speciesObservables.size();
}
Observable * System::getSpeciesObs(int index) const
{
	return speciesObservables.at(index);
}


void System::registerOutputFileLocation(string filename)
{
	if(outputFileStream.is_open()) { outputFileStream.close(); }
	if(useBinaryOutput) {
		outputFileStream.open((filename).c_str(), ios_base::out | ios_base::binary | ios_base::trunc);

		if(!outputFileStream.is_open()) {
			cerr<<"Error in System!  cannot open output stream to file "<<filename<<". "<<endl;
			cerr<<"quitting."<<endl;
			exit(1);
		}

		//ios_base::out -- Set for output only
		//ios_base::binary --  Set output to binary
		//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there

		//Also, output a header file to keep track of the number
		NFstream headerFile;
		int tabCount=0;
		headerFile.open((filename+".head").c_str());
		headerFile<<"#\tTime"; tabCount++;
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
			int oTot = (*molTypeIter)->getNumOfMolObs();
			for(int o=0; o<oTot; o++) {
				headerFile<<"\t"<<(*molTypeIter)->getMolObs(o)->getName();
				tabCount++;
			}
		}
		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ ) {
				headerFile<<"\t"<<(*functionIter)->getNiceName();
				tabCount++;
			}
		if(outputEventCounter)  headerFile<<"\tEventCounter";
		headerFile<<endl;
		for(int t=0; t<tabCount; t++) headerFile<<"\t";
		headerFile.close();

	} else {
		outputFileStream.open(filename.c_str());

		if(!outputFileStream.is_open()) {
			cerr<<"Error in System!  cannot open output stream to file "<<filename<<". "<<endl;
			cerr<<"quitting."<<endl;
			exit(1);
		}

		outputFileStream.setf(ios::scientific);
		outputFileStream.precision(8);
	}

}



void System::tagReaction(int rID) {

	if(rID<0 || rID>=this->allReactions.size() ) {
		cerr<<"!!! Error when trying to tag reaction with reaction ID "<<rID<<endl;
		cerr<<"!!! Reaction with that ID does not exist."<<endl;
		cerr<<"!!! quitting now."<<endl;
		exit(1);
	}
	allReactions.at(rID)->tag();


}


void System::addObservableForOutput(Observable *o) {
	if(o->getType()==Observable::SPECIES)
		this->speciesObservables.push_back(o);
	this->obsToOutput.push_back(o);
}


int System::addMoleculeType(MoleculeType *MoleculeType)
{
	allMoleculeTypes.push_back(MoleculeType);
	return (allMoleculeTypes.size()-1);
}


void System::addReaction(ReactionClass *reaction)
{
	if(this->universalTraversalLimit>0)
		reaction->setTraversalLimit(universalTraversalLimit);

	reaction->init();
	allReactions.push_back(reaction);
}

void System::addNecessaryUpdateReaction(ReactionClass *reaction)
{
	reaction->init();
	allReactions.push_back(reaction);
	necessaryUpdateRxns.push_back(reaction);
}

void System::setUniversalTraversalLimit(int utl) {
	this->universalTraversalLimit = utl;
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
		(*rxnIter)->setTraversalLimit(utl);

}



void System::addOutputter(Outputter *op) {
	this->allOutputters.push_back(op);
	op->outputHeader();
}
void System::dumpOutputters() {
	for(unsigned int i=0; i<allOutputters.size(); i++) {
		allOutputters.at(i)->output();
	}
}


void System::setDumpOutputter(DumpSystem *ds) {
	this->ds=ds;
}
void System::tryToDump() {
	if(ds!=0)
		ds->tryToDump(this->current_time);
}



bool System::addGlobalFunction(GlobalFunction *gf)
{
	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
	  	if(gf->getName()==(*functionIter)->getName()) return false;
	this->globalFunctions.push_back(gf);
	return true;
}




MoleculeType * System::getMoleculeTypeByName(string mName)
{
	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		//(*molTypeIter)->printDetails(); //<<endl;
		if((*molTypeIter)->getName()==mName)
		{
			return (*molTypeIter);
		}
	}
	cerr<<"!!! warning !!! cannot find molecule type name '"<< mName << "' in System: '"<<this->name<<"'"<<endl;
	exit(1);
	return 0;
}


Molecule * System::getMoleculeByUid(int uid)
{
	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		//(*molTypeIter)->printDetails(); //<<endl;
		for(int m=0; m<(*molTypeIter)->getMoleculeCount(); m++)
		{
				if( (*molTypeIter)->getMolecule(m)->getUniqueID() == uid)
					return (*molTypeIter)->getMolecule(m);
		}
	}
	cerr<<"!!! warning !!! cannot find active molecule with unique ID '"<< uid << "' in System: '"<<this->name<<"'"<<endl;
	return 0;
}

int System::getNumOfMolecules()
{
	int sum=0;
	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		sum+=(*molTypeIter)->getMoleculeCount();
	return sum;
}



int System::getMolObsCount(int moleculeTypeIndex, int observableIndex) const
{
	return allMoleculeTypes.at(moleculeTypeIndex)->getMolObsCount(observableIndex);
}






//When you are ready to run the simulation (meaning that all moleculeTypes
//all molecules, and all reactions have been created and registered with
//the system) call this function to populate all the reactant lists and
//observables.
void System::prepareForSimulation()
{
	this->selector = new DirectSelector(allReactions);

	cout<<"preparing simulation..."<<endl;
	//Note!!  : the order of preparing the system matters!  You have to prepare
	//some things before others, because certain things require other

	//First, set the observables up correctly, so when functions evaluate, they get the
	//correct values
	//for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
	//	(*molTypeIter)->addAllToObservables();
	//}

  	//First, we have to prep all the functions...
  	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
  		(*functionIter)->prepareForSimulation(this);

  	//cout<<"here 1..."<<endl;

  	for( int f=0; f<localFunctions.size(); f++)
  		localFunctions.at(f)->prepareForSimulation(this);

  	//cout<<"here 2..."<<endl;

  	for( int f=0; f<compositeFunctions.size(); f++)
  		compositeFunctions.at(f)->prepareForSimulation(this);

  	//cout<<"here 3..."<<endl;
    //this->printAllFunctions();

  	// now we prepare all reactions
	rxnIndexMap = new int * [allReactions.size()];
  	for(unsigned int r=0; r<allReactions.size(); r++)
  	{
  		rxnIndexMap[r] = new int[allReactions.at(r)->getNumOfReactants()];
  		allReactions.at(r)->setRxnId(r);
  	}

  	//cout<<"here 4..."<<endl;

	//This means we aren't going to add any more molecules to the system, so prep the rxns
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
		(*rxnIter)->prepareForSimulation();

	//cout<<"here 5..."<<endl;

	//If there are local functions to be had, make sure we set up those local function lists in the molecules
	//before we try to add molecules to reactant lists
	if(this->localFunctions.size()>0) {
	  	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	  		(*molTypeIter)->setUpLocalFunctionListForMolecules();
	}

	//cout<<"here 6..."<<endl;


  	//prep each molecule type for the simulation
  	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
  		(*molTypeIter)->prepareForSimulation();

  	//cout<<"here 7..."<<endl;

  	//add all the molecules to the appropriate observables
  	//NOT NECESSARY - molecules are added to observables when they are prepared
  	//for(obsIter=obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++)
  	//	(*obsIter)->clear();
  	//for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
  	//	(*molTypeIter)->addAllToObservables();
  	//}

  	//Add the complexes to Species observables
  	int match = 0;
  	for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++)
  	  	(*obsIter)->clear();

  	// NETGEN -- this bit replaces the commented block below
  	Complex * complex;
  	allComplexes.resetComplexIter();
  	while(  (complex = allComplexes.nextComplex()) )
  	{
  		if( complex->isAlive() )
  		{
  			for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++)
  			{
  				match = (*obsIter)->isObservable( complex );
  				for (int k=0; k<match; k++) (*obsIter)->straightAdd();
  			}
  		}
  	}
  	/*
  	for(complexIter = allComplexes.allComplexes.begin(); complexIter != allComplexes.end(); complexIter++) {
  		if((*complexIter)->isAlive()) {
  			for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++) {
  				match = (*obsIter)->isObservable((*complexIter));
  				for(int k=0; k<match; k++) (*obsIter)->straightAdd();
  			}
  		}
  	}
  	*/


  	//cout<<"here 8..."<<endl;




	//cout<<"here 9..."<<endl;





  	//if(BASIC_MESSAGE) cout<<"preparing the system...\n";
  	//printIndexAndNames();


//  if(go!=NULL)
// 	{
//  		go->writeGroupKeyFile();
// 		go->writeOutputFileHeader();
// 	}



	//finally, create the next reaction selector

	//this->selector = new LogClassSelector(allReactions);

	this->evaluateAllLocalFunctions();

  	recompute_A_tot();


}


void System::update_A_tot(ReactionClass *r, double old_a, double new_a)
{
	a_tot = selector->update(r,old_a,new_a);

	//BUILT IN DIRECT SEARCH
	//a_tot-=old_a;
	//a_tot+=new_a;
}


double System::recompute_A_tot()
{
	a_tot = selector->refactorPropensities();
	return a_tot;


//  BUILT IN DIRECT SEARCH
//	//Loop through the reactions and add up the rates
//	a_tot = 0;
//	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
//	{
//		a_tot += (*rxnIter)->update_a();
//		if(DEBUG) (*rxnIter)->printDetails();
//	}
//	return a_tot;
}






/* select the next reaction, given a_tot has been calculated */
double System::getNextRxn()
{
	nextReaction = 0;
	double x = selector->getNextReactionClass(nextReaction);
	if((int)x==-1) {
		this->printAllReactions();
		exit(1);
	}
	return selector->getNextReactionClass(nextReaction);


//  BUILT IN DIRECT SEARCH
//	double randNum = NFutil::RANDOM(a_tot);
//
//	double a_sum=0, last_a_sum=0;
//	nextReaction = 0;
//
//	//WARNING - DO NOT USE THE DEFAULT C++ RANDOM NUMBER GENERATOR FOR THIS STEP
//	// - IT INTRODUCES SMALL NUMERICAL ERRORS CAUSING THE ORDER OF RXNS TO
//	//   AFFECT SIMULATION RESULTS
//	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++)
//	{
//		a_sum += (*rxnIter)->get_a();
//		if (randNum <= a_sum && nextReaction==0)
//		{
//			nextReaction = (* rxnIter);
//			//cout<<"rNum: "<<randNum<<" last_a: "<<last_a_sum<<" a_sum "<<a_sum<<endl;
//			return (randNum-last_a_sum);
//		}
//		last_a_sum = a_sum;
//	}
//	cerr<<"Error: randNum exceeds a_sum!!!"<<endl;
//	cerr<<"randNum: "<<randNum<<"  a_sum: "<< a_sum<<" running a_tot:"<<a_tot<<endl;
//	return -1;
}


/* main simulation loop */
double System::sim(double duration, long int sampleTimes)
{
	return sim(duration,sampleTimes,true);
}


/* main simulation loop */
double System::sim(double duration, long int sampleTimes, bool verbose)
{
	System::NULL_EVENT_COUNTER=0;
	cout.setf(ios::scientific);
	cout<<"simulating system for: "<<duration<<" second(s)."<<endl;
	if(verbose) cout<<"\n";

	//First, output the header for the output of this simulation
	//outputAllObservableNames();
	//this->printAllReactions();


	//////////////////////////////
	clock_t start,finish;
	double time;
	start = clock();
	//////////////////////////////


	//Determine when to sample and print out initial setup
	double dSampleTime = duration / sampleTimes;
	double curSampleTime=current_time;

	//Do this once at the beginning, so that we start on the right page
	recompute_A_tot();

	double delta_t = 0; unsigned long long iteration = 0, stepIteration = 0;
	double end_time = current_time+duration;
	tryToDump();

	while(current_time<end_time)
	{
		//this->printAllObservableCounts(current_time);
		//2: Recompute a_tot for this time
		//cout<<" a_tot was : " << a_tot<<endl;
		//recompute_A_tot();
		//cout<<" a_tot (after recomputing) is : " << a_tot<<endl;

		//3: Select next reaction time (making sure we have something that can react)
		//   dt = -ln(rand) / a_tot;
		//Choose a random number on the OPEN interval (0,1) so that we never
		//have a dt=0 or a dt=infinity
		if(a_tot>ATOT_TOLERANCE) delta_t = -log(NFutil::RANDOM_OPEN()) / a_tot;
		else { delta_t=0; current_time=end_time; }
		if(DEBUG) cout<<"   Determine dt : " << delta_t << endl;


		//Report everything up until the next step if we have to
		if(DEBUG) cout<<"  Current Sample Time: "<<curSampleTime<<endl;
		if((current_time+delta_t)>=curSampleTime)
		{
			while((current_time+delta_t)>=(curSampleTime))
			{
				if(curSampleTime>end_time) break;
				outputAllObservableCounts(curSampleTime,globalEventCounter);
				//outputGroupData(curSampleTime);
				curSampleTime+=dSampleTime;
			}
			if(verbose) {
				cout << "Sim time: "           << (curSampleTime-dSampleTime);
				cout << "\tCPU time (total): " << ((double)(clock() - start)/(double)CLOCKS_PER_SEC) << "s";
				cout << "\t events (step): "   << stepIteration<<endl;
			}
			stepIteration=0;
			recompute_A_tot();
		}

		//cout<<"delta_t: " <<delta_t<<" atot: "<<a_tot<<endl;
		//Make sure we can react...
		if(delta_t==0) break;

		//4: Select next reaction class based on smallest j,
		//   such that sum of a_j over all j >= r2*a_tot
		double randElement = getNextRxn();
		//cout<<endl<<endl<<endl<<"-----------------------------------------------"<<endl;

		//cout<<"Fire: "<<nextReaction->getName()<<" at time "<< current_time<<endl;
		//Output selected reaction for debugging
		//cout<<"\nFiring: "<< endl;
		//nextReaction->printFullDetails();
		//cout<<endl<<endl;

		//this->printAllReactions();
		//this->printAllObservableCounts(this->current_time);
		//cout<<"\n";
		//Increment time
		iteration++;
		stepIteration++;
		globalEventCounter++;
		current_time+=delta_t;

		//5: Fire Reaction! (takes care of updates to lists and observables)
		nextReaction->fire(randElement);
		//this->printAllObservableCounts(this->current_time);
		//cout<<"\n---"<<endl;

		tryToDump();
		//	outputAllPropensities(current_time, nextReaction->getRxnId());

		//cout<<getObservableByName("Lig_free")->getCount()<<"/"<<getObservableByName("Lig_tot")->getCount()<<endl;
		//if(nextReaction->getName()=="Rule7") {
		//	cout<<getObservableByName("Lig_free")->getCount()<<"/"<<getObservableByName("Lig_tot")->getCount()<<endl;
		//	printAllReactions();
		//	exit(1);
		//}

		// TODO: debug!
		//this->getAllComplexes().printAllComplexes();
	}
	if(curSampleTime-dSampleTime<(end_time-0.5*dSampleTime)) {
		outputAllObservableCounts(curSampleTime,globalEventCounter);
	}


	finish = clock();
    time = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if(verbose) cout<<"\n";
    cout<<"   You just simulated "<< iteration <<" reactions in "<< time << "s\n";
    cout<<"   ( "<<((double)iteration)/time<<" reactions/sec, ";
    cout<<(time/((double)iteration))<<" CPU seconds/event )"<< endl;
    cout<<"   Null events: "<< System::NULL_EVENT_COUNTER;
    cout<<"   ("<<(time)/((double)iteration-(double)System::NULL_EVENT_COUNTER)<<" CPU seconds/non-null event )"<< endl;

	cout.unsetf(ios::scientific);
	return current_time;
}

double System::stepTo(double stoppingTime)
{
	double delta_t = 0;
	while(current_time<stoppingTime)
	{
		//2: Recompute a_tot for this time (this is not done here anymore!  reactions must
		//   be updated with the system as soon as a change to the propensity is made!
		//recompute_A_tot();

		//3: Select next reaction time (making sure we have something that can react)
		//   dt = -ln(rand) / a_tot;
		//Choose a random number on the closed interval (0,1) so that we never
		//have a dt=0 or a dt=infinity
		if(a_tot>ATOT_TOLERANCE) delta_t = -log(NFutil::RANDOM_CLOSED()) / a_tot;
		else
		{
			//Otherwise, we can't react for the rest of this step
			delta_t=0;
			current_time=stoppingTime;
			cout<<"Total propensity is zero, no further rxns can fire in this step."<<endl;
			break;
		}


		//Report everything up until the next step if we have to
		if((current_time+delta_t)>=stoppingTime)
		{
			//We are going to jump over the stopping time, so end the step
			break;
		}

		//4: Select next reaction class based on smallest j,
		//   such that sum of a_j over all j >= r2*a_tot
		double randElement = getNextRxn();


		//Increment time
		current_time+=delta_t;

		globalEventCounter++;

		//cout<<"Fire: "<<nextReaction->getName()<<" at time "<< current_time<<endl;

		//5: Fire Reaction! (takes care of updates to lists and observables)
		nextReaction->fire(randElement);
	}
	//cout<<"a_tot="<<a_tot;
	return current_time;
}

void System::singleStep()
{
	cout<<"  -System is at time: "<<this->current_time<<endl;
	double delta_t = 0;

	recompute_A_tot();
	cout<<"  -total propensity (a_total) calculated as: "<<a_tot<<endl;
	if(a_tot>ATOT_TOLERANCE) delta_t = -log(NFutil::RANDOM_CLOSED()) / a_tot;
	else
	{
		//Otherwise, we can't react for the rest of this step
		delta_t=0;
		cout<<"  -Total propensity is zero, no further rxns can fire."<<endl;
		return;
	}

	cout<<" -calculated time step is: "<<delta_t<<" seconds";
	double randElement = getNextRxn();

	//Increment time
	current_time+=delta_t;

	cout<<"  -Firing: "<<endl;
	nextReaction->printDetails();;

	//5: Fire Reaction! (takes care of updates to lists and observables)
	nextReaction->fire(randElement);
	cout<<"  -System time is now at time: "<<current_time<<endl;

	globalEventCounter++;
}

void System::equilibrate(double duration)
{
	double startTime = current_time;
	stepTo(duration);
	current_time = startTime;
}

void System::equilibrate(double duration, int statusReports)
{
	if(duration<=0) return;

	if(statusReports<=0) {
		equilibrate(duration);
		return;
	}
	double stepLength = duration / (double)statusReports; double eTime = 0;
	for(int i=0; i<statusReports; i++)
	{
		equilibrate(stepLength);
		eTime+=stepLength;
		cout<<"Equilibration has now elapsed for: "<<eTime<<" seconds."<<endl;
	}

}

void System::outputAllObservableNames()
{

	////////////////
	// NOTE!!!  IF YOU CHANGE ANYTHING HERE, BE SURE TO UPDATE BOTH THE GDAT FORMAT AND CSV FORMAT!!!

	if(!useBinaryOutput) {
		if(!csvFormat) {
			outputFileStream<<"#          time";
			//for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
			//	(*molTypeIter)->outputObservableNames(outputFileStream);

			int totalSpaces = 16;

			for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++) {
				string nm = (*obsIter)->getName();
				int spaces = totalSpaces-nm.length();
				if(spaces<1) { spaces = 1; }
				for(int k=0; k<spaces; k++) {
					outputFileStream<<" ";
				}
				outputFileStream<<nm;;
			}

			if(outputGlobalFunctionValues)
				for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
				{
					string nm = (*functionIter)->getNiceName();
					int spaces = totalSpaces-nm.length();
					if(spaces<1) { spaces = 1; }
					for(int k=0; k<spaces; k++) {
						outputFileStream<<" ";
					}
					outputFileStream<<nm;;
				}
			if(outputEventCounter) {
				string nm = "EventCount";
				int spaces = totalSpaces-nm.length();
				if(spaces<1) { spaces = 1; }
				for(int k=0; k<spaces; k++) {
					outputFileStream<<" ";
				}
				outputFileStream<<nm;;
			}

			outputFileStream<<endl;
		} else {

			// CSV FORMATTED OUTPUT
			outputFileStream<<"time";

			for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++) {
				string nm = (*obsIter)->getName();
				outputFileStream<<", "<<nm;;
			}

			if(outputGlobalFunctionValues)
				for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
				{
					string nm = (*functionIter)->getNiceName();
					outputFileStream<<", "<<nm;;
				}
			if(outputEventCounter) {
				string nm = "EventCount";
				outputFileStream<<", "<<nm;;
			}

			outputFileStream<<endl;
		}
	} else {
		cout<<"Warning: You cannot output observable names when outputting in Binary Mode."<<endl;
	}
}


void System::outputAllObservableCounts()
{
	outputAllObservableCounts(this->current_time,globalEventCounter);
}

void System::outputAllObservableCounts(double time)
{
	outputAllObservableCounts(time,globalEventCounter);
}



void System::outputAllObservableCounts(double cSampleTime, int eventCounter)
{
	if(!onTheFlyObservables)
	{
		for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++)
		{	(*obsIter)->clear();   }

		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		{	(*molTypeIter)->addAllToObservables(); 	}

		int match = 0;

	  	// NETGEN -- this bit replaces the commented block below
	  	Complex * complex;
	  	allComplexes.resetComplexIter();
	  	while(  (complex = allComplexes.nextComplex()) )
	  	{
	  		if( complex->isAlive() )
	  		{
	  			for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++)
	  			{
	  				match = (*obsIter)->isObservable( complex );
	  				for (int k=0; k<match; k++) (*obsIter)->straightAdd();
	  			}
	  		}
	  	}
		/*
		for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++) {
			if((*complexIter)->isAlive()) {
				for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++) {
					match = (*obsIter)->isObservable((*complexIter));
					for(int k=0; k<match; k++) (*obsIter)->straightAdd();
				}
			}
		}
		*/
	}


	if(useBinaryOutput) {
		double count=0.0; int oTot=0;

		outputFileStream.write((char *)&cSampleTime, sizeof(double));
		for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++) {
			count=((double)(*obsIter)->getCount());
			outputFileStream.write((char *) &count, sizeof(double));
		}
		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ ) {
				count=FuncFactory::Eval((*functionIter)->p);
				outputFileStream.write((char *) &count, sizeof(double));
			}

		if(outputEventCounter) {
			count=eventCounter;
			outputFileStream.write((char *) &count, sizeof(double));
		}
	}
	else {
		if(!csvFormat) {
			outputFileStream<<" "<<cSampleTime;
			for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++) {
				outputFileStream<<"  "<<((double)(*obsIter)->getCount());
			}

			if(outputGlobalFunctionValues)
				for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
					outputFileStream<<"  "<<FuncFactory::Eval((*functionIter)->p);
			if(outputEventCounter) {
				outputFileStream<<"  "<<eventCounter;
			}

			outputFileStream<<endl;
		} else {
			outputFileStream<<cSampleTime;
			for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++) {
				outputFileStream<<", "<<((double)(*obsIter)->getCount());
			}

			if(outputGlobalFunctionValues)
				for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
					outputFileStream<<", "<<FuncFactory::Eval((*functionIter)->p);
			if(outputEventCounter) {
				outputFileStream<<", "<<eventCounter;
			}

			outputFileStream<<endl;
		}
	}



}

void System::printAllObservableCounts()
{
	printAllObservableCounts(current_time,globalEventCounter);
}

void System::printAllObservableCounts(double cSampleTime)
{
	printAllObservableCounts(cSampleTime,globalEventCounter);
}

void System::printAllObservableCounts(double cSampleTime,int eventCounter)
{
	cout<<"Time";
	for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++)
		cout<<"\t"<<(*obsIter)->getName();
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			cout<<"\t"<<(*functionIter)->getNiceName();
	if(outputEventCounter) {
		cout<<"\tEventCount";
	}

	cout<<endl;

  	cout<<cSampleTime;
	for(obsIter = obsToOutput.begin(); obsIter != obsToOutput.end(); obsIter++)
		cout<<"\t"<<(*obsIter)->getCount();
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			cout<<"\t"<<FuncFactory::Eval((*functionIter)->p);
	if(outputEventCounter) {
		cout<<"\t"<<eventCounter;
	}
	cout<<endl;
}


// NETGEN  moved to ComplexList
/*
void System::printAllComplexes()
{
	cout<<"All System Complexes:"<<endl;
	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
		(*complexIter)->printDetails();
	cout<<endl;
}
*/


bool System::saveSpecies(string filename)
{
	bool debugOut = false;

	//open the output filestream
	ofstream speciesFile;
	speciesFile.open(filename.c_str());
	if(!speciesFile.is_open()) {
		cerr<<"Error in System when calling System::saveSpecies(string)!  Cannot open output stream to file "<<filename<<". "<<endl;
		cerr<<"quitting."<<endl;
		exit(1);
	}

	cout<<"\n\nsaving list of final molecular species..."<<endl;

	// create a couple data structures to store results as we go
	list <Molecule *> molecules;
	list <Molecule *>::iterator iter;
	map <int,bool> reportedMolecules;
    map <string,int> reportedSpecies;


	// loop over all the types of molecules that exist
	for( unsigned int k=0; k<allMoleculeTypes.size(); k++)
	{
		// retrieve the MoleculeType
		MoleculeType *mt = allMoleculeTypes.at(k);

		// loop over every individual molecule
		for(int j=0; j<mt->getMoleculeCount(); j++)
		{
			// check if we have looked this molecule before, and skip it if we have
			if(reportedMolecules.find(mt->getMolecule(j)->getUniqueID())!=reportedMolecules.end()) {
				continue;
			}



			//otherwise, we have not visited this particular species before, so loop over the molecules
			//that make up the species
			string speciesString = "";
			molecules.clear();
			mt->getMolecule(j)->traverseBondedNeighborhood(molecules,ReactionClass::NO_LIMIT);

			// key: partnerID1, bsiteID1, partnerID2, bsiteID2
			vector <vector <int> * > bondNumberMap;

			bool isFirst = true;
			for( iter = molecules.begin(); iter != molecules.end(); iter++ )
			{
				Molecule *m = (*iter);
				reportedMolecules.insert(pair <int,bool> (m->getUniqueID(),true));

				//Fist, output the molecule name
				if(isFirst) { speciesString += m->getMoleculeTypeName()+"("; isFirst=false; }
				else { speciesString += "."+m->getMoleculeTypeName()+"("; }

				//Go through each component of the molecule
				for(int s=0; s<m->getMoleculeType()->getNumOfComponents(); s++)
				{
					// output the component name
					string compName = m->getMoleculeType()->getComponentName(s);
					if(m->getMoleculeType()->isEquivalentComponent(s)) {
						// symmetric site, so we need to look up its sym name
						compName = m->getMoleculeType()->getEquivalenceClassComponentNameFromComponentIndex(s);
					}
					if(s==0) speciesString += compName;
					else speciesString += ","+compName;


					//output the state of the component, if it is set
					if(m->getComponentState(s)>=0) {
						speciesString += "~" + m->getMoleculeType()->getComponentStateName(s,m->getComponentState(s));
					}


					// check if the component is bound, if so we have to output a bond
					// we will label the bond incrementally, but we have to check to make
					// sure the bond wasn't declared earlier.  that's what the vector of int vectors is for.
					if(m->isBindingSiteBonded(s)) {
						if(debugOut) cout<<"binding site is bonded"<<endl;
						int partnerID = m->getBondedMolecule(s)->getUniqueID();
						int partnerSite = m->getBondedMoleculeBindingSiteIndex(s);
						int thisBondNumber = -1;

						// create the key
						vector <int> *key = new vector<int>(4);
						if(partnerID<m->getUniqueID()) {
							key->at(0) = partnerID; key->at(1) = partnerSite;
							key->at(2) = m->getUniqueID(); key->at(3)=s;
						} else {
							key->at(2) = partnerID; key->at(3) = partnerSite;
							key->at(0) = m->getUniqueID(); key->at(1)=s;
						}

						//search if that key was already inserted
						bool foundExistingBond = false;
						for(unsigned int bnmIndex =0; bnmIndex < bondNumberMap.size(); bnmIndex++) {
							if( key->at(0)==bondNumberMap.at(bnmIndex)->at(0) &&
								key->at(1)==bondNumberMap.at(bnmIndex)->at(1) &&
							    key->at(2)==bondNumberMap.at(bnmIndex)->at(2) &&
							    key->at(3)==bondNumberMap.at(bnmIndex)->at(3) ) {
								    thisBondNumber = bnmIndex+1;
									foundExistingBond = true;
									if(debugOut) cout<<"Found bond number: "<<thisBondNumber<<endl;
									break;
							}
						}

						//If it was not found, then insert it
						if(!foundExistingBond) {
							bondNumberMap.push_back(key);
							thisBondNumber = bondNumberMap.size();
							if(debugOut) cout<<"Creating bond number: "<<bondNumberMap.size()<<endl;
						}
						speciesString += "!" + NFutil::toString(thisBondNumber);
					}

				}


				speciesString += ")";
			}

			if(reportedSpecies.find(speciesString) != reportedSpecies.end()) {
				reportedSpecies[speciesString] = reportedSpecies[speciesString] + mt->getMolecule(j)->getPopulation();
			} else {
				reportedSpecies.insert(pair <string,int> (speciesString, mt->getMolecule(j)->getPopulation()));
			}

			//speciesString += "  1";
			//cout<<speciesString<<endl;

			if(debugOut) cout<<endl<<endl;

			//delete elements of the map
			while(bondNumberMap.size()>0) {
				vector <int> *v = bondNumberMap.at(bondNumberMap.size()-1);
				bondNumberMap.pop_back();
				delete v;
			}

			if(debugOut) cout<<endl<<endl;

		}
	}


	speciesFile<<"# nfsim generated species list for system: '"<< this->name <<"'\n";
	speciesFile<<"# warning! this feature is not yet fully tested! \n";
	for ( map<string,int>::iterator  it=reportedSpecies.begin() ; it != reportedSpecies.end(); it++ )
		speciesFile << (*it).first << "  " << (*it).second << "\n";
	speciesFile.flush();
	speciesFile.close();
	return true;
}


void System::printAllReactions()
{
	recompute_A_tot();
	cout<<"All System Reactions:"<<endl;
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
	{
		(*rxnIter)->printDetails();
	}
	cout<<endl;
}


void System::printAllMoleculeTypes()
{
	cout<<"All System Molecule Types:"<<endl;
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		(*molTypeIter)->printDetails();
	}
	cout<<endl;
}


// NETGEN  moved to ComplexList
/*
void System::outputComplexSizes(double cSampleTime)
{
	int size = 0;
	outputFileStream<<"\t"<<cSampleTime;
	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getComplexSize();
		if(size!=0) outputFileStream<<"\t"<<size;
	}
	outputFileStream<<endl;
}


double System::outputMeanCount(MoleculeType *m)
{
	int count = 0;
	int sum = 0;
	int allSum = 0;
	int allCount=0;
	int size=0;
	outputFileStream<<"\t"<<current_time;
	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);
		if(size>=2) { count++; sum+=size;}
		if(size>=1) { allSum+=size; allCount++; }

	}
	//cout<<sum<<"/"<<count<<"   "<<allSum<<"/"<<allCount<<endl;
	if(count!=0) {
		outputFileStream<<"\t"<<((double)sum/(double)count)<<endl;
		return ((double)sum/(double)count);
	}
	else
	{
		outputFileStream<<"\t"<<0.0<<endl;
		return 0.0;
	}

	return ((double)sum/(double)count);
}


double System::calculateMeanCount(MoleculeType *m)
{
	int count = 0;
	int sum = 0;
	int allSum = 0;
	int allCount=0;
	int size=0;

	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);
		if(size>=2) { count++; sum+=size; }
		if(size>=1) { allSum+=size; allCount++; }
	}
	return ((double)sum/(double)count);
}

void System::outputMoleculeTypeCountPerComplex(MoleculeType *m)
{
	int size = 0;
	outputFileStream<<"\t"<<current_time;
	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
	{
		size = (*complexIter)->getMoleculeCountOfType(m);

		if(size>=1) outputFileStream<<"\t"<<size;
	}
	outputFileStream<<endl;

}
*/

void System::printIndexAndNames()
{
	cout<<"All System Molecules:"<<endl;
	int idxCounter = 0;
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		cout<<idxCounter++<<"\t"<<(*molTypeIter)->getName()<<endl;
	}
	cout<<endl<<"All System Rxns:"<<endl;
	idxCounter = 0;
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
	{
		cout<<idxCounter++<<"\t"<<(*rxnIter)->getName()<<endl;
	}
	cout<<endl;
}



void System::addLocalFunction(LocalFunction *lf) {
	localFunctions.push_back(lf);
}


void System::evaluateAllLocalFunctions() {

	//Don't do all the work if we don't actually have to...
	if(localFunctions.size()==0) return;

	molList.clear();

	//loop through each moleculeType
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {

		//Loop through each molecule of that type
		for(int m=0; m<(*molTypeIter)->getMoleculeCount(); m++) {
			Molecule *mol = (*molTypeIter)->getMolecule(m);

			//Only continue if we haven't yet evaluated on this complex
			if(!mol->hasEvaluatedMolecule) {

				//First, grab the molecules in the complex
				//cout<<"in evaluate all local functions"<<endl;
				mol->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);

				//Evaluate all local functions on this complex
				for(unsigned int l=0; l<localFunctions.size(); l++) {
						//cout<<"--------------Evaluating local function on species..."<<endl;
						double val =localFunctions.at(l)->evaluateOn(mol,LocalFunction::SPECIES);
						//cout<<"     value of function: "<<val<<endl;

				}

				//Let those molecules know they've been visited
				for(molListIter=molList.begin(); molListIter!=molList.end(); molListIter++) {
					(*molListIter)->hasEvaluatedMolecule=true;
				}

				//clear the list
				molList.clear();
			}
		}
	}

	// Now go back and clear all the molecules of thier local functions...
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		for(int m=0; m<(*molTypeIter)->getMoleculeCount(); m++)
			(*molTypeIter)->getMolecule(m)->hasEvaluatedMolecule=false;
	}


}


GlobalFunction * System::getGlobalFunctionByName(string fName) {

	//First, look for the function directly in the list of global functions
	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
		if((*functionIter)->getName()==fName) {
			return (*functionIter);
		}

	//If it's not there, look up the global function reference that matches, then look up
	//the referenced function.
//	for(int i=0; i<(int)compositeFunctions.size(); i++) {
//
//	}
//
//	for( int i=0; i<(int)functionReferences.size(); i++) {
//		if(functionReferences.at(i)->name==fName) {
//			return getGlobalFunctionByName(functionReferences.at(i)->referencedFuncName);
//		}
//	}


	//cout<<"!!Warning, the system could not identify the global function: "<<fName<<".\n";
	//cout<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;
}

CompositeFunction * System::getCompositeFunctionByName(string fName)
{
	for( int i=0; i<(int)compositeFunctions.size(); i++) {
		if(compositeFunctions.at(i)->getName()==fName) {
			return compositeFunctions.at(i);
		}
	}
	//cout<<"!!Warning, the system could not identify the composite function: "<<fName<<".\n";
	//cout<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;
}

void System::finalizeCompositeFunctions()
{
	for( int i=0; i<(int)compositeFunctions.size(); i++) {
		compositeFunctions.at(i)->finalizeInitialization(this);
	}
}


LocalFunction * System::getLocalFunctionByName(string fName)
{
	for( int i=0; i<(int)localFunctions.size(); i++) {
		if(localFunctions.at(i)->getName()==fName) {
			return localFunctions.at(i);
		}
	}
	//cout<<"!!Warning, the system could not identify the local function: "<<fName<<".\n";
	//cout<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;

}


bool System::addCompositeFunction(CompositeFunction *cf) {
	this->compositeFunctions.push_back(cf);
	return true;
}




Observable * System::getObservableByName(string obsName)
{
	for(unsigned int i=0; i<obsToOutput.size(); i++) {
		if(obsToOutput.at(i)->getName().compare(obsName)==0) {
			return obsToOutput.at(i);
		}
	}

	cout.flush();
	cerr<<"!!Warning, the system could not identify the observable: "<<obsName<<".\n";
	cerr<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;
}



void System::addParameter(string name,double value) {
	this->paramMap[name]=value;
}
double System::getParameter(string name) {
	return this->paramMap.find(name)->second;
}
void System::setParameter(string name, double value) {
	if(paramMap.find(name)==paramMap.end()) {
		cout<<"Warning! System parameter: '"<<name<<"' does not exist and will not be updated."<<endl;
		return;
	}
	this->paramMap[name]=value;
}
void System::updateSystemWithNewParameters() {

	//Update all global functions
	for(unsigned int i=0; i<this->globalFunctions.size(); i++) {
		globalFunctions.at(i)->updateParameters(this);
	}

	//Update all local functions
	for(unsigned int i=0; i<this->localFunctions.size(); i++) {
		localFunctions.at(i)->updateParameters(this);
	}

	//Update all composite functions
	for(unsigned int i=0; i<this->compositeFunctions.size(); i++) {
		compositeFunctions.at(i)->updateParameters(this);
	}

	this->evaluateAllLocalFunctions();


	//Update all reactions
	for(unsigned int r=0; r<allReactions.size(); r++) {
		allReactions.at(r)->resetBaseRateFromSystemParamter();
	}


	//Update Atot (the total propensity of the system)
	this->recompute_A_tot();

}
void System::printAllParameters() {
	if(paramMap.size()==0) cout<<"no system parameters to print."<<endl;
	else cout<<"List of all system parameters:"<<endl;
	map<string,double>::iterator iter;
	for( iter = paramMap.begin(); iter != paramMap.end(); iter++ ) {
		cout << "\t" << iter->first << " = " << iter->second << endl;
	}
}

void System::printAllFunctions() {
	cout<<"System Global Functions: "<<endl;
	for(unsigned int i=0; i<this->globalFunctions.size(); i++) {
		globalFunctions.at(i)->printDetails(this);
	}

	cout<<"\nSystem Composite Functions: "<<endl;
	for(unsigned int i=0; i<this->compositeFunctions.size(); i++) {
		compositeFunctions.at(i)->printDetails(this);
	}

	cout<<"\nSystem Local Functions: "<<endl;
	for(unsigned int i=0; i<this->localFunctions.size(); i++) {
		localFunctions.at(i)->printDetails(this);
	}
}

void System::outputAllPropensities(double time, int rxnFired)
{
	if(!propensityDumpStream.is_open()) {

		string filename = this->name+"_propensity.txt";
		propensityDumpStream.open(filename.c_str());


		if(!outputFileStream.is_open()) {
				cerr<<"Error in System!  cannot open output stream to file "<<filename<<". "<<endl;
				cerr<<"quitting."<<endl;
				exit(1);
		}

		propensityDumpStream<<"time rxn";
		for(unsigned int r=0; r<allReactions.size(); r++) {
			propensityDumpStream<<" ";
			propensityDumpStream<<allReactions[r]->getName();
			for(int rl=0; rl<allReactions[r]->getNumOfReactants(); rl++) {
				propensityDumpStream<<" rL"<<NFutil::toString(rl);
			}
		}
		propensityDumpStream<<endl;
	}

	propensityDumpStream<<time<<" "<<allReactions.at(rxnFired)->getName();
	for(unsigned int r=0; r<allReactions.size(); r++) {
		propensityDumpStream<<" ";
			propensityDumpStream<<allReactions[r]->get_a();
			for(int rl=0; rl<allReactions[r]->getNumOfReactants(); rl++) {
				propensityDumpStream<<" "<<NFutil::toString((int)allReactions[r]->getReactantCount(rl));
		}
	}
	propensityDumpStream<<endl;


}






NFstream& System::getOutputFileStream()
{
    return outputFileStream;
}


// friend functions
template<class T>
NFstream& operator<<(NFstream& nfstream, const T& value)
{
    if (nfstream.useFile_)
	nfstream.file_ << value;
    else
	nfstream.str_ << value;

    return nfstream;
}

