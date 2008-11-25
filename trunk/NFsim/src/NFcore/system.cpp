
#include "NFcore.hh"

#include <math.h>
#include <fstream>

using namespace std;
using namespace NFcore;

int System::NULL_EVENT_COUNTER = 0;



System::System(string name)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = false;
//	this->go = NULL;
	this->outputGlobalFunctionValues=false;
	rxnIndexMap=0;
	useBinaryOutput=false;
	onTheFlyObservables=true;
	universalTraversalLimit=-1;
	ds=0;
}


System::System(string name, bool useComplex)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = useComplex;
	this->outputGlobalFunctionValues=false;

	rxnIndexMap=0;
	useBinaryOutput=false;
	onTheFlyObservables=true;
	universalTraversalLimit=-1;
	ds=0;
}


System::~System()
{
	if(ds!=0) delete ds;

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

	//Delete all MoleculeTypes (which deletes all molecules and templates)
	MoleculeType *s;
	while(allMoleculeTypes.size()>0)
	{
		s = allMoleculeTypes.back();
		allMoleculeTypes.pop_back();
		delete s;
	}

	//Delete all the complexes
	Complex *c;
	while(allComplexes.size()>0)
	{
		c = allComplexes.back();
		allComplexes.pop_back();
		delete c;
	}

	GlobalFunction *gf;
	while(this->globalFunctions.size()>0)
	{
		gf = globalFunctions.back();
		globalFunctions.pop_back();
		delete gf;
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

void System::registerOutputFileLocation(string filename)
{
	if(outputFileStream.is_open()) { outputFileStream.close(); }
	if(useBinaryOutput) {
		outputFileStream.open((filename).c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
		//ios_base::out -- Set for output only
		//ios_base::binary --  Set output to binary
		//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there

		//Also, output a header file to keep track of the number
		ofstream headerFile;
		int tabCount=0;
		headerFile.open((filename+".head").c_str());
		headerFile<<"#\tTime"; tabCount++;
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
			int oTot = (*molTypeIter)->getNumOfObservables();
			for(int o=0; o<oTot; o++) {
				headerFile<<"\t"<<(*molTypeIter)->getObservable(o)->getAliasName();
				tabCount++;
			}
		}
		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ ) {
				headerFile<<"\t"<<(*functionIter)->getNiceName();
				tabCount++;
			}
		headerFile<<endl;
		for(int t=0; t<tabCount; t++) headerFile<<"\t";
		headerFile.close();

	} else {
		outputFileStream.open(filename.c_str());
		outputFileStream.setf(ios::scientific);
	}

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



int System::createComplex(Molecule * m)
{
	if(!useComplex) return -1;  //Only create complexes if we intend on using them...
	int c_id = allComplexes.size();
	Complex * c = new Complex(this, c_id, m);
	allComplexes.push_back(c);
	return c_id;
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



int System::getObservableCount(int moleculeTypeIndex, int observableIndex) const
{
	return allMoleculeTypes.at(moleculeTypeIndex)->getObservableCount(observableIndex);
}

Complex * System::getNextAvailableComplex()
{
	Complex * c = allComplexes.at(nextAvailableComplex.front());
	nextAvailableComplex.pop();
	return c;
}

void System::notifyThatComplexIsAvailable(int ID_complex)
{
	nextAvailableComplex.push(ID_complex);
}

void System::purgeAndPrintAvailableComplexList()
{
	cout<<"AvailableComplexes:";
	while(	!nextAvailableComplex.empty() )
	{
		cout<<" -> "<<nextAvailableComplex.front();
		nextAvailableComplex.pop();
	}
	cout<<endl;
}


//When you are ready to run the simulation (meaning that all moleculeTypes
//all molecules, and all reactions have been created and registered with
//the system) call this function to populate all the reactant lists and
//observables.
void System::prepareForSimulation()
{
	//Note!!  : the order of preparing the system matters!  You have to prepare
	//some things before others, because certain things require other

  	//First, we have to prep all the functions...
  	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
  		(*functionIter)->prepareForSimulation(this);



  	// now we prepare all reactions
	rxnIndexMap = new int * [allReactions.size()];
  	for(unsigned int r=0; r<allReactions.size(); r++)
  	{
  		rxnIndexMap[r] = new int[allReactions.at(r)->getNumOfReactants()];
  		allReactions.at(r)->setRxnId(r);
  	}


	//This means we aren't going to add any more molecules to the system, so prep the rxns
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
		(*rxnIter)->prepareForSimulation();

	//If there are local functions to be had, make sure we set up those local function lists in the molecules
	//before we try to add molecules to reactant lists
	if(this->localFunctions.size()>0) {
	  	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	  		(*molTypeIter)->setUpLocalFunctionListForMolecules();
	}

	this->evaluateAllLocalFunctions();

  	//prep each molecule type for the simulation
  	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
  		(*molTypeIter)->prepareForSimulation();



  	//if(BASIC_MESSAGE) cout<<"preparing the system...\n";
  	//printIndexAndNames();


//  if(go!=NULL)
// 	{
//  		go->writeGroupKeyFile();
// 		go->writeOutputFileHeader();
// 	}


  	recompute_A_tot();
}


double System::recompute_A_tot()
{
	//Loop through the reactions and add up the rates
	a_tot = 0;
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
	{
		a_tot += (*rxnIter)->update_a();
		if(DEBUG) (*rxnIter)->printDetails();
	}
	return a_tot;
}






/* select the next reaction, given a_tot has been calculated */
double System::getNextRxn()
{
	double randNum = NFutil::RANDOM(a_tot);

	double a_sum=0, last_a_sum=0;
	nextReaction = 0;

	//WARNING - DO NOT USE THE DEFAULT C++ RANDOM NUMBER GENERATOR FOR THIS STEP
	// - IT INTRODUCES SMALL NUMERICAL ERRORS CAUSING THE ORDER OF RXNS TO
	//   AFFECT SIMULATION RESULTS
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++)
	{
		a_sum += (*rxnIter)->get_a();
		if (randNum <= a_sum && nextReaction==0)
		{
			nextReaction = (* rxnIter);
			//cout<<"rNum: "<<randNum<<" last_a: "<<last_a_sum<<" a_sum "<<a_sum<<endl;
			return (randNum-last_a_sum);
		}
		last_a_sum = a_sum;
	}
	cerr<<"Error: randNum exceeds a_sum!!!"<<endl;
	cerr<<"randNum: "<<randNum<<"  a_sum: "<< a_sum<<endl;
	return -1;
}

/* main simulation loop */
double System::sim(double duration, long int sampleTimes)
{
	System::NULL_EVENT_COUNTER=0;
	cout.setf(ios::scientific);
	cout<<"Simulating system for: "<<duration<<" second(s)."<<endl<<endl;

	//First, output the header for the output of this simulation
	//outputAllObservableNames();

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
		//2: Recompute a_tot for this time
		//cout<<" a_tot was : " << a_tot<<endl;
		//recompute_A_tot();
		//cout<<" a_tot (after recomputing) is : " << a_tot<<endl;

		//3: Select next reaction time (making sure we have something that can react)
		//   dt = -ln(rand) / a_tot;
		//Choose a random number on the closed interval (0,1) so that we never
		//have a dt=0 or a dt=infinity
		if(a_tot>0) delta_t = -log(NFutil::RANDOM_CLOSED()) / a_tot;
		else { delta_t=0; current_time=end_time; }
		if(DEBUG) cout<<"   Determine dt : " << delta_t << endl;


		//Report everything up until the next step if we have to
		if(DEBUG) cout<<"  Current Sample Time: "<<curSampleTime<<endl;
		if((current_time+delta_t)>=curSampleTime)
		{
			while((current_time+delta_t)>=(curSampleTime))
			{
				if(curSampleTime>end_time) break;
				outputAllObservableCounts(curSampleTime);
//				outputGroupData(curSampleTime);
				curSampleTime+=dSampleTime;
			}
			//printAllReactions();
			cout<<"Sim time: "<<current_time<<"\tCPU time: ";
			cout<<(double(clock())-double(start))/CLOCKS_PER_SEC<<"s";
			cout<<"\t events: "<<stepIteration<<endl;
			//cout<<"\tAtot:"<<a_tot<<endl;
			stepIteration=0;
		}

		//Make sure we can react...
		if(delta_t==0) break;

		//4: Select next reaction class based on smallest j,
		//   such that sum of a_j over all j >= r2*a_tot
		double randElement = getNextRxn();
		//cout<<"Fire: "<<nextReaction->getName()<<" at time "<< current_time<<endl;
		//Output selected reaction for debugging
		//cout<<"\nFiring: "<< endl;
		//nextReaction->printDetails();
		//cout<<endl<<endl;


		//Increment time
		iteration++;
		stepIteration++;
		current_time+=delta_t;

		//5: Fire Reaction! (takes care of updates to lists and observables)
		nextReaction->fire(randElement);

		tryToDump();
	}


	finish = clock();
    time = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if(BASIC_MESSAGE)
    {
    	cout<<endl<<"You just simulated "<< iteration <<" reactions in "<< time << "s\n( ";
    	cout<<((double)iteration)/time<<" reactions/sec, ";
    	cout<<(time/((double)iteration))<<" CPU seconds/event )"<< endl;
    	cout<<"Null events: "<< System::NULL_EVENT_COUNTER;
    	cout<<"   ("<<(time)/((double)iteration-(double)System::NULL_EVENT_COUNTER)<<" CPU seconds/non-null event )"<< endl;
    }

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
		if(a_tot>0) delta_t = -log(NFutil::RANDOM_CLOSED()) / a_tot;
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

		//cout<<"Fire: "<<nextReaction->getName()<<" at time "<< current_time<<endl;

		//5: Fire Reaction! (takes care of updates to lists and observables)
		nextReaction->fire(randElement);
	}
	//cout<<"a_tot="<<a_tot;
	return current_time;
}

void System::singleStep()
{
	cout<<"  -Starting at time: "<<this->current_time<<endl;
	double delta_t = 0;

	recompute_A_tot();
	cout<<"  -total propensity (a_total) calculated as: "<<a_tot<<endl;
	if(a_tot>0) delta_t = -log(NFutil::RANDOM_CLOSED()) / a_tot;
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

	cout<<"  -Firing: "<<nextReaction->getName()<<endl;

	//5: Fire Reaction! (takes care of updates to lists and observables)
	nextReaction->fire(randElement);
	cout<<"  -System time is now"<<current_time<<endl;
}

void System::equilibriate(double duration)
{
	current_time = 0;
	stepTo(duration);
	current_time = 0;
}

void System::equilibriate(double duration, int statusReports)
{
	double stepLength = duration / (double)statusReports;
	double eTime = 0;
	for(int i=0; i<statusReports; i++)
	{
		equilibriate(stepLength);
		eTime+=stepLength;
		cout<<"Equilibriation has now elapsed for: "<<eTime<<" seconds."<<endl;

	}

}

void System::outputAllObservableNames()
{
	if(!useBinaryOutput) {
		outputFileStream<<"#\tTime";
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
			(*molTypeIter)->outputObservableNames(outputFileStream);

		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
				outputFileStream<<"\t"<<(*functionIter)->getNiceName();
		outputFileStream<<endl;
	} else {
		cout<<"Warning: You cannot output observable names when outputting in Binary Mode."<<endl;
	}
}

void System::outputAllObservableCounts()
{
	outputAllObservableCounts(this->current_time);
}



void System::outputAllObservableCounts(double cSampleTime)
{
	if(!onTheFlyObservables) {
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
			(*molTypeIter)->addAllToObservables();
		}
	}


	if(useBinaryOutput) {
		double count=0.0; int oTot=0;

		outputFileStream.write((char *)&cSampleTime, sizeof(double));
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {
			oTot = (*molTypeIter)->getNumOfObservables();
			for(int o=0; o<oTot; o++) {
				count=(double)((*molTypeIter)->getObservable(o)->getCount());
				outputFileStream.write((char *) &count, sizeof(double));
			}
		}
		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ ) {
				count=FuncFactory::Eval((*functionIter)->p);
				outputFileStream.write((char *) &count, sizeof(double));
			}
	}
	else {

		outputFileStream<<"\t"<<cSampleTime;
		for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
			(*molTypeIter)->outputObservableCounts(outputFileStream);

		if(outputGlobalFunctionValues)
			for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
				outputFileStream<<"\t"<<FuncFactory::Eval((*functionIter)->p);
		outputFileStream<<endl;
	}



}

void System::printAllObservableCounts(double cSampleTime)
{
	cout<<"Time";
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		(*molTypeIter)->printObservableNames();
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			cout<<"\t"<<(*functionIter)->getNiceName();
	cout<<endl;

  	cout<<cSampleTime;
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		(*molTypeIter)->printObservableCounts();
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			cout<<"\t"<<FuncFactory::Eval((*functionIter)->p);
	cout<<endl;
}

void System::printAllComplexes()
{
	cout<<"All System Complexes:"<<endl;
	for(complexIter = allComplexes.begin(); complexIter != allComplexes.end(); complexIter++ )
		(*complexIter)->printDetails();
	cout<<endl;
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

void System::printAllGroups()
{
	cout<<"can't print all groups: groups not available in this build"<<endl;
//	cout<<"All System Groups:"<<endl;
//	for(groupIter = allGroups.begin(); groupIter != allGroups.end(); groupIter++ )
//	{
//		(*groupIter)->printDetails();
//	}
//	cout<<endl;
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

	cout<<"Evaluating all local functions here in System..."<<endl;

	list <Molecule *> molList;
	list <Molecule *>::iterator molListIter;

	//loop through each moleculeType
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ ) {

		//Loop through each molecule of that type
		for(int m=0; m<(*molTypeIter)->getMoleculeCount(); m++) {
			Molecule *mol = (*molTypeIter)->getMolecule(m);

			//evaluate all functions on this Molecule that are local to a single molecule
			for(unsigned int l=0; l<localFunctions.size(); l++) {
				if(localFunctions.at(l)->getEvaluationLevel()>0) {
					cout<<"--------------Evaluating local function on single molecule..."<<endl;
					double val = localFunctions.at(l)->evaluateOn(mol);
					cout<<"     value of function: "<<val<<endl;
				}
			}


			//Only continue if we haven't yet evaluated on this complex
			if(!mol->hasEvaluatedMolecule) {

				//First, grab the molecules in the complex
				mol->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);

				//Evaluate all local functions on this complex
				for(unsigned int l=0; l<localFunctions.size(); l++) {
					if(localFunctions.at(l)->getEvaluationLevel()==0) {
						cout<<"--------------Evaluating local function on species..."<<endl;
						double val = localFunctions.at(l)->evaluateOn(mol);
						cout<<"     value of function: "<<val<<endl;
					}

				}

				//Let those molecules know they've been visited
				for(molListIter=molList.begin(); molListIter!=molList.end(); molListIter++) {
					(*molListIter)->hasEvaluatedMolecule=true;
				}
			}
		}
	}

	// Now go back and clear all the molecules of thier local functions...
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
	{
		for(int m=0; m<(*molTypeIter)->getMoleculeCount(); m++)
			(*molTypeIter)->getMolecule(m)->hasEvaluatedMolecule=false;
	}



	//Now, since we changed things around, we have to update the molecule positions in the
	//reactant trees of DOR reactions.
}


GlobalFunction * System::getGlobalFunctionByName(string fName) {
	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
		if((*functionIter)->getName()==fName) {
			return (*functionIter);
		}
	cout<<"!!Warning, the system could not identify the global function: "<<fName<<".\n";
	cout<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;
}

Observable * System::getObservableByName(string obsName)
{
	for(unsigned mt=0; mt<allMoleculeTypes.size(); mt++) {
		for(int ob=0; ob<allMoleculeTypes.at(mt)->getNumOfObservables(); ob++ ) {
			string name = allMoleculeTypes.at(mt)->getObservableAlias(ob);
			if(obsName==name) {
				return allMoleculeTypes.at(mt)->getObservable(ob);
			}
		}
	}
	cout<<"!!Warning, the system could not identify the observable: "<<obsName<<".\n";
	cout<<"The calling function might catch this, or your program might crash now."<<endl;
	return 0;
}
