
#include "NFcore.hh"

#include <math.h>
#include <fstream>

using namespace std;
using namespace NFcore;

int System::NULL_EVENT_COUNTER = 0;

/*!
   Creates a system that does not keep track of complexes.
 */
System::System(string name)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = false;
//	this->go = NULL;
	this->outputGlobalFunctionValues=false;
}


/*!
   Constructor that creates a System that has the option of keeping track of complexes
 */
System::System(string name, bool useComplex)
{
	this->name = name;
	this->a_tot = 0;
	current_time = 0;
	nextReaction = 0;
	this->useComplex = useComplex;
//	this->go = NULL;
	this->outputGlobalFunctionValues=false;
}


/*!
  Standard deconstructor for a system that cleans everything up.
 */
System::~System()
{	
	//Need to delete reactions
  	for(unsigned int r=0; r<allReactions.size(); r++)
  		delete [] rxnIndexMap[r];
  	delete [] rxnIndexMap;
	
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
	
	//And finally delete all the groups
//	Group *g;
//	while(allGroups.size()>0)
//	{
//		g = allGroups.back();
//		allGroups.pop_back();
//		delete g;
//	}

	nextReaction = 0;
	
//	if(go!=NULL)
//		delete go;
	
	//Close our connections to output files
	outputFileStream.close();
}


void System::registerOutputFileLocation(string filename)
{
	outputFileStream.open(filename.c_str());
	outputFileStream.setf(ios::scientific);
}

void System::changeOutputFileLocation(string newFilename) {
	outputFileStream.close();
	outputFileStream.open(newFilename.c_str());
	outputFileStream.setf(ios::scientific);
}
int System::addMoleculeType(MoleculeType *MoleculeType)
{
	allMoleculeTypes.push_back(MoleculeType);
	return (allMoleculeTypes.size()-1);
}


void System::addReaction(ReactionClass *reaction)
{
	reaction->init();
	allReactions.push_back(reaction);
}


//int System::addGroup(Group * g)
//{
//	allGroups.push_back(g);
//	return (allGroups.size()-1);
//}


int System::createComplex(Molecule * m)
{
	if(!useComplex) return -1;  //Only create complexes if we intend on using them...
	int c_id = allComplexes.size();
	Complex * c = new Complex(this, c_id, m);
	allComplexes.push_back(c);
	return c_id;
}

void System::addGlobalFunction(GlobalFunction *gf)
{
	this->globalFunctions.push_back(gf);
}



void System::updateGroupProperty(char * groupName, double *value, int n_values)
{
	cout<<"Updating group property for groups named: " << groupName << endl;
	cout<<"!! Not implemented.  I just did nothing! "<<endl;
	
}

MoleculeType * System::getMoleculeTypeByName(string& mName)
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


double System::getAverageGroupValue(string groupName, int valIndex)
{
	return 0;
//	double sum = 0;
//	int count = 0;
//	for(groupIter = allGroups.begin(); groupIter != allGroups.end(); groupIter++ )
//	{
//		string name = (*groupIter)->getName();
//		if(name==groupName)
//		{
//			sum += (*groupIter)->getValue(valIndex);
//			count ++;
//		}
//	}
//	return (sum/count);
}


void System::updateAllGroupProperty(double *value, int n_values)
{
	//cout<<"Updating group property for all groups, new value[0]: " << value[0] << endl;
	
//	for(groupIter = allGroups.begin(); groupIter != allGroups.end(); groupIter++ )
//	{
//		(*groupIter)->updateGroupProperty(value, n_values);
//	}
	
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
	
	rxnIndexMap = new int * [allReactions.size()];
  	for(unsigned int r=0; r<allReactions.size(); r++)
  	{
  		rxnIndexMap[r] = new int[allReactions.at(r)->getNumOfReactants()];
  		allReactions.at(r)->setRxnId(r);
  	}
  	
  	
	//This means we aren't going to add any more molecules to the system, so prep the rxns
	for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
		(*rxnIter)->prepareForSimulation();
	
  	//prep each molecule type for the simulation
  	for( molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
  		(*molTypeIter)->prepareForSimulation();
  	
  	if(BASIC_MESSAGE) cout<<"preparing the system...\n"; //printIndexAndNames();
  	
  	
//  if(go!=NULL)
// 	{
//  		go->writeGroupKeyFile();
// 		go->writeOutputFileHeader();
// 	}
  	
  	
  	
  	//prep each molecule type for the simulation
  	for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
  		(*functionIter)->prepareForSimulation(this);
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
	
	
	double delta_t = 0; unsigned long long iteration = 0, stepIteration = 0;
	double end_time = current_time+duration;
	while(current_time<end_time)
	{
		//2: Recompute a_tot for this time
		recompute_A_tot();
		if(DEBUG) cout<<" Determine a_tot : " << a_tot;
		
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
			cout<<"\t Reactions Cycles during this step: "<<stepIteration<<endl;
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
		//2: Recompute a_tot for this time
		recompute_A_tot();
		
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
  	outputFileStream<<"#\tTime";
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		(*molTypeIter)->outputObservableNames(outputFileStream);
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			outputFileStream<<"\t"<<(*functionIter)->getNiceName();
	outputFileStream<<endl;
}

void System::outputAllObservableCounts()
{
  	outputFileStream<<"\t"<<current_time;
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		(*molTypeIter)->outputObservableCounts(outputFileStream);
	
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			outputFileStream<<"\t"<<FuncFactory::Eval((*functionIter)->p);
	outputFileStream<<endl;
}

void System::outputAllObservableCounts(double cSampleTime)
{
  	outputFileStream<<"\t"<<cSampleTime;
	for(molTypeIter = allMoleculeTypes.begin(); molTypeIter != allMoleculeTypes.end(); molTypeIter++ )
		(*molTypeIter)->outputObservableCounts(outputFileStream);
	
	if(outputGlobalFunctionValues)
		for( functionIter = globalFunctions.begin(); functionIter != globalFunctions.end(); functionIter++ )
			outputFileStream<<"\t"<<FuncFactory::Eval((*functionIter)->p);
	outputFileStream<<endl;
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



//void System::addGroupOutputter(GroupOutputter * go)
//{
//	this->go = go;
//}


//void System::outputGroupDataHeader()
//{
//	if(this->go!=NULL)
//		go->writeOutputFileHeader();
//}


//void System::outputGroupData(double cSampleTime)
//{
//	if(this->go!=NULL)
//		go->writeStateToOutputFile(cSampleTime);
//}





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
