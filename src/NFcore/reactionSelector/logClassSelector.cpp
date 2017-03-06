/*
 * logClassSelector.cpp
 *
 *  Created on: Jul 23, 2009
 *      Author: msneddon
 */


#include "reactionSelector.hh"


using namespace std;
using namespace NFcore;


LogClassSelector::LogClassSelector(vector <ReactionClass *> &rxns) :
	ReactionSelector()
{
	//First, make sure the Rxn IDs match the index, so we don't get confused later
	for(unsigned int r=0; r<rxns.size(); r++) {
		rxns.at(r)->update_a();
		if((int)r!=rxns.at(r)->getRxnId()) {
			cerr<<"Internal Error in LogClassSelector: RxnIDs do not match position in vector."<<endl;
			cerr<<"For now, just turn off the LogClassSelector."<<endl;
			exit(1);
		}
	}

	//Now we initialize the basics
	maxClassLimit = 30;
	minClassLimit = 30;
	Atot = 0;
	n_reactions = rxns.size();

	n_activeLogClasses = 0;


	//Create the log class data structures
	if(maxClassLimit>minClassLimit) totalLogClassCount=maxClassLimit*2+1;
	else totalLogClassCount=minClassLimit*2+1;

	this->logClassList = new ReactionClass **[totalLogClassCount];
	this->logClassList = &logClassList[(totalLogClassCount-1)/2];

	this->logClassSize = new int [totalLogClassCount];
	this->logClassSize = &logClassSize[(totalLogClassCount-1)/2];
	this->logClassCapacity = new int [totalLogClassCount];
	this->logClassCapacity = &logClassCapacity[(totalLogClassCount-1)/2];
	this->logClassPropensity = new double [totalLogClassCount];
	this->logClassPropensity = &logClassPropensity[(totalLogClassCount-1)/2];

	this->isLogClassActive = new bool [totalLogClassCount];
	this->isLogClassActive = &isLogClassActive[(totalLogClassCount-1)/2];


	this->mapRxnIdToLogClass = new int [n_reactions];
	this->mapRxnIdToLogClassPosition = new int [n_reactions];

	activeLogClasses = new int [totalLogClassCount];
	for(int i=0; i<totalLogClassCount; i++) {
		activeLogClasses[i]=0;
	}

	cout<<"totalLength = "<<totalLogClassCount<<endl;
	int startingClassCapacity=(int)n_reactions/10;
	for(int i=-(totalLogClassCount-1)/2; i<=(totalLogClassCount-1)/2; i++) {
		logClassSize[i]=0;
		logClassCapacity[i]=startingClassCapacity;
		logClassPropensity[i]=0;
		isLogClassActive[i]=false;

		ReactionClass ** singleLogClass = new ReactionClass *[startingClassCapacity];
		for(int k=0; k<startingClassCapacity; k++) singleLogClass[k]=0;
		logClassList[i] = singleLogClass;
	}


	//Put the reactions in their starting classes
	int currentClass = 0; double current_a=0;
	for(int r=0; r<n_reactions; r++)
	{
		mapRxnIdToLogClass[r]=0;
		mapRxnIdToLogClassPosition[r]=-1;

		current_a = rxns.at(r)->get_a();
		currentClass=calculateClass(current_a);
		place(rxns.at(r),currentClass,current_a);
	}

	//Print debug message
	cout<<endl<<endl<<endl;
	for(int i=-(totalLogClassCount-1)/2; i<=(totalLogClassCount-1)/2; i++) {
		cout<<"logClassList["<<i<<"], size= "<<logClassSize[i];
		cout<<" / "<<logClassCapacity[i]<<"  atot= "<<logClassPropensity[i];
		cout<<" is active: "<<isLogClassActive[i]<<endl;
	}

	refactorPropensities();

	//Print debug message
	cout<<endl<<endl<<endl;
	for(int i=-(totalLogClassCount-1)/2; i<=(totalLogClassCount-1)/2; i++) {
		cout<<"logClassList["<<i<<"], size= "<<logClassSize[i];
		cout<<" / "<<logClassCapacity[i]<<"  atot= "<<logClassPropensity[i];
		cout<<" is active: "<<isLogClassActive[i]<<endl;
	}
}


void LogClassSelector::setLogClassToActive(int logClass)
{
	isLogClassActive[logClass] = true;
	this->n_activeLogClasses++;

	int cPos = 0;
	for(int c=maxClassLimit; c>=-minClassLimit;c--) {
		if(isLogClassActive[c])
		{
			activeLogClasses[cPos] = c;
			cPos++;
		}
	}
}


void LogClassSelector::setLogClassToInactive(int logClass)
{
	isLogClassActive[logClass] = false;
	this->n_activeLogClasses--;

	int cPos = 0;
	for(int c=maxClassLimit; c>=-minClassLimit;c--) {
		if(isLogClassActive[c])
		{
			activeLogClasses[cPos] = c;
			cPos++;
		}
	}
}


void LogClassSelector::place(ReactionClass *r,int logClass,double a)
{
	//If capacity is exceeded, then we have to allocate more space
	if(this->logClassSize[logClass]>=this->logClassCapacity[logClass]) {
		//Determine the new capacity for the log class
		int oldCap = logClassCapacity[logClass];
		int newCap = oldCap+oldCap/2;

		//create the new log class and copy over the new data
		ReactionClass ** singleLogClass = new ReactionClass *[newCap];
		for(int k=0; k<oldCap; k++)
			singleLogClass[k]=logClassList[logClass][k];
		for(int k=oldCap; k<newCap; k++)
			singleLogClass[k]=0;

		//delete the old data
		delete [] logClassList[logClass];

		//copy over the new array
		logClassList[logClass] = singleLogClass;
		logClassCapacity[logClass] = newCap;
	}

	//Put the new data in the correct log class, and remember where we
	//put this reaction in the logClass list
	logClassList[logClass][logClassSize[logClass]] = r;
	mapRxnIdToLogClass[r->getRxnId()]=logClass;
	mapRxnIdToLogClassPosition[r->getRxnId()]=logClassSize[logClass];

	//update the size and propensities
	if(!isLogClassActive[logClass]) setLogClassToActive(logClass);
	logClassSize[logClass]++;
	logClassPropensity[logClass]+=a;
}


LogClassSelector::~LogClassSelector()
{
	//Print debug message
	cout<<endl<<endl<<endl;
	for(int i=-(totalLogClassCount-1)/2; i<=(totalLogClassCount-1)/2; i++) {
		cout<<"logClassList["<<i<<"], size= "<<logClassSize[i];
		cout<<" / "<<logClassCapacity[i]<<"  atot= "<<logClassPropensity[i];
		cout<<" is active: "<<isLogClassActive[i]<<endl;
	}
}


double LogClassSelector::refactorPropensities()
{
	//First, we have to clear all
	Atot = 0;
	int trla_index=0;
	ReactionClass **tempRxnListArray = new ReactionClass *[n_reactions];

	for(int i=-(totalLogClassCount-1)/2; i<=(totalLogClassCount-1)/2; i++) {
		for(int k=0; k<logClassSize[i]; k++) {
			tempRxnListArray[trla_index] = logClassList[i][k];
			logClassList[i][k] = 0;
			tempRxnListArray[trla_index]->update_a();
			trla_index++;
		}
		logClassSize[i]=0;
		logClassPropensity[i]=0;
	}

	//Then we can reinsert the reactions into their classes
	int currentClass = 0; double current_a=0;
	for(int r=0; r<n_reactions; r++)
	{
		current_a = tempRxnListArray[r]->get_a();
		Atot += current_a;
		currentClass=calculateClass(current_a);
		place(tempRxnListArray[r],currentClass,current_a);
	}

	return Atot;
}


double LogClassSelector::update(ReactionClass *r,double oldA, double newA)
{
	int oldClass = mapRxnIdToLogClass[r->getRxnId()];
	int newClass = calculateClass(newA);

	//If the class doesn't change, just update the propensities
	if(oldClass==newClass) {
		logClassPropensity[newClass]-=oldA;
		logClassPropensity[newClass]+=newA;
	}

	//If the class does change, we need to move the reaction
	else {

		//First remove from the original log class, by swapping in the
		//rule in the last position, if there was one.  Otherwise, we just
		//remove it.
		int oldPos = mapRxnIdToLogClassPosition[r->getRxnId()];
		if(oldPos==(logClassSize[oldClass]-1)){
			//at the last pos already, so special handling is required
			logClassList[oldClass][logClassSize[oldClass]-1] = 0;
			logClassSize[oldClass]--;

		} else {
			logClassList[oldClass][oldPos] = logClassList[oldClass][logClassSize[oldClass]-1];
			logClassList[oldClass][logClassSize[oldClass]-1] = 0;
			mapRxnIdToLogClassPosition[logClassList[oldClass][oldPos]->getRxnId()]=oldPos;
			logClassSize[oldClass]--;
		}
		if(logClassSize[oldClass]==0) {
			logClassPropensity[oldClass] = 0;
			setLogClassToInactive(oldClass);
		}
		else logClassPropensity[oldClass]-=oldA;


		//And put the rxn into the new class
		this->place(r,newClass,newA);
	}

	Atot-=oldA;
	Atot+=newA;
	return Atot;
}


double LogClassSelector::getNextReactionClass(ReactionClass *&rc)
{
	//Generate the random number based on the total propensity
	double randNum = NFutil::RANDOM(Atot);

	//First, we select the next class to fire based on the propensities
	double a_sum=0; int selectedClass=0; int c=0;

	for(int actIndex=0; actIndex<n_activeLogClasses; actIndex++)
	{
		c=activeLogClasses[actIndex];
		//if(logClassSize[c]<=0) continue;
		a_sum += logClassPropensity[c];
		if(randNum <= a_sum)
		{
			selectedClass = c;
			break;
		}
	}

	//Then, we use a rejection method to select the next rule
	 int randRule=0; double randRule_A=0; double weight;
	 do {
		 randRule = NFutil::RANDOM_INT(0,logClassSize[selectedClass]);
		 randRule_A = logClassList[selectedClass][randRule]->get_a();
		 weight = pow(2,(float)(selectedClass+1))*NFutil::RANDOM(1);
	 } while (weight > randRule_A);

	 //we have our rule
	 rc=logClassList[selectedClass][randRule];

	 return -1;
}


double LogClassSelector::getAtot()
{
	return Atot;
}



int LogClassSelector::calculateClass(double a) {

	int logClass = 0; int i=0;
	if(a==0) {
		logClass = 0;
	}

	//Calculate the class for propensities larger than 1, by
	//bit shifting to divide by 2 over and over
	else if (a>=1) {
		int e = (int)a;
		while (e>1) {
			e = e >> 1;
			i+=1;
		}
		logClass = i;
	}

	//Now calculate them for propensities less than 1 by multiplying by 2
	else if (a<1) {
		double d = a;
		 while (d < 1) {
			 d *= 2;
			 i-= 1;
		 }
		logClass = i;
	}

	//Finally, make sure we are within the bounds of the specified limits
	if(logClass>maxClassLimit) logClass = maxClassLimit;
	else if(logClass<-minClassLimit) logClass = minClassLimit;

	return logClass;
}
