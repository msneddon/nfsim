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

//	Atot = 0;
//	n_reactions = 0;
//	delete [] reactionClassList;
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
			//cout<<"selected class propensity: "<<logClassPropensity[c]<<" and the a_sum: "<<a_sum<<endl;
			//cout<<" and the randNum: "<<randNum<<endl;
			break;
		}
	}


	//cout<<"Selected class: "<<selectedClass<<endl;

	//Then, we use a rejection method to select the next rule
	 int randRule=0; double randRule_A=0; double weight;
	 do {
		 randRule = NFutil::RANDOM_INT(0,logClassSize[selectedClass]);
		 //cout<<" ** "<<logClassSize[selectedClass]<<"  and the selected rule: "<<randRule<<endl;
		 randRule_A = logClassList[selectedClass][randRule]->get_a();
		 weight = pow(2,(float)(selectedClass+1))*NFutil::RANDOM(1);
	 } while (weight > randRule_A);


	 //we have our rule
	 rc=logClassList[selectedClass][randRule];

	// rc->printDetails();
	 return -1;










	//WARNING - DO NOT USE THE DEFAULT C++ RANDOM NUMBER GENERATOR FOR THIS STEP
	// - IT INTRODUCES SMALL NUMERICAL ERRORS CAUSING THE ORDER OF RXNS TO
	//   AFFECT SIMULATION RESULTS
//	for(int c=0; r<n_reactions; r++) {
//		a_sum += reactionClassList[r]->get_a();
//		if(randNum <= a_sum)
//		{
//			rc = reactionClassList[r];
//			return (randNum-last_a_sum);
//		}
//		last_a_sum = a_sum;
//	}
//
//	cerr<<"Error in Direct Reaction Selector: randNum exceeds a_sum!!!"<<endl;
//	cerr<<"randNum: "<<randNum<<"  a_sum: "<< a_sum<<" running a_tot:"<<Atot<<endl;
//	return -1;
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


	//       double totalPropensity = sv->getTotalPropensity();
	//       int lclass = 0;
	//       if (totalPropensity == 0) {
	//               lclass = numeric_limits<int>::min();
	//       } else if (totalPropensity > 1) {
	//               int e = (int)totalPropensity;
	//               while (e>1) {
	//                       e = e >> 1;
	//                       i+=1;
	//               }
	//               lclass = i;
	//       } else if (totalPropensity < 1) {
	//               double d = totalPropensity;
	//               while (d < 1) {
	//                       d *= 2;
	//                       i-= 1;
	//               }
	//               lclass = i;
	//       } else if (totalPropensity == 1)
	//               lclass = 0;
	//
	//       return lclass;

}






//
//void System::singleSpatialStep() {
//       double v = 0;
//       double weight = 0;
//       int chosen_class;
//       int chosen;
//       int size;
//       map<int, double>::iterator iter, end;
//       // this is the new simulation step
//       current_time -= (1./a_tot)*(1.-NR::ran2(idum));
//
//       cerr << "here 0: " << a_tot << endl;
//
//       // choose a logarithmic class
//       v = a_tot*NR::ran2(idum);
//       for (iter = lweimap.begin(); iter != lweimap.end(); iter++)
//               cerr << iter->first << " " << iter->second << endl;
//       iter = lweimap.begin();
//       weight = iter->second;
//       while (weight <= v) {
//               iter++;
//               weight += iter->second;
////              if (iter == lweimap.end()) {
//       //              cerr << "fuck something is wrong" << endl;
//       //              for (iter = lweimap.begin(); iter != lweimap.end(); iter++)
//                               //cerr << iter->first << " " << iter->second << endl;
//       //              exit(1);
//               //}
//
//       }
//       chosen_class = iter->first;
//       // choose a subvolume
//       cerr << "here1 chose class " << iter->first <<  endl;
//       size = numberOfSubvolumesInLClass[chosen_class];
//       cerr << "here2 " << size << endl;
//       Subvolume *sv;
//       do {
//               chosen = (int)(size*NR::ran2(idum));
//               sv = lsubmap[chosen_class][chosen];
//               weight = pow(2,chosen_class+1)*NR::ran2(idum);
//       } while (sv->getTotalPropensity() < weight);
//
//       cerr << "here3 chose " << sv->getIndex() << endl;
//
//       // remember the propensity of the chosen subvolume
//       double oldProp = sv->getTotalPropensity();
//
//       // vector for passing to pickReaction to obtain relevant information
//       // vec[0] tells me how many compartments are affected by firing a reaction
//       // vec[1] returns the index of the first affected subvolume
//       // vec[2] returns the index of the second affected subvolume
//       // vec[3] returns the propensity of the second affected subvolume before the
//reaction fired
//       vector<double> vec(4);
//       sv->pickReaction(vec);
//       cerr << "hier4" << endl;
//       // update the membership of the subvolume(s) in the logarithmic classes
//       updateSubvolume(sv, chosen_class, oldProp);
//       // was the reaction a diffusion reaction ?
//       if (vec[0] == 2) {
////              cerr << "updating second" << endl;
//               updateSubvolume((*subvolumes)[vec[2]], subvolumesLClasses[vec[2]], vec[3]);
//       } else if (vec[0] != 1 || vec[0] != 2) {
//               cerr << "something went seriously wrong in the singleSpatialStep of the
//logarithmic classes algorithm ... leaving ... "
//                        << endl;
//               exit(1);
//       }
//}
//
//void System::updateSubvolume(Subvolume *sv, int oldClass, double oldProp) {
//       cerr << "in updateSubvolume old " << oldClass <<  " " << oldProp;
//       // calculate new logarithmic class
//
//       sv->updatePropensities();
//       int newClass = calcLClass(sv);
//       cerr << " new class and prop: " <<  newClass << " " << sv->getTotalPropensity()
//<< " " << sv->getIndex() << endl;
//
//       // change propensities in logarithmic class
//       // note that the total propensity of the system is updated by
//updateRxnMembership of the molecule class
//       // as well as removeFromSubvolume if the molecule class. hence we do not need
//to do that here
//       lweimap[oldClass] -= oldProp;
//       lweimap[newClass] += sv->getTotalPropensity();
//
//       // only do something if subvolume is in a new subvolume
//       if (newClass != oldClass) {
//               // get the number of elements that can be stored in the old lclass
//               int oldMax = lsubmap[oldClass].size();
//               // get the number of element currently stored in the old lclass
//               int oldEnd = numberOfSubvolumesInLClass[oldClass];
//               // remember position of subvolume in vector ...
//               int oldPos = subvolumesLClassesPositions[sv->getIndex()];
//
//               // get the number of elements that can be stored in the new lclass
//               int newMax = lsubmap[newClass].size();
//               // get the number of elements currently stored in the new lclass
//               int newEnd = numberOfSubvolumesInLClass[newClass];
//
//               // change logarithmic class of this subvolume in subvolumesLClasses vector
//               subvolumesLClasses[sv->getIndex()] = newClass;
//
//               // remove subvolume from current class
//               if (oldPos == oldEnd-1)
//                       numberOfSubvolumesInLClass[oldClass] -= 1;
//               else {
//                       // index of subvolume that has to be moved to the now vacant position
//                       int movedSubvolumeIndex = lsubmap[oldClass].at(oldEnd-1)->getIndex();
//                       lsubmap[oldClass].at(oldPos) = (*subvolumes)[movedSubvolumeIndex];
//                       subvolumesLClassesPositions[movedSubvolumeIndex] = oldPos;
//                       numberOfSubvolumesInLClass[oldClass] -= 1;
//               }
//
//               // insert subvolume into new class
//               if (newEnd == newMax)
//                       lsubmap[newClass].push_back(sv);
//               else
//                       lsubmap[newClass].at(newEnd) = sv;
//               subvolumesLClassesPositions[sv->getIndex()] = newEnd;
//               numberOfSubvolumesInLClass[newClass]+=1;
//       }
//       lweimap[numeric_limits<int>::min()] = 0;
//}
//
//void System::prepareForLogarithmicClasses() {
//       // initialise the logarithmic classes of the system
//       int lclass = 0;
//       for (svit = subvolumes->begin(); svit != subvolumes->end(); ++svit) {
//               lclass = calcLClass(*svit);
//               lsubmap[lclass].push_back(*svit);
//               lweimap[lclass] += (*svit)->getTotalPropensity();
//               subvolumesLClasses.push_back(lclass);
//               subvolumesLClassesPositions.push_back(numberOfSubvolumesInLClass[lclass]++);
//       }
//       // calculate the total propensity
//       a_tot = 0;
//       map<int, double>::iterator mit;
//       for (mit = lweimap.begin(); mit != lweimap.end(); ++mit) {
//               a_tot += mit->second;
//       }
////      cerr << "a tot" << a_tot << endl;
//       lweimap[numeric_limits<int>::min()] = 0;
//}
//
//// sollte das hier nicht ins system ... ein subvolumen sollte nicht wissen zur
//welchen logarithmischen klasse es gehoert.
//// in der ueberarbeitung musst du den scheiss hier rausschmeissen.
//int System::calcLClass(Subvolume *sv) {
//       int i=0;
//       double totalPropensity = sv->getTotalPropensity();
//       int lclass = 0;
//       if (totalPropensity == 0) {
//               lclass = numeric_limits<int>::min();
//       } else if (totalPropensity > 1) {
//               int e = (int)totalPropensity;
//               while (e>1) {
//                       e = e >> 1;
//                       i+=1;
//               }
//               lclass = i;
//       } else if (totalPropensity < 1) {
//               double d = totalPropensity;
//               while (d < 1) {
//                       d *= 2;
//                       i-= 1;
//               }
//               lclass = i;
//       } else if (totalPropensity == 1)
//               lclass = 0;
//
//       return lclass;
//}


















