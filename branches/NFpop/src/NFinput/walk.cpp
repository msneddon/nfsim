#include "NFinput.hh"

#include <iostream>



using namespace std;
using namespace NFcore;


//function declarations
int enterMainMenuLoop(System *s);
void enterStepLoop(System *s);

void getPrintout(System *s) ;
int getInput(int min, int max);
double getInput(double min);


void NFinput::walk(System *s)
{
	cout<<endl<<endl;
	cout<<"Let's take a walk with your system, shall we?"<<endl;
	cout<<"Please tell me what you would like to do."<<endl;
	while(true)
	{
		int exitCode = enterMainMenuLoop(s);
		if(exitCode==0) break;
	}

	cout<<endl<<"hope you had a good walk."<<endl;
}




void printMainMenu()
{
	cout<<"---------------------------"<<endl;
	cout<<" (0) equilibriate"<<endl;
	cout<<" (1) simulate"<<endl;
	cout<<" (2) enter stepper"<<endl;
	cout<<" (3) change output file"<<endl;
	cout<<" (4) print stuff"<<endl;
	cout<<" (5) exit"<<endl;
}

int enterMainMenuLoop(System *s)
{
	bool exit = false;
	int exitCode = 0;
	while(true)
	{
		cout<<endl;
		printMainMenu();
		int selection = getInput(0,5);
		double time = 0; int steps=0;
		switch(selection) {
			case 0:
			    cout<<"\nEnter equilibriation time: "<<endl;
			    time = getInput(0.0);
			    cout<<"\nequilibriating for "<< time <<" seconds.  Please wait..."<<endl;
			    s->equilibrate(time);
			    break;
			case 1:
				cout<<"\nEnter simulation time: "<<endl;
				time = getInput(0.0);
				cout<<"\nEnter number of output steps: "<<endl;
				steps = getInput(0,INT_MAX);
				s->sim(time,steps);
				break;
			case 2:
				enterStepLoop(s);
				break;
			case 3:
				cout<<"changing output file..."<<endl;
				exitCode = 2;
				break;
			case 4:
				cout<<endl<<endl;
				getPrintout(s);
				break;
			case 5:
				exit = true;
				break;
		}

		if(exit) break;
	}
	return exitCode;
}




void printStepMenu()
{
	cout<<"---------------------------"<<endl;
	cout<<"(0) step to next reaction"<<endl;
	cout<<"(1) step"<<endl;
	cout<<"(2) step 1 second"<<endl;
	cout<<"(3) step 10 seconds"<<endl;
	cout<<"(4) print stuff"<<endl;
	cout<<"(5) output observables to file"<<endl;
	cout<<"(6) return"<<endl;
}

void enterStepLoop(System *s)
{
	bool exit = false;
	while(true)
	{
		cout<<endl;
		printStepMenu();
		int selection = getInput(0,6);
		double time = 0;
		switch(selection) {
			case 0:
				s->singleStep();
			    break;
			case 1:
				cout<<"\nEnter duration of step: "<<endl;
				time = getInput(0.0);
				s->stepTo(s->getCurrentTime()+time);
				break;
			case 2:
				cout<<"Stepping for 1 second.  Please wait...  ";
				s->stepTo(s->getCurrentTime()+1);
				cout<<"done.\nCurrent Time is: "<<s->getCurrentTime()<<endl;
				break;
			case 3:
				cout<<"Stepping for 10 seconds.  Please wait...  ";
				s->stepTo(s->getCurrentTime()+10);
				cout<<"done.\nCurrent Time is: "<<s->getCurrentTime()<<endl;
				break;
			case 4:
				cout<<endl<<endl;
				getPrintout(s);
				break;
			case 5:
				cout<<"Writing to file...  ";
				s->outputAllObservableCounts();
				cout<<"done."<<endl;
				break;
			case 6:
				exit = true;
				break;
		}

		if(exit) break;
	}
}





void printSpecificMolecule(System *s)
{
	while(true)
	{
		cout<<"Select the MoleculeType:"<<endl;
		cout<<" (-1) none"<<endl;
		for(int m=0; m<s->getNumOfMoleculeTypes(); m++)
			cout<<" ("<<m<<") "<<s->getMoleculeType(m)->getName()<<" - has "<< s->getMoleculeType(m)->getMoleculeCount()<<" molecules."<<endl;
		int selection = getInput(-1,s->getNumOfMoleculeTypes()-1);
		if(selection==-1) break;

		while(true)
		{
			cout<<endl<<"Select a molecule (0 to "<<s->getMoleculeType(selection)->getMoleculeCount()-1<<", or -1 to exit):"<<endl;
			int selection2 = getInput(-1,s->getMoleculeType(selection)->getMoleculeCount()-1);
			if(selection2==-1) break;
			cout<<endl<<endl;
			s->getMoleculeType(selection)->getMolecule(selection2)->printDetails();
		}
	}
}
void printSpecificMoleculeByUid(System *s)
{
	while(true)
	{
		cout<<"Enter the molecule's unique id (or -1 to return):"<<endl;
		int selection = getInput(-1,Molecule::getUniqueIdCount()-1);
		if(selection==-1) break;

		cout<<endl;
		Molecule *mol = s->getMoleculeByUid(selection);
		if(mol!=0) mol->printDetails();
	}
}

void getPrintout(System *s)
{
	cout<<"What would you like to print out?"<<endl;
	cout<<" (0) all reactions"<<endl;
	cout<<" (1) all MoleculeTypes"<<endl;
	cout<<" (2) all index values and names"<<endl;
	cout<<" (3) all complexes"<<endl;
	cout<<" (4) all observables"<<endl;
	cout<<" (5) specific reaction"<<endl;
	cout<<" (6) specific molecule by MoleculeType"<<endl;
	cout<<" (7) specific molecule by unique ID"<<endl;
	int selection = getInput(0,8);
	switch(selection) {
		case 0:
			cout<<endl<<endl;
			s->printAllReactions();
		    break;
		case 1:
			cout<<endl<<endl;
			s->printAllMoleculeTypes();
			break;
		case 2:
			cout<<endl<<endl;
			s->printIndexAndNames();
			break;
		case 3:
			cout<<endl<<endl;
			// NETGEN -- redirect this call to the ComplexList object at s->allComplexes
			(s->getAllComplexes()).printAllComplexes();
			break;
		case 4:
			cout<<endl<<endl;
			s->printAllObservableCounts(s->getCurrentTime());
			break;
		case 5:
			cout<<endl<<endl;
			s->getReaction(0)->printDetails();
			break;
		case 6:
			cout<<endl<<endl;
			printSpecificMolecule(s);
			break;
		case 7:
			cout<<endl<<endl;
			printSpecificMoleculeByUid(s);
			break;
	}
}








int getInput(int min, int max)
{
	bool validInput = false;
	int intVal = -1;
	while(!validInput)
	{
		cout<<">";
		string input; cin>>input;
		try {
			intVal = NFutil::convertToInt(input);
			if(intVal>=min&&intVal<=max) validInput=true;
			else cout<<"   ---not an option, try again."<<endl;
		} catch (std::runtime_error e) {
			cout<<"   ---not an option, try again."<<endl;
		}
		cin.clear();
	}
	return intVal;
}

double getInput(double min)
{
	bool validInput = false;
	double doubVal = -1;
	while(!validInput)
	{
		cout<<">";
		string input; cin>>input;
		try {
			doubVal = NFutil::convertToDouble(input);
			if(doubVal>=min) validInput=true;
			else cout<<"   ---not an option, try again."<<endl;
		} catch (std::runtime_error e) {
			cout<<"   ---not an option, try again."<<endl;
		}
		cin.clear();
	}
	return doubVal;
}




