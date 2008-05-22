/*! \mainpage NFsim: The Network Free Simulator
 *
 * \section intro_sec Overview
 *
 * The network free simulator is...
 *
 * \section install_sec Installation
 *
 * 
 * 
 * \subsection key Key Features
 *  NFsim can ...
 *  and it can also...
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 *  \section install_sec Developers
 * To begin developing and extending NFsim, the best place to start looking is in
 * the src/NFtest/simple_system directory. Here you'll find two files, simple_system.hh
 * and simple_system.cpp.  Together, this code specifies a simple enzymatic type reaction
 * that is completely hard coded.  This will give you an idea of the basic classes and
 * functions used to define, initialize, run, and output a simulaiton.  From there, you
 * can dive into the specific classes and functions that you need to work with.  Details
 * about how to run the simple_system example are given in these files.
 * 
 * All of the other main classes are defined in the NFcore namespace and are found in the NFcore
 * directory and the NFreactions directory.  The NFcore directory contains the basic structure
 * of the simulation engine while the NFreactions directory contains the classes associated with
 * actually executing rules and transforming molecules.  NFinput contains what's needed for
 * the xml parser (built using the TinyXML package) and the command line parser.  NFutil primarily
 * contains a nice implementation of the Mersenne Twister random number generator which should
 * be used for all random number generation in NFsim.  NFoutput is more sparse as it deals only
 * with handling the more complicated output required of groups and complexes.  (Basic outputting
 * is handled easily with the Observable class in the NFcore namespace).
 * 
 *  \section install_sec Authors & Acknowledgments
 * The NFsim code was written and developed by Michael Sneddon with help from James Faeder and
 * Thierry Emonet.  The core simulation algorithm is based on 
 * 
 * Special thanks to other members of the Emonet lab, particularly William Pontius, Garrit Jentsch, 
 * and Oleksii Sliusarenko for helpful feedback.  For questions or assistance with the code, please contact
 * michael.sneddon@yale.edu.
 * 
 * 
 */
#include "NFsim.hh"



#include <iostream>
#include <string>
using namespace std;


//!  Outputs an Ascii NFsim logo.
/*!
  @author Michael Sneddon
*/
void printLogo(int indent, string version);


//!  Outputs a friendly help message.
/*!
  @author Michael Sneddon
*/
void printHelp(string version);


//!  Main executable for the NFsim program.
/*!
  @author Michael Sneddon
*/
int main(int argc, const char *argv[])
{
	string versionNumber = "0.7";
	
	cout<<"starting NFsim v"+versionNumber+"..."<<endl<<endl;
	clock_t start,finish;
	double time;
	start = clock();
	///////////////////////////////////////////////////////////
	
	
	bool parsed = false;
	bool verbose = false;
	map<string,string> argMap;
	if(NFinput::parseArguments(argc, argv, argMap))
	{
		//First, find the arguments that we might use in any situation
		if(argMap.find("v")!=argMap.end()) verbose = true;
		
		
		//Handle the case of no parameters
		if(argMap.empty()) {
			cout<<endl<<"\tNo parameters given, so I won't do anything."<<endl;
			cout<<"\tIf you'd like help, pass me the -help flag."<<endl;
			parsed = true;
		}
		
		//Handle when the user asks for help!	
		else if (argMap.find("help")!=argMap.end()) 
		{
			printHelp(versionNumber);
			parsed = true;					
		}
		
		
		//Handle the case of reading from an xml file		
		else if (argMap.find("xml")!=argMap.end()) 
		{
			string filename = argMap.find("xml")->second;
			if(!filename.empty())
			{
				
				System *s = NFinput::initializeFromXML(argv[2], verbose);
			
				if(s!=NULL)
				{
					//Here we just run some stuff for testing... The output is just
					s->registerOutputFileLocation((s->getName()+".gdat").c_str());
					s->outputAllObservableNames();
					
					s->equilibriate(0,10);
					s->sim(50,50);
					s->stepTo(600);
					
					
					//s->sim(200,20);  // sim for 200 seconds, outputting 20 times
					s->printAllReactions();
					delete s;
				}
			}
			else
			{
				cout<<"You must specify an xml file to read."<<endl;
			}
			parsed = true;
		}
		
		
		//Handle the case of running a test		
		else if (argMap.find("test")!=argMap.end()) 
		{
			string test = argMap.find("test")->second;
			if(!test.empty())
			{
				cout<<"running test: '"+test+"'"<<endl;
				if(test=="simple_system")
				{
					NFtest_ss::run();
				}
				if(test=="transcription")
				{
					NFtest_transcription::run();
				}
			}
			else
			{
				cout<<"You must specify a test to run."<<endl;				
			}
					
			parsed = true;
		}

		
		
		//Finally, give the logo to anyone who wants it
		if (argMap.find("logo")!=argMap.end()) 
		{
			cout<<endl<<endl;
			printLogo(15,versionNumber);
			cout<<endl<<endl;
			cout<<"wow. that was awesome."<<endl;
			parsed = true;
		}
		
		
	}
	
	if(!parsed)
	{
		cout<<"Could not identify what you want to do."<<endl;
	}
	
	
	///////////////////////////////////////////////////////////
	// Finish and check the run time;
    finish = clock();
    time = (double(finish)-double(start))/CLOCKS_PER_SEC;
    cout<<endl<<"done.  Total run time: "<< time << "s"<<endl<<endl;
    return 0;
}














void printLogo(int indent, string version)
{
	string s;
	for(int i=0; i<indent; i++) s.append(" ");
	
	int space = 9-version.length();
	if(space<0) { 
		cout<<"\n\nCome on!!! you don't even know how to print out the NFsim logo!"<<endl;
		cout<<"What kind of code developer are you!!\n\n"<<endl;
	}
	string s2;
	for(int i=0; i<space; i++) s2.append(" ");
	cout<<s<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	cout<<s<<"%                                   %"<<endl;
	cout<<s<<"%     @@    @  @@@@@      v"<<version<<s2<<"%"<<endl;
	cout<<s<<"%     @ @   @  @                    %"<<endl;
	cout<<s<<"%     @  @  @  @@@@  ___            %"<<endl;
	cout<<s<<"%     @   @ @  @    /__  | |\\ /|    %"<<endl;
	cout<<s<<"%     @    @@  @    ___\\ | | v |    %"<<endl;
	cout<<s<<"%                                   %"<<endl;
	cout<<s<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
}



void printHelp(string version)
{
	cout<<"To run NFsim at the command prompt, use flags to specify what you want"<<endl;
	cout<<"to do.  Flags are given in this format in any order: \"-[flagName]\"."<<endl;
	cout<<"Some of the flags require an additional parameter.  For instance, the"<<endl;
	cout<<"-xml flag requires the filename of the xml file.  The format would look"<<endl;
	cout<<"something like: \"-xml modelFile.xml\"."<<endl;
	cout<<""<<endl;
	cout<<"Here are the list of possible flags:"<<endl;
	cout<<""<<endl;
	cout<<"  -help          well, you already know what this one does."<<endl;
	cout<<""<<endl;
	cout<<"  -xml           used to specify the input xml file to read.  the xml file"<<endl;
	cout<<"                 must be given directly after this flag."<<endl;
	cout<<""<<endl;
	cout<<"  -v             specify verbose output and print all kinds of extra things."<<endl;
	cout<<""<<endl;
	cout<<"  -test          used to specify a given preprogrammed test."<<endl;
	cout<<""<<endl;
	cout<<"  -logo          prints out the ascii NFsim logo."<<endl;
	cout<<""<<endl;
	cout<<""<<endl;
}










