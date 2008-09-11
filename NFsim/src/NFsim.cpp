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
 * arguements accecpted:
 * -v = verbose output
 * -sim = length of time in seconds of a simulation
 * -ogf = output global function values at each point of the output
 * -help = output a help message
 * -xml = read an xml file
 * -test = run a predefined test
 * -logo = print the nfsim logo
 * -eq = equilibriate for a length of time (in seconds)
 * -oSteps = number of times 
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
 * Thierry Emonet.
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
#include <time.h>

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
	string versionNumber = "0.71";
	
	
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
		else if (argMap.find("help")!=argMap.end())  {
			printHelp(versionNumber);
			parsed = true;					
		}
		
		
		//Handle the case of reading from an xml file		
		else if (argMap.find("xml")!=argMap.end()) 
		{
			string filename = argMap.find("xml")->second;
			if(!filename.empty())
			{
				//Create the system from the XML file
				System *s = NFinput::initializeFromXML(filename, verbose);
			
				if(s!=NULL)
				{
					//If requested, be sure to output the values of global functions
					if (argMap.find("ogf")!=argMap.end()) {
						s->turnOnGlobalFuncOut();
					}
					
					if (argMap.find("b")!=argMap.end()) {
						s->setOutputToBinary();
					}
					
					//Here we just run some stuff for testing... The output is just
					if (argMap.find("o")!=argMap.end()) {
						string outputFileName = argMap.find("o")->second;
						s->registerOutputFileLocation(outputFileName);
					} else {
						if(s->isOutputtingBinary())
							s->registerOutputFileLocation(s->getName()+"_nf.dat");
						else {
							s->registerOutputFileLocation(s->getName()+"_nf.gdat");
							s->outputAllObservableNames();
						}
					}
					
					
					//If requested, walk through the simulation
					if (argMap.find("walk")!=argMap.end()) {
						NFinput::walk(s);
					}
					//Otherwise, run as normal
					else
					{
						//Parameters (assigned first to thier default values if these parameters
						//are not explicitly given...
						double eqTime = 0;
						double sTime = 10;
						int oSteps = 10;
						
						eqTime = NFinput::parseAsDouble(argMap,"eq",eqTime);
						sTime = NFinput::parseAsDouble(argMap,"sim",sTime);
						oSteps = NFinput::parseAsInt(argMap,"oSteps",(int)sTime);
						
						cout<<endl<<endl<<endl<<"Equilibriating for :"<<eqTime<<"s.  Please wait."<<endl<<endl;
						s->equilibriate(eqTime);
						s->sim(sTime,oSteps);
						
						cout<<endl<<endl;
						s->printAllReactions();
					}
					delete s;
				}
				else  {
					cout<<"Couldn't create a system from your XML file.  I don't know what you did."<<endl;
				}
			}
			else  {
				cout<<"You must specify an xml file to read."<<endl;
			}
			parsed = true;
		}
		
		
		//Handle the case of running a test		
		else if (argMap.find("test")!=argMap.end()) 
		{
			string test = argMap.find("test")->second;
			bool foundATest = false;
			if(!test.empty())
			{
				cout<<"running test: '"+test+"'"<<endl;
				if(test=="simple_system") {
					NFtest_ss::run();
					foundATest=true;
				}
				if(test=="transcription") {
					NFtest_transcription::run();
					foundATest=true;
				}
				if(test=="tlbr") {
					NFtest_tlbr::run(argMap);
					foundATest=true;
				}
				if(test=="mathFuncParser") {
					FuncFactory::test();
					foundATest=true;
				}
				
				if(!foundATest) {
					cout<<"  That test could not be identified!!  Skipping!"<<endl;
				}
					
			}
			else {
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
	
	
	if(!parsed) {
		cout<<"Could not identify what you want to do.  Try running the -help flag for advice."<<endl;
	}
	
	
	///////////////////////////////////////////////////////////
	// Finish and check the run time;
    finish = clock();
    time = (double(finish)-double(start))/CLOCKS_PER_SEC;
    cout<<endl<<"done.  Total CPU time: "<< time << "s"<<endl<<endl;
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
	cout<<"something like: \"-xml modelFile.xml\".  Simulation output is dumped to"<<endl;
	cout<<"a file named: \"[modelName]_nf.gdat\" in the current directory by default."<<endl;
	cout<<""<<endl;
	cout<<"Here is a list of most of the possible flags:"<<endl;
	cout<<""<<endl;
	cout<<"  -help          well, you already know what this one does..."<<endl;
	cout<<""<<endl;
	cout<<"  -xml           used to specify the input xml file to read.  the xml file"<<endl;
	cout<<"                 must be given directly after this flag."<<endl;
	cout<<""<<endl;
	cout<<"  -o             used to specify the output file name."<<endl;
	cout<<""<<endl;
	cout<<"  -sim           used to specify the length (in seconds) of a simulation when"<<endl;
	cout<<"                 running an xml file.  Fractional seconds are valid.  For"<<endl;
	cout<<"                 instance, you could use: -sim 525.50"<<endl;
	cout<<""<<endl;
	cout<<"  -eq            used to specify the length (in seconds) to equilibriate the"<<endl;
	cout<<"                 system before running a simulation."<<endl;
	cout<<""<<endl;
	cout<<"  -oSteps        used to specify the number of times throughout the simulation"<<endl;
	cout<<"                 that observables will be outputted.  Must be an integer value."<<endl;
	cout<<"                 Default is to output once per simulation second."<<endl;
	cout<<""<<endl;
	cout<<"  -v             specify verbose output and print all kinds of extra things."<<endl;
	cout<<""<<endl;
	cout<<"  -test          used to specify a given preprogrammed test. Some tests"<<endl;
	cout<<"                 include \"tlbr\" and \"simple_system\".  Tests do not read in"<<endl;
	cout<<"                 other command line flags, so your -sim flag won't do anything."<<endl;
	cout<<""<<endl;
	cout<<"  -logo          prints out the ascii NFsim logo, for your viewing pleasure."<<endl;
	cout<<""<<endl;
	cout<<""<<endl;
}










