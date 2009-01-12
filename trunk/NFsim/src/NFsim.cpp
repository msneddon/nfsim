////////////////////////////////////////////////////////////////////////////////
//
//    NFsim: The Network Free Stochastic Simulator
//    A software framework for efficient simulation of biochemical reaction
//    systems with a large or infinite state space.
//
//    Copyright (C) 2008  Michael Sneddon, James Faeder, Thierry Emonet
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
//    For more information on NFsim, see http://emonet.biology.yale.edu/nfsim
//
////////////////////////////////////////////////////////////////////////////////

/*! \mainpage NFsim: The Network Free Stochastic Simulator
 *
 * \section intro_sec Overview
 *
 * NFsim is fully generalized stochastic reaction network simulator designed
 * to handle systems with a large (or even infinite) state space.  It has a
 * number of features that make it ideal for handling large and complex
 * biochemical systems, such as functionally defined rate laws and reactions
 * that depend on a local context.  NFsim is designed to operate with the BioNetGen
 * Language (http://bionetgen.org/).  The new version of BNG is able to
 * generate an XML encoded form of the BNG Language, which NFsim can take as input.
 *
 * For more details on setting up, running, and getting output from an NFsim simulation
 * see the User Manual.  The User Manual also has additional information for new
 * developers.  The Manual is available online along with examples here:
 * http://emonet.biology.yale.edu/NFsim
 *
 *
 *
 * \section key Command Line Argument List
 *
 *  Arguments can be provided to NFsim through the command line.  Below is a partial list
 *  of the available commands and a brief description of what they do.  For more details,
 *  see the NFsim user manual.
 *
 *  -help = outputs a helpful message to the console
 *
 *  -xml [filename] = specifies the xml file to read
 *
 *  -sim [Duration in sec] = specifies the length of time to simulate the system
 *
 *  -oSteps [num of steps] = specifies the number of times to output during the simulation
 *
 *  -eq [Duration in sec] = specifies the length of time to equilibriate before simulating
 *
 *  -o [filename] = specifies the name of the output file
 *
 *  -v = verbose output when reading an xml file and building a system
 *
 *  -b = output in binary (faster, but output is not human readable)
 *
 * -utl [integer] = universal traversal limit, see manual
 *
 * -notf = disables On the Fly Observables, see manual
 *
 *
 *  \section devel_sec Developers
 * To begin developing and extending NFsim, the best place to start looking is in
 * the src/NFtest/simple_system directory. Here you'll find two files, simple_system.hh
 * and simple_system.cpp.  Together, this code specifies a simple enzymatic type reaction
 * that is completely hard coded.  This will give you an idea of the basic classes and
 * functions used to define, initialize, run, and output a simulation.  From there, you
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
 * is handled easily with the System and Observable classes in the NFcore namespace).
 *
 * Another note for developers: class functions and member variables are generally well
 * commented in the header file in which they are declared.  So if you are lost in some source
 * file, and you think there aren't any comments, be sure to check the header file before
 * you ask for help.
 *
 *  \section author_sec Authors & Acknowledgments
 * The NFsim code was written and developed by Michael Sneddon with help from James Faeder and
 * Thierry Emonet.  James Faeder wrote the extended BioNetGen code that can output XML encodings
 * of the BNGL.
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

bool runRNFscript(map<string,string> argMap, bool verbose);
System *initSystemFromFlags(map<string,string> argMap, bool verbose);
bool runFromArgs(System *s, map<string,string> argMap, bool verbose);


//!  Main executable for the NFsim program.
/*!
  @author Michael Sneddon
*/
int main(int argc, const char *argv[])
{
	string versionNumber = "1.01";


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
		if(argMap.find("v")!=argMap.end()) {
			verbose = true;
		}
		if(argMap.find("seed")!= argMap.end()) {
			int seed = abs(NFinput::parseAsInt(argMap,"seed",0));
			NFutil::SEED_RANDOM(seed);
			cout<<"seeding random number generator with: "<<seed<<endl;
		}


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

		//If we are running from a RNF script file
		else if(argMap.find("rnf")!=argMap.end()) {
			cout<<"handling RNF file"<<endl;
			runRNFscript(argMap,verbose);
			parsed = true;
		}

		//Handle the case of reading from an xml file directly
		else if (argMap.find("xml")!=argMap.end())
		{
			System *s = initSystemFromFlags(argMap, verbose);
			if(s!=NULL) {
				runFromArgs(s,argMap,verbose);
			}
			parsed = true;
			delete s;
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



bool runRNFscript(map<string,string> argMap, bool verbose)
{
	//Step 1: open the file and initialize the argMap
	vector<string> commands;
	if(!NFinput::readRNFfile(argMap, commands, verbose)) {
		cout<<"Error when running the RNF script."<<endl;
		return false;
	}
	if(argMap.find("v")!=argMap.end()) verbose = true;


	//Step 2: using the argMap, set up the system
	System *s=initSystemFromFlags(argMap,verbose);
	if(s!=0) {
		s->prepareForSimulation();
		//Step 3: provided the system is set up correctly, run the RNF script
		bool output = NFinput::runRNFcommands(s,argMap,commands,verbose);
		delete s;
		return output;
	}

	return false;
}


System *initSystemFromFlags(map<string,string> argMap, bool verbose)
{
	if (argMap.find("xml")!=argMap.end())
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

				if (argMap.find("dump")!=argMap.end()) {
					if(!NFinput::createSystemDumper(argMap.find("dump")->second, s, verbose)) {
						cout<<endl<<endl<<"Error when creating system dump outputters.  Quitting."<<endl;
						delete s;
						return 0;
					}
				}

				if (argMap.find("utl")!=argMap.end()) {
					int utl = -1;
					utl = NFinput::parseAsInt(argMap,"utl",utl);
					s->setUniversalTraversalLimit(utl);
				}

				if (argMap.find("b")!=argMap.end()) {
					s->setOutputToBinary();
				}

				//Here we just run some stuff for testing... The output is just
				if (argMap.find("o")!=argMap.end()) {
					string outputFileName = argMap.find("o")->second;
					s->registerOutputFileLocation(outputFileName);
					s->outputAllObservableNames();
				} else {
					if(s->isOutputtingBinary())
						s->registerOutputFileLocation(s->getName()+"_nf.dat");
					else {
						s->registerOutputFileLocation(s->getName()+"_nf.gdat");
						s->outputAllObservableNames();
					}
				}


				//turn off on the fly calculation of observables
				if(argMap.find("notf")!=argMap.end()) {
					s->turnOff_OnTheFlyObs();
				}


				//Finally, return the system if we made it here without problems
				return s;
			}
			else  {
				cout<<"Couldn't create a system from your XML file.  I don't know what you did."<<endl;
				return 0;
			}
		}
	} else {
		cout<<"Couldn't create a system from your XML file.  No -xml [filename] flag given."<<endl;
	}
	return 0;
}


bool runFromArgs(System *s, map<string,string> argMap, bool verbose)
{

	//If requested, walk through the simulation
	if (argMap.find("walk")!=argMap.end()) {
		NFinput::walk(s);
	}
	//Otherwise, run as normal
	else
	{
		double eqTime = 0;
		double sTime = 10;
		int oSteps = 10;

		eqTime = NFinput::parseAsDouble(argMap,"eq",eqTime);
		sTime = NFinput::parseAsDouble(argMap,"sim",sTime);
		oSteps = NFinput::parseAsInt(argMap,"oSteps",(int)sTime);


		//Prepare the system for simulation!!
		s->prepareForSimulation();

		//Output some info on the system if we ask for it
		if(verbose) {
			cout<<"\n\nparse appears to be successful.  Here, check your system:\n";
			s->printAllMoleculeTypes();
			s->printAllReactions();
			cout<<"-------------------------\n";
		}

		cout<<endl<<endl<<endl<<"Equilibrating for :"<<eqTime<<"s.  Please wait."<<endl<<endl;
		s->equilibrate(eqTime);
		s->sim(sTime,oSteps);

		cout<<endl<<endl;
		s->printAllReactions();
	}
	return true;
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
	cout<<"to do.  Flags are given in this format in any order: \"-flagName\"."<<endl;
	cout<<"Some of the flags require an additional parameter.  For instance, the"<<endl;
	cout<<"-xml flag requires the filename of the xml file.  The format would look"<<endl;
	cout<<"something like: \"-xml modelFile.xml\".  Simulation output is dumped to"<<endl;
	cout<<"a file named: \"[modelName]_nf.gdat\" in the current directory by default."<<endl;
	cout<<""<<endl;
	cout<<"Here is a list of most of the possible flags:"<<endl;
	cout<<""<<endl;
	cout<<"  -help             well, you already know what this one does..."<<endl;
	cout<<""<<endl;
	cout<<"  -xml [filename]   used to specify the input xml file to read.  the xml"<<endl;
	cout<<"                    file must be given directly after this flag."<<endl;
	cout<<""<<endl;
	cout<<"  -o [filename]     used to specify the output file name."<<endl;
	cout<<""<<endl;
	cout<<"  -sim [time]       used to specify the length (in seconds) of a simulation"<<endl;
	cout<<"                    when running an xml file.  Fractional seconds are valid."<<endl;
	cout<<"                    for instance, you could use: -sim 525.50"<<endl;
	cout<<""<<endl;
	cout<<"  -eq [time]        used to specify the length (in seconds) to equilibrate the"<<endl;
	cout<<"                    system before running the simulation."<<endl;
	cout<<""<<endl;
	cout<<"  -oSteps [steps]   used to specify the number of times throughout the"<<endl;
	cout<<"                    simulation that observables will be outputted.  Must"<<endl;
	cout<<"                    be an integer value.  Default is to output once per"<<endl;
	cout<<"                    simulation second."<<endl;
	cout<<""<<endl;
	cout<<"  -v                specify verbose output to the console."<<endl;
	cout<<""<<endl;
	cout<<"  -b                use this flag to tell NFsim to output in binary (not ascii)"<<endl;
	cout<<""<<endl;
	cout<<"  -notf             tells NFsim to Not use On The Fly output.  Normally,"<<endl;
	cout<<"                    observables are computed On The Fly - that is they are"<<endl;
	cout<<"                    updated after every simulation step.  This is good if you"<<endl;
	cout<<"                    output frequently or have many molecules in your system."<<endl;
	cout<<"                    However, it can be faster to recompute observable counts"<<endl;
	cout<<"                    right before you output especially if you don't output"<<endl;
	cout<<"                    too often.  Use this flag to switch to recomputing at "<<endl;
	cout<<"                    every output step instead of using On The Fly output."<<endl;
	cout<<""<<endl;
	cout<<"  -ogf              output the value of all global functions."<<endl;
	cout<<""<<endl;
	cout<<"  -utl [integer]    sets the universal traversal limit"<<endl;
	cout<<""<<endl;
	cout<<"  -test             used to specify a given preprogrammed test. Some tests"<<endl;
	cout<<"                    include \"tlbr\" and \"simple_system\".  Tests do not read"<<endl;
	cout<<"                    in other command line flags"<<endl;
	cout<<""<<endl;
	cout<<"  -seed             used to specify the seed for the random number generator."<<endl;
	cout<<"                    This allows you to run the same simulation and get the"<<endl;
	cout<<"                    exact same results perhaps to compare performance"<<endl;
	cout<<""<<endl;
	cout<<"  -logo             prints out the ascii NFsim logo, for your viewing pleasure."<<endl;
	cout<<""<<endl;
	cout<<""<<endl;
}























//string filename = argMap.find("xml")->second;
//if(!filename.empty())
//{
//	//Create the system from the XML file
//	System *s = NFinput::initializeFromXML(filename, verbose);
//
//	if(s!=NULL)
//	{
//		//If requested, be sure to output the values of global functions
//		if (argMap.find("ogf")!=argMap.end()) {
//			s->turnOnGlobalFuncOut();
//		}
//
//		if (argMap.find("dump")!=argMap.end()) {
//			if(!NFinput::createSystemDumper(argMap.find("dump")->second, s, verbose)) {
//				cout<<endl<<endl<<"Error when creating system dump outputters.  Quitting."<<endl;
//				delete s;
//				exit(1);
//			}
//		}
//
//
//		if (argMap.find("utl")!=argMap.end()) {
//			int utl = -1;
//			utl = NFinput::parseAsInt(argMap,"utl",utl);
//			s->setUniversalTraversalLimit(utl);
//		}
//
//		if (argMap.find("b")!=argMap.end()) {
//			s->setOutputToBinary();
//		}
//
////		string mtName="mRNA";
////		MoleculeType *mt = s->getMoleculeTypeByName(mtName);
////		Outputter *o = new DumpMoleculeType("complex.out",s,mt);
////		s->addOutputter(o);
//
//		//Here we just run some stuff for testing... The output is just
//		if (argMap.find("o")!=argMap.end()) {
//			string outputFileName = argMap.find("o")->second;
//			s->registerOutputFileLocation(outputFileName);
//		} else {
//			if(s->isOutputtingBinary())
//				s->registerOutputFileLocation(s->getName()+"_nf.dat");
//			else {
//				s->registerOutputFileLocation(s->getName()+"_nf.gdat");
//				s->outputAllObservableNames();
//			}
//		}
//
//
//		//turn off on the fly calculation of observables
//		if(argMap.find("notf")!=argMap.end()) {
//			s->turnOff_OnTheFlyObs();
//		}
//
//
//		//If requested, walk through the simulation
//		if (argMap.find("walk")!=argMap.end()) {
//			NFinput::walk(s);
//		}
//		//Otherwise, run as normal
//		else
//		{
//			//Parameters (assigned first to their default values if these parameters
//			//are not explicitly given...
//			double eqTime = 0;
//			double sTime = 10;
//			int oSteps = 10;
//
//			eqTime = NFinput::parseAsDouble(argMap,"eq",eqTime);
//			sTime = NFinput::parseAsDouble(argMap,"sim",sTime);
//			oSteps = NFinput::parseAsInt(argMap,"oSteps",(int)sTime);
//
//
////						//Test for local functions ...
////						if(s->getName()=="localFunc2"){
////							cout<<"\n\n\n-------\nEntering local function test!!!"<<endl;
////							DORRxnClass::test1(s);
////							cout<<"yada"<<endl; exit(0);
////
////
////							MoleculeType *rec = s->getMoleculeTypeByName("Receptor");
////							MoleculeType *cheR = s->getMoleculeTypeByName("CheR");
////
////							Observable *obs_rCheR = s->getObservableByName("r_Cher");
////							vector <TemplateMolecule *> tmList;
////							TemplateMolecule::traverse(obs_rCheR->getTemplateMolecule(),tmList);
////							cout<<"\n----\n"<<tmList.size()<<endl;
////
////
////							vector <Observable *> obs;
////							TemplateMolecule * rec2 = new TemplateMolecule(rec);
////							rec2->addStateValue("m","2");
////							Observable * rec2obs = new Observable("RecM2", rec2);
////							obs.push_back(rec2obs);
////
////							obs.push_back(obs_rCheR);
////
////							vector <StateCounter *> sc;
////							StateCounter *scRecM = new StateCounter("RecMSum", rec, "m");
////							sc.push_back(scRecM);
////
////							vector <string> paramConstNames;
////							vector <double> paramConstValues;
////
////							LocalFunction *lf = new LocalFunction(s,
////									"simpleFunc",
////									"RecM2+(RecMSum*(5*5)/sqrt(10))",
////									obs,sc,paramConstNames,paramConstValues);
////
////							LocalFunction *lf2 = new LocalFunction(s,
////								"simpleFunc2",
////								"RecM2+(RecMSum+5)",
////								obs,sc,paramConstNames,paramConstValues);
////
////							lf->setEvaluationLevel(1);
////
////							//prepare!
////							lf->addTypeIMoleculeDependency(rec);
////
////							lf->printDetails();
////							//cout<<"reevaluating function on molecule"<<endl;
////							//lf->evaluateOn(rec->getMolecule(0));
////							//lf->printDetails();
////
////
////
////							///////////////////Testing DOR reactions....
////							TemplateMolecule *cherTemp = new TemplateMolecule(cheR);
////							TemplateMolecule *recTemp = new TemplateMolecule(rec);
////							recTemp->addStateValue("m","3");
////							vector <TemplateMolecule *> templates;
////							templates.push_back( recTemp );
////							templates.push_back( cherTemp );
////
////							TransformationSet *ts = new TransformationSet(templates);
////							ts->addLocalFunctionReference(recTemp,"Pointer1",LocalFunctionReference::SPECIES_FUNCTION);
////							ts->addLocalFunctionReference(recTemp,"Pointer2",LocalFunctionReference::SINGLE_MOLECULE_FUNCTION);
////							ts->addStateChangeTransform(recTemp,"m","1");
////							ts->finalize();
////
////
////
////							vector <LocalFunction *> lfList;
////							lfList.push_back(lf);
////							vector <string> lfPointerNameList;
////							lfPointerNameList.push_back("Pointer1");
////
////
////							DORRxnClass *r = new DORRxnClass("DorTest",0.5,ts,lfList,lfPointerNameList);
////							s->addReaction(r);
////
////							s->prepareForSimulation();
////
////
////
////							cout<<"\n\n\n\n\n------------**********-------------\n\n\n\n\n"<<endl;
////							rec->printDetails();
////							cheR->printDetails();
////
////							cout<<endl<<endl<<endl;
////							s->sim(10,10);
////
////
////							//rec->printAllMolecules();
////							//lf->printDetails();
////							//rec->getMolecule(0)->printDetails();
////							//cheR->getMolecule(0)->printDetails();
////
////
//////							cout<<"\n\n\n\n\n----\n\n";
//////							cout<<"rec mol count: "<<rec->getMoleculeCount()<<endl;
//////							for(int i=0; i<rec->getMoleculeCount(); i++) {
//////								r->printTreeForDebugging();
//////								cout<<"adding..."<<endl;
//////								rec->getMolecule(i)->printDetails();
//////								rec->getMolecule(i)->setComponentState("m",3);
//////								rec->getMolecule(i)->printDetails();
//////								r->directAddForDebugging(rec->getMolecule(i));
//////							}
//////							r->printTreeForDebugging();
////
////							cout<<"yada"<<endl; exit(0);
////							cout<<"ending test."<<endl;
////
////						}
////						else{
//
//
//			//Prepare the system for simulation!!
//			s->prepareForSimulation();
//
//			//Output some info on the system if we ask for it
//			if(verbose) {
//				cout<<"\n\nparse appears to be succussful.  Here, check your system:\n";
//				s->printAllMoleculeTypes();
//				s->printAllReactions();
//				cout<<"-------------------------\n";
//			}
//
//			cout<<endl<<endl<<endl<<"Equilibriating for :"<<eqTime<<"s.  Please wait."<<endl<<endl;
//			s->equilibriate(eqTime);
//			s->sim(sTime,oSteps);
//
//			cout<<endl<<endl;
//			s->printAllReactions();
////			s->dumpOutputters();
////						}
//		}
//		delete s;
//	}
//	else  {
//		cout<<"Couldn't create a system from your XML file.  I don't know what you did."<<endl;
//	}
//}
//else  {
//	cout<<"You must specify an xml file to read."<<endl;
//}













