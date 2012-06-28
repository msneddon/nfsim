


#include "NFinput.hh"





using namespace NFinput;
using namespace std;




// reads in the RNF file, and overwrites any argument it finds
bool NFinput::readRNFfile(map<string,string> &argMap, vector<string> &commands, bool verbose)
{
	if (argMap.find("rnf")!=argMap.end())
	{
		string filename = argMap.find("rnf")->second;
		if(!filename.empty())
		{
			//open up the file
			ifstream rnfFile(filename.c_str());
			string line; int lineCounter=0;
			bool started = false; bool ended = false;

			//If the file is open, start reading in the details...
			if (rnfFile.is_open()) {

				//Loop over the lines
				while (! rnfFile.eof() ) {

					//First, grab the line
					lineCounter++;
					getline(rnfFile,line);
					NFutil::trim(line);

					//identify and remove comments
					string::size_type pos = line.find_first_of("#");
					if(pos!=string::npos) {
						//cout<<"identified the comment."<<endl;
						line = line.substr(0,pos);
						NFutil::trim(line);
					}

					//Skip empty lines (which will also skip entire comment lines)
					if(line.size()==0) { continue; }

					//Identify the start tag...
					pos = line.find("begin");
					if(pos!=string::npos) {
					//	cout<<"found begin, all parameters have been found."<<endl;
						started = true;
						continue;
					}
					//Stop once we get to the end tag
					pos = line.find("end");
					if(pos!=string::npos) {
					//	cout<<"found end, ending the parse."<<endl;
						ended = true;
						break;
					}

					//If we are reading in the execution block, just add it to
					//the list - it will be handled later by the system
					if(started) {
					//	cout<<"["<<lineCounter<<"]: '"<<line<<"'"<<endl;
						commands.push_back(line);
						continue;
					}



					//Make sure that we have an argument
					pos = line.find_first_of("-");
					if(pos==string::npos || pos!=0) {
						cout<<"\nSyntax error in rnf file: '"<<filename<<"' on line ["<<lineCounter<<"]\n";
						cout<<"   >> "+line+"\n";
						cout<<"   This line gives no argument and is not commented out."<<endl;
						break;
					}

					string::size_type firstWhiteSpace = line.find_first_of(" \t");

					string argName = "";
					string argValue = "";
					if(firstWhiteSpace==string::npos) {
						argName = line.substr(1,line.size()-1);
						argValue = "";
					//	cout<<"got argname: '"<<argName<<"'"<<endl;
					} else {
						argName = line.substr(1,firstWhiteSpace-1);
						argValue = line.substr(firstWhiteSpace,line.size()-1);
						NFutil::trim(argValue);
					//	cout<<"got argname: '"<<argName<<"'"<<endl;
					//	cout<<"got argvalue: '"<<argValue<<"'"<<endl;
					}


					//Throw out some warnings if we are setting arguments in an improper way
					//that don't do anything
					if(argName=="sim" || argName=="eq") {
						cout<<"\nWarning in RNF file: '"<<filename<<"' on line ["<<lineCounter<<"]\n";
						cout<<"   >> "+line+"\n";
						cout<<"   The "<<argName<<" argument cannot be used in an RNF script.  Instead, add the\n";
						cout<<"   "<<argName<<" command to the execution block.  See the NFsim manual for details.\n"<<endl;
						continue;
					} else if(argName=="oSteps") {
						cout<<"\nWarning in RNF file: '"<<filename<<"' on line ["<<lineCounter<<"]\n";
						cout<<"   >> "+line+"\n";
						cout<<"   The oSteps argument cannot be used in an RNF script.  Instead, the\n";
						cout<<"   output steps should be specified with the sim command in the \n";
						cout<<"   execution block.  See the NFsim manual for details.\n"<<endl;
						continue;
					}

					//Make sure we haven't defined the argument yet
					if(argMap.find(argName)!=argMap.end()) {
						cout<<"\nSyntax error in rnf file: '"<<filename<<"' on line ["<<lineCounter<<"]\n";
						cout<<"   >> "+line+"\n";
						cout<<"   This argument was previously defined either in the actual command line\n";
						cout<<"   or directly in this RNF file."<<endl;
						break;
					}

					argMap[argName] = argValue;
				}

				if(!ended && started) {
					cout<<"\nSyntax error in RNF file: '"<<filename<<"' on line ["<<lineCounter<<"]\n";
					cout<<"   Execution block starting with 'begin' tag is not properly closed with an 'end' tag."<<endl;
					return false;
				}

				//Properly close the input stream
				rnfFile.close();
			}
			else {
				cout<<"\n   Error when opening the rnf file named: '"<<filename.c_str()<<"'"<<endl;
				cout<<"   File not found or is not accessible."<<endl;
				return false;
			}
		}
	} else {
		cout<<"\nCouldn't open the rnf file.  No -rnf [filename] flag given."<<endl;
		return false;
	}

	//if we got here, then we did something right...
	return true;
}



void echo(string command,System *s)
{
	int id1=command.find("echo");
	string message = command.substr(id1+4);
	//we could add outputting of parameters/variables if we parse them here
	NFutil::trim(message);
	cout<<message<<endl;
}

void print(string com,System *s)
{
	bool success=false;

	if(com.find("moleculeTypes")!=string::npos) {
		cout<<"\n"<<endl; s->printAllMoleculeTypes();
		success=true;
	}
	if(com.find("molecules")!=string::npos) {
		cout<<"\n"<<endl;
		for(int m=0; m<s->getNumOfMoleculeTypes(); m++)
			s->getMoleculeType(m)->printAllMolecules();
		success=true;
	}
	if(com.find("reactions")!=string::npos || com.find("rxns")!=string::npos ) {
		cout<<"\n"<<endl; s->printAllReactions();
		success=true;
	}
	if(com.find("functions")!=string::npos || com.find("funcs")!=string::npos) {
		cout<<"\n"<<endl; s->printAllFunctions();
		success=true;
	}
	if(com.find("parameters")!=string::npos || com.find("params")!=string::npos) {
		s->printAllParameters();
		success=true;
	}
	if(com.find("observables")!=string::npos || com.find("obs")!=string::npos) {
		cout<<"\n"<<endl; s->printAllObservableCounts(s->getCurrentTime());
		success=true;
	}


	if(!success) {
		cout<<"\nWarning in RNF execution command. \n";
		cout<<"   >> "+com+"\n";
		cout<<"   Could not figure out what you wanted to print.\n"<<endl;
	}
	cout<<"\n"<<endl;
}


void simulate(string command,System *s, bool verbose)
{
	int id1=command.find("sim");
	string times = command.substr(id1+3);
	NFutil::trim(times);

	string::size_type firstWhiteSpace = times.find_first_of(" \t");

	string simTime_str = "0"; double simTime = 0;
	string oSteps_str = "0";  int oSteps = 0;
	if(firstWhiteSpace==string::npos) {
		simTime_str = times.substr(0);
	} else {
		simTime_str = times.substr(0,firstWhiteSpace);
		oSteps_str = times.substr(firstWhiteSpace,times.size()-1);
	}

	NFutil::trim(simTime_str);
	NFutil::trim(oSteps_str);

	try {
		simTime=NFutil::convertToDouble(simTime_str);
		oSteps=NFutil::convertToInt(oSteps_str);
		cout<<">> ";
		s->sim(simTime,oSteps,verbose);
	} catch (std::runtime_error e) {
		cout<<"\nError in RNF execution command. \n";
		cout<<"   >> "+command+"\n";
		cout<<"   Could not convert simulation times or output steps to numbers.\n"<<endl;
		cerr<<e.what()<<endl;
	}

	cout<<"\n"<<endl;
}

void equilibrate(string command,System *s)
{
	int id1=command.find("eq");
	string times = command.substr(id1+3);
	NFutil::trim(times);

	string::size_type firstWhiteSpace = times.find_first_of(" \t");

	string simTime_str = "0"; double simTime = 0;
	string oSteps_str = "0";  int oSteps = 0;
	if(firstWhiteSpace==string::npos) {
		simTime_str = times.substr(0);
	} else {
		simTime_str = times.substr(0,firstWhiteSpace);
		oSteps_str = times.substr(firstWhiteSpace,times.size()-1);
	}

	NFutil::trim(simTime_str);
	NFutil::trim(oSteps_str);

	try {
		simTime=NFutil::convertToDouble(simTime_str);
		oSteps=NFutil::convertToInt(oSteps_str);
		cout<<">> equilibrating the system for "<<simTime<<" seconds"<<endl;
		s->equilibrate(simTime,oSteps);

	} catch (std::runtime_error e) {
		cout<<"\nError in RNF execution command. \n";
		cout<<"   >> "+command+"\n";
		cout<<"   Could not convert eq times or output steps to numbers.\n"<<endl;
		cerr<<e.what()<<endl;
	}

	cout<<"\n"<<endl;
}


void setParameter(string command, System *s) {

	int id1=command.find("set");
	string paramString = command.substr(id1+3);
	NFutil::trim(paramString);

	string::size_type firstWhiteSpace = paramString.find_first_of(" \t");

	string paramName = "";
	string paramValue_str = "0";  double paramValue = 0;
	if(firstWhiteSpace==string::npos) {
		cout<<"Could not update parameter: '"<<paramName<<"'! No value given!"<<endl;
		return;
	} else {
		paramName = paramString.substr(0,firstWhiteSpace);
		paramValue_str = paramString.substr(firstWhiteSpace,paramString.size()-1);
	}

	NFutil::trim(paramName);
	NFutil::trim(paramValue_str);
	try {
		paramValue=NFutil::convertToDouble(paramValue_str);
		cout<<"Trying to set paramater: '"<<paramName<<"' to value: "<<paramValue<<endl;
		s->setParameter(paramName,paramValue);

	} catch (std::runtime_error e) {
		cout<<"\nError in RNF execution command. \n";
		cout<<"   >> "+command+"\n";
		cout<<"   Could not convert parameter value to a valid number.\n"<<endl;
		cerr<<e.what()<<endl;
		return;
	}
}



bool NFinput::runRNFcommands(System *s, map<string,string> &argMap, vector<string> &commands, bool verbose)
{
	cout<<"\n\nrunning RNF commands\n-----------------"<<endl;
	string com = "";
	for(int c=0; c<(int)commands.size(); c++) {
		com = commands.at(c);
	//	cout<<"Parsing command: "<<com<<endl;


		if(com.find("echo")!=string::npos) {
			echo(com,s);
			continue;
		} else if(com.find("print")!=string::npos) {
			print(com,s);
			continue;
		} else if(com.find("eq")!=string::npos) {
			equilibrate(com,s);
			continue;
		} else if(com.find("sim")!=string::npos) {
			simulate(com,s,verbose);
			continue;
		} else if(com.find("set")!=string::npos) {
			setParameter(com,s);
			continue;
		}  else if(com.find("update")!=string::npos) {
			s->updateSystemWithNewParameters();
			continue;
		} else {
			cout<<"could not figure out what you wanted to do for command:\n";
			cout<<com<<endl;
		}

	}


	return false;
}


























