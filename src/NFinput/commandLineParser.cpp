/*
 * commandLineParser.cpp
 *
 *  Created on: Oct 21, 2008
 *      Author: msneddon
 */

#include "NFinput.hh"





using namespace NFinput;
using namespace std;


bool NFinput::parseArguments(int argc, const char *argv[], map<string,string> &argMap)
{
	for(int a=1; a<argc; a++)
	{
		string s(argv[a]);

		//Find the strings that start with a minus, these are the flags
		if(s.compare(0,1,"-")==0)
		{
			string sFlag = s.substr(1,s.size()-1);
			if(sFlag.empty()) {
				cout<<"   Error in the command line arguments.  You gave a '-' with a space"<<endl;
				cout<<"   directly following. This is not a valid argument."<<endl<<endl;
				return false;
			}

			if(sFlag.compare(0,1,"-")==0){
				sFlag = s.substr(1,sFlag.size()-1);
			}
			if(sFlag.empty()) {
				cout<<"   Error in the command line arguments.  You gave a '--' with a space"<<endl;
				cout<<"   directly following. This is not a valid argument."<<endl<<endl;
				return false;
			}


			//See if the flag has some other input value that follows
			string sVal;
			if((a+1)<argc) {
				sVal=argv[a+1];
			}
			if(sVal.compare(0,1,"-")==0) {
				sVal = "";
			} else {
				a++;
			}

			//cout<<"found:  '"<<sFlag<<"' with arg: '"<<sVal<<"' "<<endl;
			if(argMap.find(sFlag)!=argMap.end()) {
				cout<<"Found two values for the same command line flag: '"<<sFlag<<"' so I am stopping"<<endl;
				return false;
			}

			argMap[sFlag] = sVal;
		}
		else
		{
			cout<<"   Warning when parsing command line arguments.  Valid arguments are preceded by a standard dash, as in '-logo'."<<endl;
			cout<<"\n   This argument: '"<< s <<"'"<<endl;
			cout<<"   did not begin with a proper dash and was ignored."<<endl<<endl;
		}
	}
	return true;
}



int NFinput::parseAsInt(map<string,string> &argMap,string argName,int defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			int intVal = NFutil::convertToInt(strVal);
			return intVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as an integer."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}


void NFinput::parseAsCommaSeparatedSequence(map<string,string> &argMap,string argName,vector<int> &sequence)
{
	if(argMap.find(argName)!=argMap.end()) {
		string argString = argMap.find(argName)->second;
		try {

			vector <string> numberStrings;
			numberStrings.push_back("");
			for(unsigned int i=0; i<argString.length(); i++)
			{
				if(argString.at(i)=='\"' || argString.at(i)==' ') { continue; }
				if(argString.at(i)==',') { numberStrings.push_back(""); continue; }
				numberStrings.at(numberStrings.size()-1) = numberStrings.at(numberStrings.size()-1) + argString.at(i);
			}

			for(unsigned int i=0; i<numberStrings.size(); i++) {
				sequence.push_back(NFutil::convertToInt(numberStrings.at(i)));
			}

		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+argString+"' as a comma separated "+
					"integer sequence, so I'm quitting."<<endl;
			exit(1);
		}
	}
}

double NFinput::parseAsDouble(map<string,string> &argMap,string argName,double defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			double doubleVal = NFutil::convertToDouble(strVal);
			return doubleVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as a double."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}



bool NFinput::parseSequence(string numString, vector <double> &outputTimes)
{
	double startVal=0, stepVal=1, endVal=0;
	try {

		string::size_type c1 = numString.find_first_of(':');
		if(c1!=string::npos) {
			string::size_type c2 = numString.find_first_of(':',c1+1);
			if(c2!=string::npos) {
					startVal= NFutil::convertToDouble(numString.substr(0,c1));
					stepVal= NFutil::convertToDouble(numString.substr(c1+1,c2-c1-1));
					endVal= NFutil::convertToDouble(numString.substr(c2+1));

			} else {
					startVal= NFutil::convertToDouble(numString.substr(0,c1));
					endVal= NFutil::convertToDouble(numString.substr(c1+1));
			}
		}

	} catch(std::runtime_error e) {
		return false;
	}

	if(startVal>endVal) {
		cout<<"Error: start value of sequence must be <= end value."<<endl;
		return false;
	} else if(stepVal<=0) {
		cout<<"Error: step value of sequence must be >0."<<endl;
		return false;
	} else {
		if(outputTimes.size()>=1)
			if(startVal<=outputTimes.at(outputTimes.size()-1)) {
				cout<<"\n\nError in NFinput::creatComplexOutputDumper: output times given \n";
				cout<<"must be monotonically increasing without any repeated elements.";
				return false;
			}
		//Only if everything went as planned to we then add the output steps accordingly
		for(double d=startVal; d<=endVal; d+=stepVal) {
			outputTimes.push_back(d);
		}
		return true;
	}
	return true;
}




bool NFinput::createSystemDumper(string paramStr, System *s, bool verbose)
{
	if(verbose) cout<<"Parsing system dump flag: "<<paramStr<<"\n";

	//First, extract out the times the user wants to dump the observables
	int b1 = paramStr.find_first_of('[');
	int b2 = paramStr.find_first_of(']');
	if(b1>b2) { cout<<"Error in NFinput::createSystemDumper:, ']' was found before '['."<<endl; return false; }
	if(b1<0 && b2<0) { cout<<"Error in NFinput::createSystemDumper:, enclosing brackets '[' and ']' were not found."<<endl; return false; }
	if(b1<0) { cout<<"Error in NFinput::createSystemDumper:, '[' was not found."<<endl; return false; }
	if(b2<0) { cout<<"Error in NFinput::createSystemDumper:, ']' was not found."<<endl; return false; }
	string timesStr=paramStr.substr(b1+1,b2-b1-1);



	//Now find path to the folder we are dumping to
	string pathToFolder = paramStr.substr(b2+1);
	if(pathToFolder.size()!=0) {
		string::size_type arrowPos = pathToFolder.find("->");
		if(arrowPos!=string::npos) {
			pathToFolder = pathToFolder.substr(arrowPos+2);
		} else {
			cout<<"Warning: path to folder ("+pathToFolder+") is not written correctly."<<endl;
			cout<<"Should be written as: [t1;t2;]->/path/to/folder/"<<endl;
			cout<<"no system dumps were scheduled."<<endl;
			return true;
		}
	}


	//Create a vector storing the output times
	string numString=""; vector <double> outputTimes;
	if(verbose) { cout<<"  scheduling system dumps at simulation times:"; }
	if(pathToFolder.size()>0) { cout<<"scheduling system dumps to directory ("+pathToFolder+")"<<endl; }
	else { cout<<"scheduling system dumps to directory (.)"<<endl; }
	for(unsigned int i=0; i<timesStr.length(); i++) {
		if(timesStr.at(i)==';') {
			if(numString.size()==0) continue;
			try {
				double doubleVal = NFutil::convertToDouble(numString);
				if(outputTimes.size()>0) {
					if(doubleVal<=outputTimes.at(outputTimes.size()-1)) {
						cout<<"\n\nError in NFinput::creatComplexOutputDumper: output times given \n";
						cout<<"must be monotonically increasing without any repeated elements.";
						return false;
					}
				}
				outputTimes.push_back(doubleVal);
			} catch (std::runtime_error e) {
				//could not parse it directly as a double, so first we also have to try parsing
				//it as a matlab style sequence...
				bool success = parseSequence(numString, outputTimes);
				if(!success) {
					cout<<"\nWarning in NFinput::creatComplexOutputDumper: could not parse time: '"<<numString<<"'"<<endl;
					cout<<"Ignoring that element."<<endl;
				}
			}
			numString="";
			continue;
		}
		numString += timesStr.at(i);
	}
	if(numString.size()!=0) {
		try {
			double doubleVal = NFutil::convertToDouble(numString);
			if(outputTimes.size()>0) {
				if(doubleVal<=outputTimes.at(outputTimes.size()-1)) {
					cout<<"\nError in NFinput::creatComplexOutputDumper: output times given ";
					cout<<"must be monotonically increasing without any repeated elements.";
					return false;
				}
			}
			outputTimes.push_back(doubleVal);
		} catch (std::runtime_error e) {
			//try to parse it as a matlab style sequence
			bool success = parseSequence(numString, outputTimes);
			if(!success) {
				cout<<"\nWarning in NFinput::creatComplexOutputDumper: could not parse time: '"<<numString<<"'"<<endl;
				cout<<"Ignoring that element."<<endl;
			}
		}
	}
	if(verbose) {
		if(outputTimes.size()==0) {
			cout<<" none given.";
		} else {
			for(unsigned int i=0; i<outputTimes.size(); i++)
				cout<<" "<<outputTimes.at(i)<<";";
		}
		cout<<endl;
	}




	//Here is where we actually create the system dumper
	DumpSystem *ds = new DumpSystem(s, outputTimes, pathToFolder, verbose);
	s->setDumpOutputter(ds);
	return true;

}

