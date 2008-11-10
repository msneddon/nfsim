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
				cout<<"Warning: possible error in arguments.  You gave a '-' with a space"<<endl;
				cout<<"directly following. This is not a valid argument."<<endl<<endl;
				continue;
			}


			//See if the flag has some other input value that follows
			string sVal;
			if((a+1)<argc) {
				sVal=argv[a+1];
			}
			if(sVal.compare(0,1,"-")==0) {
				sVal = "";
			}

			//cout<<"found:  '"<<sFlag<<"' with arg: '"<<sVal<<"' "<<endl;
			if(argMap.find(sFlag)!=argMap.end()) {
				cout<<"Found two values for the same flag: '"<<sFlag<<"' so I am stopping"<<endl;
				return false;
			}

			argMap[sFlag] = sVal;
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



bool NFinput::createComplexOutputDumper(string paramStr, System *s, bool verbose)
{
	if(verbose) cout<<"Parsing dump complex (dc) flag: "<<paramStr<<endl;

	//First, extract out the times the user wants to dump the observables
	int b1 = paramStr.find_first_of('[');
	int b2 = paramStr.find_first_of(']');
	if(b1>=b2) { cout<<"Error in NFinput::creatComplexOutputDumper:, ']' was not found."<<endl; return false; }
	if(b1<0) { cout<<"Error in NFinput::creatComplexOutputDumper:, '[' was not found."<<endl; return false; }
	if(b2<0) { cout<<"Error in NFinput::creatComplexOutputDumper:, ']' was not found."<<endl; return false; }
	string timesStr=paramStr.substr(b1+1,b2-b1-1);

	//Create a vector storing the output times
	string numString=""; vector <double> outputTimes;
	if(verbose) cout<<"  scheduling complex output for times:";
	for(unsigned int i=0; i<timesStr.length(); i++) {
		if(timesStr.at(i)==';') {
			if(numString.size()==0) continue;
			try {
				double doubleVal = NFutil::convertToDouble(numString);
				if(outputTimes.size()>0) {
					if(doubleVal<=outputTimes.at(outputTimes.size()-1)) {
						cout<<"Error in NFinput::creatComplexOutputDumper: output times given \n";
						cout<<"must be monotonically increasing without any repeated elements.";
						return false;
					}
				}
				if(verbose) cout<<" "<<doubleVal<<";";
				outputTimes.push_back(doubleVal);
			} catch (std::runtime_error e) {
				cout<<"Warning in NFinput::creatComplexOutputDumper: could not parse time: '"<<numString<<"'"<<endl;
				cout<<"Ignoring that element."<<endl;
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
					cout<<"Error in NFinput::creatComplexOutputDumper: output times given ";
					cout<<"must be monotonically increasing without any repeated elements.";
					return false;
				}
			}
			if(verbose) cout<<" "<<doubleVal<<";";
			outputTimes.push_back(doubleVal);
		} catch (std::runtime_error e) {
			cout<<"Warning in NFinput::creatComplexOutputDumper: could not parse time: '"<<numString<<"'"<<endl;
			cout<<"Ignoring that element."<<endl;
		}
	}

	//Now we have the list of times we want to output.  So let's create the dumper and add it to the system.



	return true;
}



