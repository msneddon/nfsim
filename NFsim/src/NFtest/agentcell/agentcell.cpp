
#include <iostream>

#include "agentcell.hh"
#include "../../NFcore/NFcore.hh"
#include "../../NFinput/NFinput.hh"
#include "cell/cell.hh"
#include "cell/environment.hh"

using namespace NFcore;
using namespace std;


void runAgentCell(map<string,string> argMap, bool verbose)
{
	clock_t acstart,acfinish;
	double actime;


	cout<<"Running AgentCell Emulator."<<endl;
	if (argMap.find("xml")!=argMap.end())
	{
			System *s=0;
			string filename = argMap.find("xml")->second;
			if(!filename.empty())
			{
				cout<<"Reading file:"+ filename<<endl;
				int suggestedTraversalLimit = ReactionClass::NO_LIMIT;
				s=NFinput::initializeFromXML(filename,false,100000,verbose,suggestedTraversalLimit);
				s->setUniversalTraversalLimit(suggestedTraversalLimit);
				if (argMap.find("b")!=argMap.end()) s->setOutputToBinary();
			}
			if(s!=0) { cout<<" success."<<endl; }
			else { cout<<" failed, and quiting."<<endl; return; }


			//Make sure we have a CheY we can use for output
			Observable *cheY=s->getObservableByName("Yp");
			if(cheY==0) {
				cout<<"\n AgentCell Emulator requires an observable named 'Yp'"<<endl;
				cout<<" in order to accurately swim the cell.  Yp was not defined"<<endl;
				cout<<" in this file, so program will now quit."<<endl;
				return;
			}

			double eqTime = 10;
			double simTime = 10;
			double dt = 0.01;
			int nCells = 1;
			double motorThresh = 1800;


			eqTime = NFinput::parseAsDouble(argMap,"eq",eqTime);
			simTime = NFinput::parseAsDouble(argMap,"sim",simTime);

			dt = NFinput::parseAsDouble(argMap,"dt",dt);
			nCells = NFinput::parseAsInt(argMap,"nCells",(int)nCells);
			motorThresh = NFinput::parseAsInt(argMap,"motorThresh",(int)motorThresh);

			string oDir = "";
			if (argMap.find("oDir")!=argMap.end())
			{
				oDir= argMap.find("oDir")->second;
				if(oDir.empty()) {
					cout<<"Could not parse output directory, with flag 'oDir', given."<<endl;
					cout<<"quiting."<<endl; return;
				}
			} else {
				cout<<"No valid output directory, with flag 'oDir', given."<<endl;
				cout<<"quiting."<<endl; return;
			}

			cout<<"\n\n----- Run Options -----"<<endl;
			Environment *e;
			if(argMap.find("constEnvironment")!=argMap.end()) {
				cout<<"Environment:   "<<"Constant"<<endl;
				e= new ConstantEnvironment(0);
			} else if(argMap.find("linearEnvironment")!=argMap.end()) {
				double slope = pow(10,-8.0);
				slope = NFinput::parseAsDouble(argMap,"linearEnvironment",slope);
				double intercept = 0;
				intercept = NFinput::parseAsDouble(argMap,"zIntercept",intercept);
				cout<<"Environment:   "<<"Linear Gradient (slope: "<<slope<<", intercept: "<<intercept<<")"<<endl;
				e = new LinearEnvironment(slope,intercept);
			} else {
				double slope = pow(10,-8.0);
				cout<<"Environment:   "<<"Linear Gradient (default, slope: "<<slope<<", intercept: 0)"<<endl;
				e = new LinearEnvironment(slope,0);
			}

			cout<<"Eq time:       "<<eqTime<<"s"<<endl;
			cout<<"Sim time:      "<<simTime<<"s"<<endl;
			cout<<"dt:            "<<dt<<"s"<<endl;
			cout<<"Motor thresh:  "<<motorThresh<<" CheYp"<<endl;
			cout<<"Cell Count:    "<<nCells<<endl;
			cout<<"Output Dir:    "<<oDir<<endl;


			// More options
			//Motion:
			//PARAMETERS:
			// speed = 20 [um/sec]
			// sphere radius=1 [um = 1E-6 m], T=305 [K]
			// viscosity = 0.027 [g/(cm sec)], boltzman = 1.3807E-16 [erg/deg  (= g cm^2/sec^2)]
			// rotationaDiffusion = boltzman * T / (8*PI*viscosity*(1e-4*radius)^3)   = 0.0620577 [sec^-1]
			double cellSpeed = 20;
			double rotDiffusionConstant = 0.0620577;


			cout<<endl<<endl;
			string outputFolder;
			for(int i=0; i<nCells; i++) {
				acstart = clock();
				outputFolder = oDir+"/c" + NFutil::toString((i+1))+"/";

				string molTrajFileName = outputFolder+"molTraj.out";
				s->registerOutputFileLocation(molTrajFileName);
				s->prepareForSimulation();
				if (argMap.find("b")==argMap.end()) {
					s->outputAllObservableNames();
				}



				//Create the cell
				AgentCell *ac = new AgentCell(s,e,
						cellSpeed,
						rotDiffusionConstant,
						motorThresh,
						0,
						outputFolder);

				//Equilibrate first
				cout<<"Equilibrating Cell "<<i+1<<" ..."<<endl;
				ac->equilibriate(eqTime,dt);

				//Create a cell and run it
				cout<<"Simulating Cell "<<i+1<<" ..."<<endl;
				ac->stepTo(simTime,dt);


				acfinish = clock();
				actime = (double(acfinish)-double(acstart))/CLOCKS_PER_SEC;
				cout<<"done.  Elapsed CPU time: "<< actime << "s"<<endl<<endl;

				delete s;

				if(i<(nCells-1)) {
					cout<<endl<<endl<<"Reinitializing system..."<<endl;
					int suggestedTraversalLimit = ReactionClass::NO_LIMIT;
					s=NFinput::initializeFromXML(filename,false,100000,false,suggestedTraversalLimit);
					if (argMap.find("b")!=argMap.end()) {
						s->setOutputToBinary();
					}
					s->setUniversalTraversalLimit(suggestedTraversalLimit);
				}
			}

	}
}

