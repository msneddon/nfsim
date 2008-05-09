/*! \mainpage NFsim: The Network Free Simulator (beta v7)
 *
 * \section intro_sec Introduction
 *
 * The network free simulator is...
 *
 * \section install_sec Installation
 *
 * \subsection key Key Features
 *  NFsim can ...
 *  and it can also...
 * 
 */
#include "NFsim.hh"



#include <iostream>
#include <string>
using namespace std;




void printLogo(int indent, string version);




int main(int argc, char *argv[])
{
	cout<<"starting NFsim V7..."<<endl<<endl;
	clock_t start,finish;
	double time;
	start = clock();
	///////////////////////////////////////////////////////////
	
	
	//First check and parse the parameters
		if(argc==1)
			cout<<endl<<"No parameters given, so I won't do anything."<<endl;	
		else if(argc>1)
		{
			//Check if we are running a test
			if(strcmp(argv[1],"-test")==0)	
			{
				if(argc>2)
				{
					//Determine which test to run
//					if(strncmp(argv[2],"TLBR",4)==0)
//						NFtest_TLBR::run(argc, argv);
//					else if(strncmp(argv[2],"testCompare",11)==0)
//						NFtest_compare::run();
//					else if(strcmp(argv[2],"transformation")==0)
//						NFtest_transformations::run();
					if(strcmp(argv[2],"simple_system")==0)
						NFtest_simple_system::run();
					else
						cout<<"Could not identify test: "<<argv[2]<<endl;
				} 
				else
					cout<<"You must specify which test to run."<<endl;	
			}
			
			//Check if we are running an actual system run
			else if(strcmp(argv[1],"-run")==0)	
			{
				if(argc>2)
				{
					if(strcmp(argv[2],"an")==0)
						cout<<"an..  not in v6"<<endl;
					//	run_AN_system(argc, argv);
					else if(strcmp(argv[2],"ng")==0)
					{
							//NG::run(argc, argv);
					}
					else
						cout<<"Could not identify run: "<<argv[2]<<endl;
				}
				else
					cout<<"You must specify which model to run."<<endl;	
			}
			
			
			else if(strcmp(argv[1],"-swim")==0)	
			{
				cout<<"swimming.."<<endl;
//				char * cellFileName1 = "/home/msneddon/Desktop/NF_output/chemotaxisOut/testRun2/cell_0/cellout.txt";
//				char * cellFileName2 = "/home/msneddon/Desktop/NF_output/chemotaxisOut/testRun2/cell_1/cellout.txt";
//				char * cellFileName3 = "/home/msneddon/Desktop/NF_output/chemotaxisOut/testRun2/cell_2/cellout.txt";
//				char * cellFileName4 = "/home/msneddon/Desktop/NF_output/chemotaxisOut/testRun2/cell_3/cellout.txt";
//				char * cellFileName5 = "/home/msneddon/Desktop/NF_output/chemotaxisOut/testRun2/cell_4/cellout.txt";
//				
//				double startTime = 0;
//				double startX=50, startY=0, startZ=0;
//				bool shouldOutputStats = true;
//				
//				ChemotacticCell *c1 = new ChemotacticCell(0,"", "" ,"","",cellFileName1,startTime, startX, startY,startZ,false);
//				ChemotacticCell *c2 = new ChemotacticCell(1,"", "" ,"","", cellFileName2,startTime, startX, startY,startZ,false);
//				ChemotacticCell *c3 = new ChemotacticCell(2,"", "" ,"","", cellFileName3,startTime, startX, startY,startZ,false);
//				ChemotacticCell *c4 = new ChemotacticCell(3,"", "" ,"","", cellFileName4,startTime, startX, startY,startZ,false);
//				ChemotacticCell *c5 = new ChemotacticCell(4,"", "" ,"","", cellFileName5,startTime, startX, startY,startZ,false);
//				
//				double eTime = 500;
//				c1->equilibriate(eTime);
//				c2->equilibriate(eTime);
//				c3->equilibriate(eTime);
//				c4->equilibriate(eTime);
//				c5->equilibriate(eTime);
//				
//				double step = 100.01;
//				double n_steps = 20;
//				for(int i=1; i<=n_steps; i++)
//				{
//					c1->stepTo(step*i,0.01);
//					c2->stepTo(step*i,0.01);
//					c3->stepTo(step*i,0.01);
//					c4->stepTo(step*i,0.01);
//					c5->stepTo(step*i,0.01);
//					cout<<endl;
//				}
//				delete c1;
//				delete c2;
//				delete c3;
			}
			
			else if(strcmp(argv[1],"-xml")==0)	
			{
				if(argc>2)
				{
					System *s = NFinput::initializeFromXML(argv[2]);
					
					//Here we just run some stuff for testing... The output is just
					s->registerOutputFileLocation((s->getName()+".gdat").c_str());
					s->outputAllObservableNames();
					s->sim(200,100);  // sim for 20 seconds, outputting 100 times
					s->printAllReactions();
					delete s;
				}
				else
					cout<<"You must specify an xml file to read."<<endl;
				
				
			}
			
			
			else if(strcmp(argv[1],"-logo")==0)	
			{
				cout<<endl<<endl;
				printLogo(15,"0.7");
				cout<<endl<<endl;
				cout<<"wow. that was awesome."<<endl;
			}
			else if(strcmp(argv[1],"-help")==0)	
			{
				cout<<"Welcome to NFsim: The Network Free Simulator"<<endl;
				cout<<"help message to be inserted here."<<endl;
			}
			
			
			else
			{
				cout<<"Cannot identify what you want to do."<<endl;
				
				
			}
			
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












