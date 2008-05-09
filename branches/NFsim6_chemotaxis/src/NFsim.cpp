/*! \mainpage The Network Free Simulator (beta v6)
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
#include <stdio.h> 
#include <exception>

using namespace std;


//#include "NFmodel/chemotaxis/util/util.hh";
#include "NFmodel/chemotaxis/ng_swim/cell/cell.hh"
#include "NFmodel/chemotaxis/ng_swim/environment/environment.hh"
#include "NFmodel/chemotaxis/ng_swim/population/population.hh"
using namespace NG_swim;





#include <math.h>







class puffer {
	
	
public:
	
	puffer(double threshold) {};
	~puffer() {};
	void addPuff(double N, double x, double y, double z, double t0, double D);
	double getConc(double x1, double y1, double z2);
	
protected:
	//vector <puff> puffs;
	
	
};


void puffer::addPuff(double N, double x, double y, double z, double t0, double D) {}



class puff {
	
	public:
		puff(double N, double x, double y, double z, double t0, double D);
		~puff() {};
		
		double conc(double x1, double y1, double z1, double t);
		
		double x;
		double y;
		double z;
		double t0;
		double dt;
		double N;
		double D;
		
};

puff::puff(double N, double x, double y, double z, double t0, double D)
{
	this->N=N;
	this->x=x;
	this->y=y;
	this->z=z;
	this->t0=t0;
	this->D=D;
}

double puff::conc(double x1, double y1, double z1, double t)
{
	double dx = x1-x; double dy = y1-y; double dz = z1-z;
	double r_squared = dx*dx+dy*dy+dz*dz;
	double prefactor = (N/(pow(4.0*NFutil::PI*D*(t-t0),3.0/2.0)));
	return prefactor*exp((-r_squared)/(4*D*t-t0));
}




vector < puff * > puffs;


double gaus(double r, double t, double N)
{
	double D = 5;
	double c = (N/(pow(4*NFutil::PI*D*t,3./2.))) *exp((-r*r)/(4*D*t));
	return c;
}







#include "NFmodel/chemotaxis/ng_creator/ng_creator.hh"


int main(int argc, char *argv[])
{
/*	puff p(5000, 0,0,0, 0, 890);
	
	double t=0.04;
	for (double x=0; x<10;x+=0.1)
	{	
		cout<<(p.conc(x,0,0,t)/(NFutil::NA*1e-15))/1000<<endl;
	}
	
	NFutil::NA;
	exit(1);
	*/
	
		/*	NGparam p;
			p.setLiteSystem();
			p.setToActivityOutput();
			System * s = create("liteSys",p);
				
			//s->equilibriate(200,5);
			
			int n_values = 1;
				double * ligandConcentation = new double[n_values];
				ligandConcentation[0] = 0.00001; //0.01mM
				
				int ti = 600; int samples = 60;
				
				//s->equilibriate(1000,100);
				//s->sim(1000,1000);
				s->sim(300,30);
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples);
				    
				ligandConcentation[0] = 0.0001; //0.1mM
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples); 
				    
				ligandConcentation[0] = 0.001;  //1mM
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples);
				    
				ligandConcentation[0] = 0.01;   //10mM
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples);
				    
				ligandConcentation[0] = 0.1; // 100 mM
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples);
				
				ligandConcentation[0] = 1; 
				s->updateAllGroupProperty(ligandConcentation, n_values);
				s->sim(ti,samples);
				
			s->printAllReactions();
			delete s;
			exit(0);
	
	*/
	if(false) {
	bool noise = false;
	
	if(noise)
	{
		NGparam p;
		p.setFullLite();
	
		p.setToNoiseOutput();
	
		System * s = create("testSys",p);
	
		s->equilibriate(2700,30);
		s->sim(10000,200000);
		s->printAllReactions();
		delete s;
	}
	else
	{
		NGparam p;
		p.setFullLite();
		
		p.setInitReceptorMethToZero();
		p.setToActivityOutput();
		
		
		
		
		
		System * s = create("testSys",p);
			
		//s->equilibriate(200,5);
		
		s->sim(20000,2000);
			
			
		s->printAllReactions();
		delete s;
	}
	
	cout<<"done."<<endl;
	exit(0);
	
	}
	

	
//	//shape = 4 scale = 18.32045567939674
//	ChemotaxisUtil::GammaSampler *g = new ChemotaxisUtil::GammaSampler(4.,1./18.32045567939674,-4.606176605609249);
//	for(int i=0; i<5000; i++)
//	{
//		cout<<g->gen(0,180)<<endl;
//		
//	}
	
	
	
//	const char * outputDirectory = "/home/msneddon/Desktop/swim_output/pop2_cap";
//	
//	int cellCount=1;
//	
//	double equilibriateTime = 500;
//	double simTime = 3000;
//	double popDt = 5;
//	double cellDt = 0.01;
//	
//	Environment *e = new CapillaryEnvironment();
//	
//	Population *p = new HomogenousPopulation(e,outputDirectory,cellCount);
//	p->equilibriate(equilibriateTime);
//	p->simulate(simTime,popDt,cellDt);
//	exit(0);
	
	
	cout<<"starting NFsim V6..."<<endl<<endl;
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
					if(strncmp(argv[2],"TLBR",4)==0)
						NFtest_TLBR::run(argc, argv);
					else if(strncmp(argv[2],"testCompare",11)==0)
						NFtest_compare::run();
					else if(strcmp(argv[2],"transformation")==0)
						NFtest_transformations::run();
					else if(strcmp(argv[2],"simple_system")==0)
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
							NG::run(argc, argv);
					}
					else
						cout<<"Could not identify run: "<<argv[2]<<endl;
				}
				else
					cout<<"You must specify which model to run."<<endl;	
			}
			
			
			else if(strcmp(argv[1],"-swim")==0)	
			{
				CapillaryEnvironment *e = new CapillaryEnvironment();
				
				
				double speed = 20.0; // um/s
				double rotDifConstant = 0.06205;
				double cellInternalStartTime = 0.0;
				
				NGparam p;
				p.setFullLite();
				p.setLiteSystem();
				
				Cell *c = new Cell(e, speed, rotDifConstant, cellInternalStartTime, "",p);
				for(int i=0; i<2000; i++)
				{
					c->swimToNewPosition(0.1);
					cout<<c->getXposition()<<","<<c->getYposition()<<","<<c->getZposition()<<";\n";
				}
				
				//e->checkPos();
				delete e;
			}
			
			else if(strcmp(argv[1],"-capillary")==0)	
			{
				//CapillaryEnvironment *e = new CapillaryEnvironment();
				//e->checkPos();
				
				
				if(argc>2)
				{
					const char * outputDirectory = argv[2];
					
//					char *p; int c = strtoi(argv[3],&p);
//					if(*p != '\n' && *p != '\0')
//					{
//						cout<<"invalid number of cells given..."<<endl; exit(1);
//					}
					
					
					int cellCount=10;
					
					double equilibriateTime = 300;
					double simTime = 2400;
					double popDt = 5;
					double cellDt = 0.01;
					
					Environment *e = new CapillaryEnvironment();
					Population *pop;
					if(argc>3)
					{
						const char * cheRlevel = argv[3];
						pop = new HomogenousPopulation(e,outputDirectory,cellCount, cheRlevel);
					}
					else
						pop = new HomogenousPopulation(e,outputDirectory,cellCount);
					
					
					cout<<"Population created with "<<cellCount<<" cells."<<endl;
					pop->equilibriate(equilibriateTime);
					pop->simulate(simTime,popDt,cellDt);
					
					
					delete pop;
					delete e;
				} else {
					cout<<"Specify the filename for the capillary simulation and the number of cells"<<endl;
				}
			}
			
			else if(strcmp(argv[1],"-xml")==0)	
			{
				if(argc>2)
				{
					System *s = NFinput::initializeFromXML(argv[2]);
					delete s;
				}
				else
					cout<<"You must specify an xml file to read."<<endl;	
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



