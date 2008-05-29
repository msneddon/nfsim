#include "population.hh"

#include <string>





using namespace std;
using namespace NG_swim;



Population::Population(Environment *e)
{
	this->env = e;
	this->globalTime = 0;
}

Population::~Population()
{
	Cell *c;
	while(cells.size()>0)
	{
		c = cells.back();
		cells.pop_back();
		delete c;
	}
}
		

void Population::addCell(Cell *c)
{
	cells.push_back(c);
}


void Population::equilibriate(double duration)
{
  	for( cellIter = cells.begin(); cellIter != cells.end(); cellIter++ )
  	{
  		(*cellIter)->equilibriate(duration);
  	}
  	
}


void Population::simulate(double duration, double singleCellStep, double dt)
{
	clock_t start,stop;
	double elapsedTime=0;
	start = clock();
	
	while(globalTime+singleCellStep<=duration)
	{
		stop = clock();
		elapsedTime += (double(stop)-double(start))/CLOCKS_PER_SEC;
		start = stop;
		cout<<"Real Time:\t"<<elapsedTime;
		
		cout<<"\tAdvancing pop to sim time:\t"<<globalTime+singleCellStep<<"\t(";
		for( cellIter = cells.begin(); cellIter != cells.end(); cellIter++ )
		{
			cout<<"-"; cout.flush();
			(*cellIter)->stepTo(globalTime+singleCellStep,dt);
		}
		globalTime+=singleCellStep;
		cout<<")"<<endl;
	}
}



HomogenousPopulation::HomogenousPopulation(Environment *e, const char * outputDirectory, int cellCount) : Population(e)
{
	//Default parameters
	double speed = 20.0; // um/s
	double rotDifConstant = 0.06205;
	double cellInternalStartTime = 0.0;
	
	NGparam p;
	p.setFullLite();
	p.setLiteSystem();
	
	for(int i=0; i<cellCount; i++)
	{
		Cell *c = new Cell(e, speed, rotDifConstant, cellInternalStartTime, outputDirectory,p);
		this->addCell(c);
	}
}

HomogenousPopulation::HomogenousPopulation(Environment *e, const char * outputDirectory, int cellCount, const char * cheRlevel) : Population(e)
{
	//Default parameters
	double speed = 20.0; // um/s
	double rotDifConstant = 0.06205;
	double cellInternalStartTime = 0.0;
	
	NGparam p;
	p.setFullLite();
	p.setLiteSystem();
	cout<<"here"<<endl;
	cout<<"Cher Level:  '"<<cheRlevel<<"'"<<endl;
	if(strcmp(cheRlevel,"r_2x")==0) p.setCheR_2x();
	else if(strcmp(cheRlevel,"r_4x")==0) p.setCheR_4x();
	else if(strcmp(cheRlevel,"r_8x")==0) p.setCheR_8x();
	
	for(int i=0; i<cellCount; i++)
	{
		Cell *c = new Cell(e, speed, rotDifConstant, cellInternalStartTime, outputDirectory,p);
		this->addCell(c);
	}
}


HomogenousPopulation::~HomogenousPopulation()
{
}
	
	