#ifndef NG_POPULATION_HH_
#define NG_POPULATION_HH_


#include "../environment/environment.hh"
#include "../cell/cell.hh"


namespace NG_swim {

class Cell;
class Environment;

class Population {
	
	public:
		Population(Environment *e);
		virtual ~Population();
		
		void equilibriate(double duration);
		void simulate(double duration, double singleCellStep, double dt);
		void addCell(Cell *c);
	
	protected:
		vector <Cell *> cells;
		Environment *env;
		vector<Cell *>::iterator cellIter;
		double globalTime;
};



class HomogenousPopulation:public Population {
	
	public:
		HomogenousPopulation(Environment *e, const char * outputDirectory, int cellCount);
		HomogenousPopulation(Environment *e, const char * outputDirectory, int cellCount, const char * cheRlevel);
		virtual ~HomogenousPopulation();
		
		
	
};



}

#endif /*NG_POPULATION_HH_*/
