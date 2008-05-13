


#include "simple_system.hh"


using namespace NFtest_simple_system;
using namespace NFcore;



void NFtest_simple_system::run()
{
	cout<<"Running the simple system"<<endl;

	/**
	 * 
	 * 
	 * As a simple demonstration of how to create and run a simulation, here is an example
	 * where we have a system that looks like:
	 * 
	 * X(p~1) -> X(p~0)
	 * X(y,p~0) + Y(x) <-> X(y!1,p~0).Y(x!1)
	 * X(y!1,p~0).Y(x!1) -> X(y,p~1) + Y(x)
	 * 
	 * This is basically a simple enzymatic reaction where Y is the enzyme which can
	 * phosphorylate X, and X can auto-dephosphorylate.  This system has basic binding,
	 * unbinding, and state change reactions.  Below are the numbered steps of getting
	 * this system together.  Notice the functions defined at the end that do the actual work.
	 * 
	 * To run this example, call the NFsim program as follows (either from the Release or Debug
	 * directories, or using eclipse's launch):
	 * 
	 *        ./NFsim6 -test simple_system
	 * 
	 */
	
	
	//  1)  Create the reaction system by creating an object called System with the given name
	//      This will be the base unit for running simulations and outputting results
	System *s = new System("Simple System");
	
	
	//  2)  Create the types of molecules that are in the system (see functions below)
	//      MoleculeTypes contain all the information about the types of molecules that can exist.
//	MoleculeType *molX = createX(s);
//	MoleculeType *molY = createY(s);
	
	
	//  3)  Instantiate the actual molecules (this populate function is the easiest way, but you can do it
	//      manually as well by creating each molecule separately - see the populate function for details
	//      on how this can be done).
//	molY->populateWithDefaultMolecules(500);
//	molX->populateWithDefaultMolecules(5000);
	
	
	//  4)  Create the reactions and add them to the system.  These are calls to specific functions
	//      below where I set up the details of the reactions.  The numbers are in rates and are in
	//      arbitrary units for now.
//	ReactionClass * x_dephos = createReactionXDephos(molX, 0.4);
//	ReactionClass *rXbindY = createReactionXYbind(molX, molY, 10.0);
//	ReactionClass *rXunbindY = createReactionXYunbind(molX, molY, 5.0);
//	ReactionClass *rYphosX = createReactionYphosX(molX, molY, 0.5);
	
//	s->addReaction(x_dephos);
//	s->addReaction(rXbindY);
//	s->addReaction(rXunbindY);
//	s->addReaction(rYphosX);
	
	
	//  5)  Add the observables that we want to track throughout the simulation.  Again, to 
	//      see how this is done, see the function below.
//	addObs(s, molX, molY);
	
	
	
	//  6)  Prepare the system for simulation (this adds molecules to reactionLists
	//      and counts up the observables)
	s->prepareForSimulation();
	s->printAllReactions();
	
	
	//  7)  Register the output file name (This will put the file in your working directory, which
	//      if you run from eclipse is in the base NFsim directory) )
	//      Here, you also want to output the headers to the file, which is not done automatically
	//      because you can run a simulation with multiple calls to the sim functions.
	s->registerOutputFileLocation("simple_system_output.txt");
	s->outputAllObservableNames();
	
	
	
	
	//8)  Run the simulation!
	
	//You can optionally equilibriate the system for a set amount of time where output is not
	//recorded (all times are in seconds) using this function where the first parameter is the
	//length of the equilibriation, and the second is the number of times we want to print a
	//message that says how we are doing.  After it equilibriates, the simulation time is reset to zero.
	//s->equilibriate(50,10);
	
	//There are two ways to run a simulation.  First, you can just call the function sim as in:
	s->sim(500,500);
	
	//Calling this sim function is the easist way to run a simulation.  The first parameter is the
	//number of seconds you want to run for, the second is the number of times you want to output
	//the results to the outputfile.  So in this example, we run for 500 seconds and output 500 times,
	//so we output once per second.
	
	
	//The second way to run a simulation is to call this stepTo function:
	cout<<"Calling the stepTo function and stepping to 600 seconds"<<endl;
	double stoppingTime = s->stepTo(600);
	cout<<"The last reaction was fired at simulation time: "<<stoppingTime<<endl;
	
	//This function runs the simulation until the given time is reached, and it also returns the
	//time of the simulation when the last reaction fired.  This function does not output to a file
	//automatically, so to output using this function, you will have to use a call to this function:
	s->outputAllObservableCounts();
	
	//The stepTo function will have to be written as part of a loop where you decide when it should
	//output.  This gives you more control, but of course requires more work.
	s->printAllReactions();
	
	
	
	//Be nice and clean up when we are done.  Deleting the system (should) delete everything else
	//associated with the system, so we don't ever have to remember to delete individual reactions, molecules
	//moleculeTypes, etc...
	delete s;
}






/*MoleculeType * NFtest_simple_system::createX(System *s)
//{
	// create MoleculeType X  with one binding site and one state 
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "y";
	
	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = "p";
	
	//This is the default state value that new molecules are created with
	int * stateValues = new int [numOfStates];
	stateValues[0] = 0;
	
	//When we create a molecule, it automatically adds itself to the system, 
	MoleculeType *molX = new MoleculeType("MolX",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	return molX;
}

MoleculeType * NFtest_simple_system::createY(System *s)
{
	// create MoleculeType Y 
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "x";
	
	int  numOfStates = 0;
	string * stateNames = new string [numOfStates];
	int * stateValues = new int [numOfStates];
	
	MoleculeType *molY = new MoleculeType("MolY",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	return molY;
}

*/



//ReactionClass * NFtest_simple_system::createReactionXDephos(MoleculeType *molX, double rate)
//{
	
/*	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addStateValue("p",1);
	
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp );
	
	ReactionClass *r = new ReactionClass("X_dephos",templates, rate);
	NFcore::Transformation::genStateChangeTransform(xTemp,"p",0,r);
	return r;*/
//}

//ReactionClass * NFtest_simple_system::createReactionYphosX(MoleculeType *molX, MoleculeType *molY, double rate)
//{
/*	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addStateValue("p",0);
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	TemplateMolecule::bind(xTemp,"y",yTemp,"x");
	
	vector <TemplateMolecule *> templates;
	templates.push_back( yTemp );
	ReactionClass *r = new ReactionClass("Y_phos_X",templates,rate);
	NFcore::Transformation::genUnbindingTransform(yTemp,"x",r);
	NFcore::Transformation::genStateChangeTransform(xTemp,"p",1,r);
	return r;*/
//}


//ReactionClass * NFtest_simple_system::createReactionXYbind(MoleculeType *molX,MoleculeType *molY, double rate)
//{
/*	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addEmptyBindingSite("y");
	xTemp->addStateValue("p",0);
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	yTemp->addEmptyBindingSite("x");
		
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp );
	templates.push_back( yTemp );
	ReactionClass *r = new ReactionClass("Y_bind_X",templates,rate);
	NFcore::Transformation::genBindingTransform(yTemp,xTemp,"x","y",r);
	return r;*/
//}

//ReactionClass * NFtest_simple_system::createReactionXYunbind(MoleculeType *molX, MoleculeType *molY, double rate)
//{
/*	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	TemplateMolecule::bind(xTemp,"y",yTemp,"x");
		
	vector <TemplateMolecule *> templates;
	templates.push_back( yTemp );
	ReactionClass *r = new ReactionClass("Y_unbind_X",templates,rate);
	NFcore::Transformation::genUnbindingTransform(yTemp,"x",r);
	return r;*/
//}





/*void NFtest_simple_system::addObs(System * s, MoleculeType *molX, MoleculeType *molY)
{
	//To create an observable, we must first create a TemplateMolecule that we would
	//like to match.  For instance, to create an observable for X, for create a template
	//molecule like this:
	TemplateMolecule *xNotPhos = new TemplateMolecule(molX);
	
	//Then, we would like to set some constraints.  For this, let us set the constraint
	//that X has an open binding site at 'y' and it is not phosphorylated.
	xNotPhos->addStateValue("p",0);
	xNotPhos->addEmptyBindingSite("y");
	
	
	//Now, we create an observable from the templateMolecule and give it a name
	//that will be used in the output.
	Observable * obsxNotPhos = new Observable("X(p~0)_free",xNotPhos);
	
	
	//Finally, we have to add the observable to the MoleculeType that is being observed.  If you
	//don't add to the correct moleculeType, the count will always be zero!  And that's it!  Adding
	//to the moleculeType will record the observable with the system, and this observable will
	//always output correctly.
	molX->addObservable(obsxNotPhos);
	
	
	
	//Below I do the same for some other species...
	
	//Xp
	TemplateMolecule *xPhos = new TemplateMolecule(molX);
	xPhos->addStateValue("p",1);
	xPhos->addEmptyBindingSite("y");
	Observable * obsxPhos = new Observable("X(p~1)_free",xPhos);
	molX->addObservable(obsxPhos);
	
	
	//X-Y
	TemplateMolecule *xBoundP = new TemplateMolecule(molX);
	TemplateMolecule *yBound2 = new TemplateMolecule(molY);
	TemplateMolecule::bind(xBoundP,"y",yBound2,"x");
	Observable * obsXBoundP = new Observable("XY",xBoundP);
	molX->addObservable(obsXBoundP);
	
	//Yfree
	TemplateMolecule *yFree = new TemplateMolecule(molY);
	yFree->addEmptyBindingSite("x");
	Observable * obsyFree = new Observable("Y_free",yFree);
	molY->addObservable(obsyFree);
	
	
	TemplateMolecule *xTot = new TemplateMolecule(molX);
	Observable * obsxTot = new Observable("Xtot",xTot);
	molX->addObservable(obsxTot);
	
	TemplateMolecule *yTot = new TemplateMolecule(molY);
	Observable * obsyTot = new Observable("Ytot",yTot);
	molY->addObservable(obsyTot);
	
}
*/








