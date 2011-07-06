


//All we need in this file is to include the header file
//that has the function prototypes, the other includes, and
//declaration that we are using the NFcore namespace.
#include "simple_system.hh"




void NFtest_ss::run()
{
	cout<<"Running the simple system"<<endl;

	/**
	 *
	 * This example is a good starting point to learn the basics of the underlying NFsim
	 * code and what happens behind the scenes.  Generally, you will specify the system
	 * through BioNetGen and NFsim will just parse the resulting xml model specification
	 * file.  However, for learning the code, there is nothing better than hardcoding a
	 * simple system yourself and learning which functions get called when.
	 *
	 * So, as a simple demonstration of how to create and run a simulation, here is an example
	 * where we have a system that looks like this:
	 *
	 * X(p~1) -> X(p~0)
	 * X(y,p~0) + Y(x) <-> X(y!1,p~0).Y(x!1)
	 * X(y!1,p~0).Y(x!1) -> X(y,p~1) + Y(x)
	 *
	 * This is basically a simple enzymatic reaction where Y is the enzyme which can
	 * phosphorylate X, and X can auto-dephosphorylate.  This system has basic binding,
	 * unbinding, and state change reactions along with a reaction that unbinds and has
	 * a state change.  Below are the numbered steps of getting this system together by
	 * hardcoding the reaction rules.  Notice the functions defined at the end that do
	 * the actual work.
	 *
	 * To run this example, call the NFsim program as follows:
	 *
	 *        ./NFsim6 -test simple_system
	 *
	 */

	//First we define some parameters for rates and counts
	int numOfMoleculeY = 3011;
	int numOfMoleculeX = 6022;
	double dephosRate = 0.2;
	double kOn = 0.0003;
	double kOff = 0.2;
	double kCat = 0.1;




	//  1)  Create the reaction system by creating an object called System with the given name
	//      This will be the base unit for running simulations and outputting results
	System *s = new System("Simple System");


	//  2)  Create the types of molecules that are in the system (see functions below)
	//      MoleculeTypes contain all the information about the types of molecules that can exist.
	MoleculeType *molX = createX(s);
	MoleculeType *molY = createY(s);


	//  3)  Instantiate the actual molecules (this populate function is the easiest way, but you can do it
	//      manually as well by creating each molecule separately - see the MoleculeType::populate function for
	//      an example on how this can be done).
	molY->populateWithDefaultMolecules(numOfMoleculeY);
	molX->populateWithDefaultMolecules(numOfMoleculeX);


	//  4)  Create the reactions and add them to the system.  These are calls to specific functions
	//      below where I set up the details of the reactions.  The numbers are the rates and are in
	//      arbitrary units here.  In general, the rates should be in units of per second.
	ReactionClass * x_dephos = createReactionXDephos(molX, dephosRate);
	ReactionClass *rXbindY = createReactionXYbind(molX, molY, kOn);
	ReactionClass *rXunbindY = createReactionXYunbind(molX, molY, kOff);
	ReactionClass *rYphosX = createReactionYphosX(molX, molY, kCat);

	s->addReaction(x_dephos);
	s->addReaction(rXbindY);
	s->addReaction(rXunbindY);
	s->addReaction(rYphosX);


	//  5)  Add the observables that we want to track throughout the simulation.  Again, to
	//      see how this is done, see the function below.
	addObs(s, molX, molY);



	//  6)  Prepare the system for simulation (this adds molecules to reactionLists
	//      and counts up the observables)
	s->prepareForSimulation();
	s->printAllReactions();


	//  7)  Register the output file name (This will put the file in your working directory)
	//      Here, you also want to output the header to the file, which is not done automatically
	//      because you can run a simulation with multiple calls to the sim functions.
	s->registerOutputFileLocation("simple_system_output.txt");
	s->outputAllObservableNames();




	//8)  Run the simulation!

	//You can optionally equilibriate the system for a set amount of time where output is not
	//recorded (all times are in seconds) using this function where the first parameter is the
	//length of the equilibriation, and the second is the number of times we want to print a
	//message that says how we are doing.  After it equilibriates, the simulation time is reset to zero.
	//The first parameter is the number of seconds we want to equilibriate.  The second (optional) parameter
	//is the number of times you want to print an 'ok' output message.  If the second parameter is not
	//given, nothing is outputted to the console.
	//s->equilibriate(50,10);

	//There are two ways to run a simulation.  First, you can just call the function sim as in:
	s->sim(500,500);

	//Calling this sim function is the easist way to run a simulation.  The first parameter is the
	//number of seconds you want to run for, the second is the number of times you want to output
	//the results to the outputfile.  So in this example, we run for 500 seconds and output 500 times,
	//so we output once per second.


	//The second way to run a simulation is to call this stepTo function:
	cout<<endl<<endl<<"Calling the stepTo function and stepping to the system time t=600 seconds"<<endl;
	//double stoppingTime = s->stepTo(600);
	//cout<<"The last reaction was fired at simulation time: "<<stoppingTime<<endl<<endl;

	//This function runs the simulation until the given time is reached, and it also returns the
	//time of the simulation when the last reaction fired (which is now the current time in the
	//system.  This function does not output to a file automatically, so to output using this function,
	//so to get results, you will have to use a call to this function:
	s->outputAllObservableCounts();

	//The stepTo function will have to be written as part of a loop where you decide when it should
	//output.  This gives you more control, but of course requires more work.


	//Here we can print out the list of reactions.  This is nice because it tells us how many times
	//a reaction fired, and how many molecules are in each list at the end of the simulation
	s->printAllReactions();



	//Be nice and clean up when we are done.  Deleting the system (should) delete everything else
	//associated with the system, so we don't ever have to remember to delete individual reactions, molecules
	//moleculeTypes, etc...
	delete s;
}








MoleculeType * NFtest_ss::createX(System *s)
{
	//In NFsim, there are several ways to declare a new MoleculeType.  First, we will
	//go through the extensive way, where we specify every option for a molecule of
	//type X.  This includes the specification of all components of X (components
	//can bind to other molecules and can have some state attribute) in addition
	//to default state values and the complete list of all the possible state
	//values for each state.

	//First, we declare three vectors which we will use as we create the molecule
	//type.  they are:
	vector <string> compName;  //Vector of the names of each component
	vector <string> defaultCompState;  //Vector of the default state of each component
	vector < vector <string> > possibleCompStates;  //Vector to hold the possible states for each component

	//First, let's create the binding site for molecule Y.
	compName.push_back("y");  //This declares the name of Y
	defaultCompState.push_back("No State");  //This is arbitrary - there is no state associated with y
	vector <string> possibleYstates; //declare a vector that holds possible states for y (which is empty)
	possibleCompStates.push_back(possibleYstates); //add the possible states for y to the vector
	                                               //we will pass to the constructor

	//Next, let's create component p, a phosphorylation site
	compName.push_back("p");  //again, here is the name of the component
	defaultCompState.push_back("Unphos");  //It begins in the Unphos state
	vector <string> possibleXstates;  //Here is a vector now of all the possible states for component p
	possibleXstates.push_back("Unphos");  //It can be Unphos
	possibleXstates.push_back("Phos");  //Or it can also be Phos.
	possibleCompStates.push_back(possibleXstates);  //Remember the possible states

	//Now, we create the MoleculeType object
	MoleculeType *molX = new MoleculeType("MolX", compName, defaultCompState, possibleCompStates, s);
	return molX;
}

MoleculeType * NFtest_ss::createY(System *s)
{
	vector <string> compName;
	vector <string> defaultCompState;
	vector < vector <string> > possibleCompStates;


	compName.push_back("x");
	defaultCompState.push_back("No State");
	vector <string> possibleXstates;
	possibleCompStates.push_back(possibleXstates);

	//Now, we create the MoleculeType object
	MoleculeType *molY = new MoleculeType("MolY", compName, defaultCompState, possibleCompStates, s);
	return molY;
}










ReactionClass * NFtest_ss::createReactionXDephos(MoleculeType *molX, double rate)
{
	//Here is your first glimpse at defining a reaction rule.  This is the simplist rule
	//as it only has a single state change operation.  There are several straightforward steps.
	//First, you have to define the set of TemplateMolecules that represent the possible reactants.
	//Here, we create a templateMolecule representing molecules of type X, and set the state to
	//be phosphorylated.
	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addComponentConstraint("p","Phos");


	//We have to create a vector (basically a storage array) for the Template Molecules that we
	//want to add to our reaction.  We do this using the standard library class std::vector.  We
	//only have one reactant in this reaction, so we just put it on the vector.
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp );


	//Once we have our set of templateMolecules defined, we can specify the transformations on
	//those molecules. To do this, we create a TransformationSet object passing in the templates
	//Next, we use the TransformationSet functionality to specify the transform.  Here we specify
	//that the transformation should apply to the template molecule (defined above) and it should
	//change the state of "p" to a value of zero (which means unphosphorylated in our system).
	//Finally, once we have added all operations we want (there is no limit!), we have to finalize
	//our transformationSet.
	TransformationSet *ts = new TransformationSet(templates);
	ts->addStateChangeTransform(xTemp,"p","Unphos");
	ts->finalize();


	//Now we can create our reaction.  This is simple: just give it a name, a rate, the parameter name
	//of the rate (which is the empty string here, because we don't have any parameters in the system
	//declared) and the transformationset that you just created.  It will take care of the rest!
	ReactionClass *r = new BasicRxnClass("X_dephos",rate,"",ts,molX->getSystem());
	return r;
}


ReactionClass * NFtest_ss::createReactionXYbind(MoleculeType *molX,MoleculeType *molY, double rate)
{
	//Now we want to create a binding reaction.  This is a little trickier, but still not too bad.
	//First we need to create TemplateMolecules to represent our two reactants, Molecule X and
	//Molecule Y.  Remember to specify that Mol X has an empty binding site for y and that
	//Y has an empty binding site for X!  This is not done for you and can lead to errors if
	//you try to bind a site that is already bound!  Also, we want to bind only if X is not
	//phosphorylated.
	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addEmptyComponent("y");
	xTemp->addComponentConstraint("p","Unphos");
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	yTemp->addEmptyComponent("x");


	//Again, we create the vector of templates to store our reactants.  There are two reactants
	//involved, so we have to add both of them to our templates vector.  Remember that only the
	//primary reactants should be added to this vector!  If we had another molecule, say Z, that
	//is bound to Y, we would have to choose to add either Y or Z to the second reactant, not both!
	//In general, there will only be one or two templates added to this vector (for unimolecular and
	//bimolecular reactions) no matter how many templates are actually involved!  Reactions are still
	//able to access all of them, however, by traversing the bonds of these two "head" molecules.
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp );
	templates.push_back( yTemp );

	//Again, create a TransformationSet with the templates, add the binding operation between X and
	//Y by specifying the sites that are binding, (ordering between x and y do not matter in this
	//function) and finalize the TransformationSet.
	TransformationSet *ts = new TransformationSet(templates);
	ts->addBindingTransform(xTemp,"y", yTemp, "x");
	ts->finalize();

	//Create and return the reaction!
	ReactionClass *r = new BasicRxnClass("Y_bind_X",rate,"",ts,molX->getSystem());
	return r;
}

ReactionClass * NFtest_ss::createReactionXYunbind(MoleculeType *molX, MoleculeType *molY, double rate)
{
	//Now its time for the unbinding reaction.  For a reactant to match the unbinding reaction,
	//the molecules X and Y must be bound.  To do this, again create the two template molecules for
	//X and Y, and now we must specify that they are bound through the named binding sites.  We
	//do that with the bind function in the TemplateMolecule class.
	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	TemplateMolecule::bind(xTemp,"y","",yTemp,"x","");

	//Like before, we create the vector of templates.  Notice that this is a unimolecular reaction!
	//even though there are two templates, only one "species" or reactant is involved.  Therefore, we
	//have to only insert one of the templates into this list.  It doesn't matter whether it is
	//the template for X or the template for Y.  Just pick one!
	vector <TemplateMolecule *> templates;
	templates.push_back( yTemp );

	//Again, here we go with the transformationSet.  We only have to specify one unbinding transform.
	//Either Y unbinds, or X unbinds, but not both.  Here, I arbitrarily choose molecule X to unbind
	//its binding site at site "y".  The templates take care of the fact that Molecule Y is on the other
	//end and also has to be updated when this is called.
	TransformationSet *ts = new TransformationSet(templates);
	ts->addUnbindingTransform(xTemp,"y",yTemp,"y");
	ts->finalize();

	//Create the reaction in the usual way.
	ReactionClass *r = new BasicRxnClass("Y_unbind_X",rate,"",ts,molX->getSystem());
	return r;
}


ReactionClass * NFtest_ss::createReactionYphosX(MoleculeType *molX, MoleculeType *molY, double rate)
{
	//This function demonstrates how you can create reactions with multiple transformations.  Just
	//as we did for the unbinding reaction, we create the templates and bind them.  We would also
	//like to specify that molecule X must be dephosphorylated for this reaction to fire, so add
	//that constraint too.
	TemplateMolecule *xTemp = new TemplateMolecule(molX);
	xTemp->addComponentConstraint("p",0);
	TemplateMolecule *yTemp = new TemplateMolecule(molY);
	TemplateMolecule::bind(xTemp,"y","",yTemp,"x","");

	//Again, just like in the unbinding reaction, just add one of the templates to the vector
	vector <TemplateMolecule *> templates;
	templates.push_back( xTemp );

	//Create the transformation set.  Add all the operations you want and finalize once you are done.
	TransformationSet *ts = new TransformationSet(templates);
	ts->addUnbindingTransform(xTemp,"y",yTemp,"x");
	ts->addStateChangeTransform(xTemp,"p","Phos");
	ts->finalize();

	//Return the Reaction.
	ReactionClass *r = new BasicRxnClass("Y_phos_X",rate, "",ts,molX->getSystem());
	return r;
}











void NFtest_ss::addObs(System * s, MoleculeType *molX, MoleculeType *molY)
{
	//To create an observable, we must first create a TemplateMolecule that we would
	//like to match.  For instance, to create an observable for X, first create a template
	//molecule like this:
	TemplateMolecule *xNotPhos = new TemplateMolecule(molX);

	//Then, we would like to set some constraints.  For this, let us set the constraint
	//that X has an open binding site at 'y' and it is not phosphorylated.
	xNotPhos->addComponentConstraint("p",0);
	xNotPhos->addEmptyComponent("y");


	//Now, we create an observable from the templateMolecule and give it a name
	//that will be used in the output.
	Observable * obsxNotPhos = new MoleculesObservable("X(p~0,y)",xNotPhos);


	//Finally, we have to add the observable to the MoleculeType that is being observed.  If you
	//don't add to the correct moleculeType, the count will always be zero!  And that's it!  Adding
	//to the moleculeType will record the observable with the system, and this observable will
	//always output correctly.
	s->addObservableForOutput(obsxNotPhos);

}








