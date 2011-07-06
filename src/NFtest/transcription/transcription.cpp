


//All we need in this file is to include the header file
//that has the function prototypes, the other includes, and
//declaration that we are using the NFcore namespace.
#include "transcription.hh"




void NFtest_transcription::run()
{
	cout<<"Running the transcription system"<<endl;




	//  1)
	System *s = new System("Transcription System");


	//  2)
	MoleculeType *molRNA = createRNA(s);


	//  3)  Instantiate the actual molecules (this populate function is the easiest way, but you can do it
	//      manually as well by creating each molecule separately - see the MoleculeType::populate function for
	//      an example on how this can be done).
	molRNA->populateWithDefaultMolecules(500);


	//  4)  Create the reactions and add them to the system.  These are calls to specific functions
	//      below where I set up the details of the reactions.  The numbers are the rates and are in
	//      arbitrary units here.  In general, the rates should be in units of per second.
	ReactionClass * rna_degrade = createReactionRNAdegrades(molRNA, 0.5);
	ReactionClass *rna_transcribe = createReactionRNAtranscribed(molRNA, 100.0);

	s->addReaction(rna_degrade);
	s->addReaction(rna_transcribe);


	//  5)  Add the observables that we want to track throughout the simulation.  Again, to
	//      see how this is done, see the function below.
	addObs(s, molRNA);



	//  6)  Prepare the system for simulation (this adds molecules to reactionLists
	//      and counts up the observables)
	s->prepareForSimulation();
	s->printAllReactions();


	//  7)  Register the output file name (This will put the file in your working directory)
	//      Here, you also want to output the header to the file, which is not done automatically
	//      because you can run a simulation with multiple calls to the sim functions.
	s->registerOutputFileLocation("transcription_system_output.txt");
	s->outputAllObservableNames();




	//8)  Run the simulation!

	//You can optionally equilibrate the system for a set amount of time where output is not
	//recorded (all times are in seconds) using this function where the first parameter is the
	//length of the equilibration, and the second is the number of times we want to print a
	//message that says how we are doing.  After it equilibrates, the simulation time is reset to zero.
	//The first parameter is the number of seconds we want to equilibrate.  The second (optional) parameter
	//is the number of times you want to print an 'ok' output message.  If the second parameter is not
	//given, nothing is outputted to the console.
	s->equilibrate(0,10);

	//There are two ways to run a simulation.  First, you can just call the function sim as in:
	s->sim(50,50);

	//Calling this sim function is the easist way to run a simulation.  The first parameter is the
	//number of seconds you want to run for, the second is the number of times you want to output
	//the results to the outputfile.  So in this example, we run for 500 seconds and output 500 times,
	//so we output once per second.


	//The second way to run a simulation is to call this stepTo function:
	cout<<endl<<endl<<"Calling the stepTo function and stepping to the system time t=600 seconds"<<endl;
	double stoppingTime = s->stepTo(600);
	cout<<"The last reaction was fired at simulation time: "<<stoppingTime<<endl<<endl;

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








MoleculeType * NFtest_transcription::createRNA(System *s)
{
	// create MoleculeType X  with one binding site and one state.  In NFsim, you have to
	//specify arrays that name the binding sites, the states, and give a default state value.
	//The default binding site value, of course, is unbound.
	int numOfBsites = 1;
	string * bSiteNames = new string [numOfBsites];
	bSiteNames[0] = "ribosome";

	int numOfStates = 1;
	string * stateNames = new string [numOfStates];
	stateNames[0] = "activated";

	//This is the default state value that new molecules are created with
	int * stateValues = new int [numOfStates];
	stateValues[0] = 0;


	vector <string> compName;
	vector <string> defaultCompState;
	vector <vector <string> > possibleCompStates;

	compName.push_back("ribosome");
	defaultCompState.push_back("");
	vector <string> p1;
	possibleCompStates.push_back(p1);


	compName.push_back("activiated");
	defaultCompState.push_back("inactive");
	vector <string> p2;
	p2.push_back("inactive");
	p2.push_back("active");
	possibleCompStates.push_back(p2);



	MoleculeType *molRNA = new MoleculeType("RNA",compName,defaultCompState,possibleCompStates,s);


	//When we create a molecule, it automatically adds itself to the system, so all we have to
	//do here is create it, and return it so we can use it to add reactions to the system
	//MoleculeType *molRNA = new MoleculeType("RNA",stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	return molRNA;
}











ReactionClass * NFtest_transcription::createReactionRNAdegrades(MoleculeType *molRNA, double rate)
{
	//Here is your first glimpse at defining a reaction rule.  This is the simplist rule
	//as it only has a single state change operation.  There are several straightforward steps.
	//First, you have to define the set of TemplateMolecules that represent the possible reactants.
	//Here, we create a templateMolecule representing molecules of type X, and set the state to
	//be phosphorylated.
	TemplateMolecule *rnaTemp = new TemplateMolecule(molRNA);
	//xTemp->addStateValue("activated",1);


	//We have to create a vector (basically a storage array) for the Template Molecules that we
	//want to add to our reaction.  We do this using the standard library class std::vector.  We
	//only have one reactant in this reaction, so we just put it on the vector.
	vector <TemplateMolecule *> templates;
	templates.push_back( rnaTemp );

	TransformationSet *ts = new TransformationSet(templates);
	ts->addDeleteMolecule(rnaTemp,0);
	ts->finalize();

	//Now we can create our reaction.  This is simple: just give it a name, a rate, and the transformation
	//set that you just created.  It will take care of the rest!
	ReactionClass *r = new BasicRxnClass("RNA_degradation",rate,"",ts,molRNA->getSystem());
	return r;
}


ReactionClass * NFtest_transcription::createReactionRNAtranscribed(MoleculeType *molRNA, double rate)
{
	vector <TemplateMolecule *> templates;
	vector <TemplateMolecule *> newProduct;

	TemplateMolecule *newRnaTemp = new TemplateMolecule(molRNA);
	newProduct.push_back(newRnaTemp);



	TransformationSet *ts = new TransformationSet(templates);
	SpeciesCreator *sc = new SpeciesCreator(newProduct);
	ts->addAddSpecies(sc);
	ts->finalize();

	ReactionClass *r = new BasicRxnClass("RNA_transcription",rate,"",ts,molRNA->getSystem());
	return r;
}





void NFtest_transcription::addObs(System * s, MoleculeType *molRNA)
{
	//To create an observable, we must first create a TemplateMolecule that we would
	//like to match.  For instance, to create an observable for X, first create a template
	//molecule like this:
	TemplateMolecule *totalRNA = new TemplateMolecule(molRNA);



	//Now, we create an observable from the templateMolecule and give it a name
	//that will be used in the output.
	MoleculesObservable * obsTotalRNA = new MoleculesObservable("RNA",totalRNA);


	//Finally, we have to add the observable to the MoleculeType that is being observed.  If you
	//don't add to the correct moleculeType, the count will always be zero!  And that's it!  Adding
	//to the moleculeType will record the observable with the system, and this observable will
	//always output correctly.
	molRNA->addMolObs(obsTotalRNA);
}








