/*! \file NFcore.hh
    \brief Contains declarations of core NFsim classes
*/
#ifndef NFCORE_HH_
#define NFCORE_HH_

//Include stl IO and string functionality
#include <iostream>
#include <fstream>
#include <string>

//Include stl containers
#include <vector>
#include <list>
#include <queue>
#include <map>

// Include various NFsim classes from other files
#include "../NFscheduler/NFstream.h"
#include "../NFutil/NFutil.hh"
#include "../NFreactions/NFreactions.hh"
#include "moleculeLists/moleculeList.hh"
#include "../NFfunction/NFfunction.hh"
#include "../NFoutput/NFoutput.hh"
#include "reactionSelector/reactionSelector.hh"


#include "templateMolecule.hh"
#include "observable.hh"

#define DEBUG 0   			// Set to 1 to display all debug messages
#define BASIC_MESSAGE 0		// Set to 1 to display basic messages (eg runtime)

// debug nauty?
#define DEBUG_NAUTY 0

using namespace std;

//!  Contains the primary classes and functions of NFsim
/*!
  This namespace contains all the basic components needed for a typical simulation
  in NFsim.  The main classes are  defined in the NFcore.hh header file for convenience, but
  the actual definitions are in the appropriately named cpp files.  Additional classes in
  the NFcore namespace for handling reactions, functions, modules, and other features are
  spread throughout the code.  When developing NFsim, it is advisable to work in this namespace
  if you are working on core NFsim functionality.  If, instead, you are merely putting in
  an add on or manipulating a System, you should create a new namespace.
  @author Michael Sneddon
*/
namespace NFcore
{

	//Forward declarations to deal with cyclic dependencies
	class MapGenerator;
	class MappingSet;
	class ReactantList;
	class TransformationSet;
	class MoleculeList;

	class GlobalFunction;
	class CompositeFunction;
	//class FunctionReference;
	class LocalFunction;

	class Outputter;
	class DumpMoleculeType;
	class DumpSystem;

	class TemplateMolecule;
	class Observable;
	class MoleculesObservable;
	class SpeciesObservable;

	/*****************************************
	 * Class declarations
	 *    Here is the list of classes that are defined in this
	 *    header file and are used throughout the NFsim program
	 */

	class System;  /* the system that contains necessary information to set up,
	                  run, and output results of the simulation */

	class MoleculeType;  /* indicates the "type" of molecule */
	class Molecule;  /* the actual class that represents instances of a molecule */
	//class TemplateMolecule; /* used to determine if a molecule matches a particular pattern
	//                           think of these as regular expressions for molecule comparison */

	class ReactionClass; /* defines a reaction class, (in other words, a rxn rule) */

	class Observable;  /* object that moniters counts of things we want to keep track of */

	class Complex;  /* collection of molecules that are bonded to each
						other and whose membership dynamically can change*/

	class ComplexList;  /* a container class to organize complexes */

	class ReactantList;

	class ReactionSelector;



	//!  Container to organize all system complexes.
	/*!
	    @author Justin Hogg
	*/
	class ComplexList
	{
		public:
		    // constructor, destructor..
			ComplexList();
			~ComplexList();

			// initialize, queries, etc
			bool isUsingComplex ( ) const { return useComplex; }
            void setUseComplex ( bool _useComplex ) { useComplex = _useComplex; }
            void setSystem ( System * _sys ) { sys = _sys; }

            // core methods of this class
			int createComplex(Molecule * m);
			Complex * getComplex(int ID_complex) const { return allComplexes.at(ID_complex); };
			Complex * getNextAvailableComplex();
			void notifyThatComplexIsAvailable(int ID_complex);

			// output and printing
			void printAllComplexes();
			void purgeAndPrintAvailableComplexList(); /*< ONLY USE FOR DEBUG PURPOSES, AS THIS DELETES ALL COMPLEX BOOKKEEPING */
			void outputComplexSizes(double cSampleTime);
			void outputMoleculeTypeCountPerComplex(MoleculeType *m);
			double outputMeanCount(MoleculeType *m);
			double calculateMeanCount(MoleculeType *m);

            // a public iterator:
			// sets an iternal iterator to the beginning of the vector
			void resetComplexIter() {  complexIter_public = allComplexes.begin();  };

			// returns the next complex ptr on the vector and increments iterator.  Returns 0 if at the end of the vector.
			Complex * nextComplex() {  return (complexIter_public < allComplexes.end() ? *complexIter_public++ : 0);  };


		protected:
			vector <Complex * > allComplexes;         /*!< container of all complexes in the simulation */
			queue <int> nextAvailableComplex;         /*!< queue tells us which complexes can be used next */

			System * sys;                             /* pointer to the system which this ComplexList belongs to */
			bool useComplex;                          /* true if the system is tracking complexes */

		private:
			vector <Complex *>::iterator  complexIter;         /* to iterate over allComplexes */
			vector <Complex *>::iterator  complexIter_public;  /* to iterator over allComplexes, associated with public iterator */
	};



	//!  The main class that begins and runs simulations
	/*!
	   This class is what runs the actual simulation.  It keeps lists of
	   all MoleculeTypes that exist and all Reactions that can operate
	   on those moleculeTypes.  It also contains the function (named sim)
	   which runs the main simulation loop.  To create a new System,
	   first add all MoleculeTypes, all reactions, then call the function
	   prepareForSimulation().  Then call the sim function as many times
	   as you like.  Output is printed to the given output file and lists
	   all observables per time according to how often you want values
	   recorded.
	   @author Michael Sneddon
	 */
	class System
	{
	
		// _NETGEN_
	    // Netgen is a system wrapper for network generation and
		// needs access to protected elements of System.
		friend class Netgen;

		public:

			/*!
				 creates a system with the given name that does not
				 keep track dynamically of complex formation
			 */
			System(string name);

			/*!
				 creates a system that keeps track of complex formation if
				 the setComplex parameter is set to true
			 */
			System(string name, bool useComplex);

			/*!
				creates a system that keeps track of complex formation if
				the setComplex parameter is set to true
			*/
			System(string name, bool useComplex, int globalMoleculeLimit);

			/*!
				 destroys the system and cleans up all memory associated with it
			 */
			~System();

			// Basic functions to get the properties and objects of the system
			string getName() const { return name; };
			bool isUsingComplex() { return useComplex; };   // NETGEN -- is this needed?
			bool isOutputtingBinary() { return useBinaryOutput; };
			double getCurrentTime() const { return current_time; };
			int getGlobalMoleculeLimit() const { return globalMoleculeLimit; };

			int getMolObsCount(int moleculeTypeIndex, int observableIndex) const;
			Observable * getObservableByName(string obsName);
			double getAverageGroupValue(string groupName, int valIndex);

			ReactionClass *getReaction(int rIndex) { return allReactions.at(rIndex); };

			MoleculeType * getMoleculeType(int mtIndex) { return allMoleculeTypes.at(mtIndex); };
			MoleculeType * getMoleculeTypeByName(string name);
			int getNumOfMoleculeTypes() { return allMoleculeTypes.size(); };
			Molecule * getMoleculeByUid(int uid);
		    int getNumOfMolecules();

			//Functions used when setting up the system.  Note that most of these methods
			//are called automatically when you create an object (notable exception are
			//reactions), so use these methods with caution!  You probably don't need to
			//use them!
			int addMoleculeType(MoleculeType *moleculeType);
			void addReaction(ReactionClass *reaction);
			void addNecessaryUpdateReaction(ReactionClass *reaction);
			// NETGEN
			//int createComplex(Molecule * m);

			bool addGlobalFunction(GlobalFunction *gf);
			GlobalFunction * getGlobalFunctionByName(string fName);
			bool addCompositeFunction(CompositeFunction *cf);
			CompositeFunction * getCompositeFunctionByName(string fName);
			void finalizeCompositeFunctions();

			void printAllFunctions();

			bool saveSpecies() { return saveSpecies(string(name+"_nf.species")); };
			bool saveSpecies(string filename);


			LocalFunction * getLocalFunctionByName(string fName);
			//bool addFunctionReference(FunctionReference *fr);

			/* once all elements are added, we need to prepare and initialize for simulations */
			void prepareForSimulation();


			void setUniversalTraversalLimit(int utl);


			void setEvaluateComplexScopedLocalFunctions( bool val ) { evaluateComplexScopedLocalFunctions = val; };
			bool getEvaluateComplexScopedLocalFunctions( ) const { return evaluateComplexScopedLocalFunctions; };

			/* tell the system where to ouptut results*/
			void setOutputToBinary();
			void registerOutputFileLocation(string filename);


			void setDumpOutputter(DumpSystem *ds);
			void tryToDump();

			void turnOnGlobalFuncOut() { this->outputGlobalFunctionValues=true; };
			void turnOffGlobalFuncOut() { this->outputGlobalFunctionValues=false; };

			/* once all elements are added, we need to prepare and initialize for simulations */
			void tagReaction(int rID);



			void addLocalFunction(LocalFunction *lf);
			void getLocalFunction(string funcName) const { cout<<"getLocalFunction not yet implemented."<<endl;exit(1);};


			/* functions to output simulation properties to the registered file*/
			void outputAllObservableNames();
			void outputAllObservableCounts();
			void outputAllObservableCounts(double cSampleTime);
			void outputAllObservableCounts(double cSampleTime,int eventCounter);
			int getNumOfSpeciesObs() const;
			Observable * getSpeciesObs(int index) const;

			/* functions that print out other information to the console */
			// NETGEN
			//void printAllComplexes();
			void printAllReactions();
			void printIndexAndNames();
			void printAllMoleculeTypes();

			void printAllObservableCounts();
			void printAllObservableCounts(double cSampleTime);
			void printAllObservableCounts(double cSampleTime,int eventCounter);

			// NETGEN
			//void purgeAndPrintAvailableComplexList(); /*< ONLY USE FOR DEBUG PURPOSES, AS THIS DELETES ALL COMPLEX BOOKKEEPING */
			//void outputComplexSizes(double cSampleTime);
			//void outputMoleculeTypeCountPerComplex(MoleculeType *m);
			//double outputMeanCount(MoleculeType *m);
			//double calculateMeanCount(MoleculeType *m);

			void update_A_tot(ReactionClass *r, double old_a, double new_a);




			void evaluateAllLocalFunctions();


			void addObservableForOutput(Observable *o);

			void addOutputter(Outputter *op);
			void dumpOutputters();


			/* functions needed while running the simulation */
			// NETGEN
			//Complex * getComplex(int ID_complex) const { return allComplexes.at(ID_complex); };
			//Complex * getNextAvailableComplex();
			//void notifyThatComplexIsAvailable(int ID_complex);

			double sim(double time, long int sampleTimes);

			/* run the simulation for a given length of time, output all results to the registered
			 * file.  The number of sample times is the number of times the function will output to
			 * a file divided equally throughout the elapsed time  */
			double sim(double time, long int sampleTimes, bool verbose);

			/* run the simulation up until the stopping time (but not exceding the stopping time. This
			 * will not output anything to file (so must be done manually) and returns the current time
			 * of the simulation, which will always be less than the stopping time */
			double stepTo(double stoppingTime);

			void singleStep();

			/* runs the simulation without output for the given time.  After, the clock is reset to
			 * to the start time, so it is as if no time has elapsed */
			void equilibrate(double duration);
			void equilibrate(double duration, int statusReports);


			//For moleculeTypes to do a quick lookup and see where a particular reaction
			//is mapped to.  This is an optimization....
			void registerRxnIndex(int rxnId, int rxnPos, int rxnIndex) {rxnIndexMap[rxnId][rxnPos]=rxnIndex;};
			int getRxnIndex(int rxnId, int rxnPos) const { return rxnIndexMap[rxnId][rxnPos]; };

			void turnOff_OnTheFlyObs();
			void turnOnOutputEventCounter() { outputEventCounter=true; };

			void addParameter(string name,double value);
			double getParameter(string name);
			void setParameter(string name, double value);
			void updateSystemWithNewParameters();
			void printAllParameters();

	        NFstream& getOutputFileStream();

	        // NETGEN -- method to access allComplexes
	        ComplexList & getAllComplexes( )  {  return allComplexes;  };

			/*! keeps track of null events (ie binding events that have
			    been rejected because molecules are on the same complex)
			 */
			static int NULL_EVENT_COUNTER;

			/*!
				turns on csv format, so that instead of a gdat file, a comma delimited
				file is generated.
			*/
			void turnOnCSVformat() { this->csvFormat = true; };

		protected:

			///////////////////////////////////////////////////////////////////////////
			// The invariant system properties, created when the system is created and before
			// the system is prepared
			string name;         /*!< arbitrary name of the system  */
			// NETGEN -- is this needed?
			bool useComplex;     /*!< sets whether or not to dynamically track complexes */
			bool useBinaryOutput; /*!< set to true to turn on binary output of data */
			bool evaluateComplexScopedLocalFunctions; /*!< set to true to turn on enable complex-scoped local functions */
			int universalTraversalLimit; /*!< sets depth to traverse molecules when updating reactant lists */
			bool onTheFlyObservables;    /*!< sets whether or not observables are calculated on the fly */
		    bool outputGlobalFunctionValues; /*< set to true to output the value of all global functions at each output step */
		    int globalMoleculeLimit; /*< total number of any particular molecule that can be created, default=100,000 */
		    bool outputEventCounter; /*< set to true to output the cumulative number of events at each output step */

		    int globalEventCounter;

		    ///////////////////////////////////////////////////////////////////////////
			// The container objects that maintain the core system configuration
			vector <MoleculeType *> allMoleculeTypes;  /*!< container of all MoleculeTypes in the simulation */
			vector <ReactionClass *> allReactions;    /*!< container of all Reactions in the simulation */
            // NETGEN
			//vector <Complex * > allComplexes;         /*!< container of all complexes in the simulation */
			//queue <int> nextAvailableComplex;         /*!< queue tells us which complexes can be used next */
			vector <Outputter *> allOutputters;    /*!< manages the outputters of the system */

			ComplexList  allComplexes;                /*!< a container to track all complexes in the system */

			vector <Observable *> obsToOutput; /*!< keeps ordered list of pointers to observables for output */
			vector <Observable *> speciesObservables;

			DumpSystem *ds;


			///////////////////////////////////////////////////////////////////////////
			// The container objects that maintain the functional expressions
			//vector <FunctionReference *> functionReferences;
			vector <GlobalFunction *> globalFunctions;    /*!< container of all global functions available to the system */
			vector <LocalFunction *> localFunctions;      /*!< container of all local functions available to the system */
			vector <ReactionClass *> necessaryUpdateRxns; /*!< list of all  reactions that need to update propensity after each step*/

			vector <CompositeFunction *> compositeFunctions;


			///////////////////////////////////////////////////////////////////////////
			// Properties of the system that update in time
			double a_tot;        /*< the sum of all a's (propensities) of all reactions */
			double current_time; /*< keeps track of the simulation time */
			ReactionClass * nextReaction;  /*< keeps track of the next reaction to fire */

			///////////////////////////////////////////////////////////////////////////
			// protected functions needed only by the system while running a simulation
			double get_A_tot() const { return a_tot; };
			double recompute_A_tot();
			double getNextRxn();


			///////////////////////////////////////////////////////////////////////////
			// Neccessary variables and methods for outputting
			//ofstream outputFileStream; /* the stream to a file to write out the results */
			NFstream outputFileStream; /* NFstream is a smart stream that uses ofstream or stringstream depending on whether NF_MPI is defined */
			void outputGroupDataHeader();


			void outputAllPropensities(double time, int rxnFired);
			ofstream propensityDumpStream;

			bool csvFormat;


			///////////////////////////////////////////////////////////////////////////
			//random data structures and variables used for optimization....
			int **rxnIndexMap; /*!< maps reaction index values to a reaction, used for MoleculeTypes to
			                        quickly lookup a reaction */


			map <string,double> paramMap;



			//Data structure that performs the selection of the next reaction class
			ReactionSelector * selector;


		private:
			list <Molecule *> molList;
			list <Molecule *>::iterator molListIter;

			///////////////////////////////////////////////////////////////////////////
			//Iterators that allow fast traversal of the object containers

			vector<Observable *>::iterator obsIter;
			vector<MoleculeType *>::iterator molTypeIter;  /* to iterate over allMoleculeType */
			vector <ReactionClass *>::iterator rxnIter;    /* to iterate over allReactions */
			// NETGEN -- moved to ComplexList class
			//vector <Complex *>::iterator complexIter;      /* to iterate over allComplexes */
			vector <GlobalFunction *>::iterator functionIter; /* to iterate over Global Functions */
	};



	//!  Keeps track of the types of molecules that can exist.
	/*!
      This class maintains information about a "type" of molecule in the System.
	  It keeps track of all reactions that this "type" of molecule can possibly
	  be a part of and a list of the default binding sites and states.  It also
	  keeps track and updates all Observables that pertain to this MoleculeType.
	  And perhaps most importantly, it keeps track of all Molecules that exist in
	  the simulation of this molecule "type."  It has functions which allow it
	  also to populate itself with a number of default molecules, making initializations
	  of the simulation easier.
	    @author Michael Sneddon
	 */
	class MoleculeType
	{
		public:

			MoleculeType(
					string name,
					vector <string> &compName,
					System *s);

			MoleculeType(
					string name,
					vector <string> &compName,
					vector <string> &defaultCompState,
					System *s);

			MoleculeType(
					string name,
					vector <string> &compName,
					vector <string> &defaultCompState,
					vector < vector<string> > &possibleCompStates,
					System *system);

			MoleculeType(
					string name,
					vector <string> &compName,
					vector <string> &defaultCompState,
					vector < vector<string> > &possibleCompStates,
					vector <bool> isIntegerComponent,
					System *system);

			MoleculeType(
					string name,
					vector <string> &compName,
					vector <string> &defaultCompState,
					vector < vector<string> > &possibleCompStates,
					vector <bool> isIntegerComponent,
					bool pop_type,
					System *system);

			~MoleculeType();

			string getName() const { return name; };
			int getTypeID() const { return type_id; };
			System * getSystem() const { return system; };

			//Function to access component information
			int getNumOfComponents() const { return numOfComponents; };
			string getComponentName(int cIndex) const;
			void getPossibleComponentStates(int cIndex, list <string> &nameList);
			int getDefaultComponentState(int cIndex) const { return defaultCompState[cIndex]; };

			int getCompIndexFromName(string cName) const;
			string getComponentStateName(int cIndex, int cValue);
			int getStateValueFromName(int cIndex, string stateName) const;




			// set of functions that deal with equivalent (aka symmetric) components
		    // Here is how this works:
			// Every set of sites with an equivalent name is termed an "equivalency class".  For instance,
			// given a molecule such as: L(r,r,r,t,t), there are two equivalency classes, one for site 'r'
			// and one for site 't'.  Internally to NFsim, the sites are relabeled 'r1','r2','r3', etc.
			//
			//
			int getNumOfEquivalencyClasses() const { return this->n_eqComp; };
			// returns a string array with length numberOfEquivalencyClasses giving the generic component names
			string *getEquivalencyClassCompNames() const { return this->eqCompOriginalName; };
			void addEquivalentComponents(vector <vector <string> > &identicalComponents);
			bool isEquivalentComponent(string cName) const;
			bool isEquivalentComponent(int cIndex) const;
			void getEquivalencyClass(int *&components, int &n_components, string cName) const;

			// given a generic component name or specific component index, return the equivalence class number
			int getEquivalencyClassNumber(string cName) const;
			int getEquivalenceClassNumber(int cIndex) const;

			// given a component index, return the generic component name
			string getEquivalenceClassComponentNameFromComponentIndex(int cIndex) const;





			// query or set population type
			bool isPopulationType() const { return population_type; };

			bool isIntegerComponent(string cName) const;
			bool isIntegerComponent(int cIndex) const;

			//functions that handle the observables
			int getNumOfMolObs() const { return (int)molObs.size(); };
			string getMolObsName(int obsIndex) const;
			MoleculesObservable * getMolObs(int obsIndex) const { return molObs.at(obsIndex); };
			int getMolObsCount(int obsIndex) const;
			void removeFromObservables(Molecule * m);
			void addToObservables(Molecule * m);
			void outputMolObsNames(NFstream &fout);
			void outputMolObsCounts(NFstream &fout);
			void printMolObsNames();
			void printMolObsCounts();



			void addAllToObservables();


			//function to access particular molecules or reactions (these are really only
			//used when debugging or running the walker...
			Molecule * getMolecule(int ID_molecule) const;
			int getMoleculeCount() const;

			int getReactionCount() const { return reactions.size(); };
			int getRxnIndex(ReactionClass * rxn, int rxnPosition);



			//Functions to generate molecules, remove molecules at the beginning
			//or during a running simulation
			Molecule *genDefaultMolecule();

			void addMoleculeToRunningSystem(Molecule *&mol);
			void addMoleculeToRunningSystemButDontUpdate(Molecule *&mol);
			void removeMoleculeFromRunningSystem(Molecule *&m);
			void removeMoleculeFromRunningSystemButDontUpdate(Molecule *&m);
			void removeFromRxns(Molecule * m);



			//Adds the basic components that this MoleculeType needs to reference
			void addReactionClass(ReactionClass * r, int rPosition);
			void addMolObs(MoleculesObservable * mo) { molObs.push_back(mo); }; //could add check here to make sure observable is of this type
			int createComplex(Molecule *m) { return (system->getAllComplexes()).createComplex(m); };
			void addTemplateMolecule(TemplateMolecule *t);

			/* handle DOR reactions */
			//void addDORrxnClass(ReactionClass * r, int rPosition);
			//int getDORrxnCount() const { return indexOfDORrxns.size(); };
			//ReactionClass * getDORrxn(int index) const { return reactions.at(indexOfDORrxns.at(index)); };
			//int getDORreactantPosition(int index) const { return reactionPositions.at(indexOfDORrxns.at(index)); };
			//int getDORreactantIndex(int index) const { return indexOfDORrxns.at(index); };


			/* updates a molecules membership (assumes molecule is of type this) */
			void updateRxnMembership(Molecule * m);

			/* auto populate with default molecules */
			void populateWithDefaultMolecules(int moleculeCount);

			/* this method assumes all molecules in the simulation
			 * have been defined, and all reaction classes and observables
			 * have been added.  Then, this function will add those
			 * molecules to rxn lists and instantiate the counts of the observables.
			 * In general, you do not need to worry about this function because
			 * it automatically gets called by the System when you prepare the System*/
			void prepareForSimulation();


			//Debugging function that prints some useful information about this MoleculeType
			void printDetails() const;
			void printAllMolecules();


			//Functionality for local functions

			/* Type I local functions are functions that this molecule type needs to have
			 * always updated because it can react in a DOR reaction which needs the specified
			 * function.  DOR reactions identify these MoleculeTypes and add their function to
			 * the list of Type I functions for this MoleculeType.
			 *
			 * Type II local functions are those that require this MoleculeType as an observable
			 * or counter to be evaluated.  When this molecule gets updated, it has to automatically
			 * update all of its type II local functions.
			 *
			 * These functions return the index in the array indicating where these functions were
			 * added.  This allows local functions to quickly update molecules that need the local
			 * function information.
			 */
			int addLocalFunc_TypeI(LocalFunction *lf);
			int addLocalFunc_TypeII(LocalFunction *lf);
			vector <LocalFunction *> locFuncs_typeI;
			vector <LocalFunction *> locFuncs_typeII;

			int getNumOfTypeIFunctions() const {return locFuncs_typeI.size(); };
			LocalFunction *getTypeILocalFunction(int index) { return locFuncs_typeI.at(index); };
			int getNumOfTypeIIFunctions() const {return locFuncs_typeII.size(); };
			LocalFunction *getTypeIILocalFunction(int index) { return locFuncs_typeII.at(index); };

			int getNumOfDORrxns() const { return indexOfDORrxns.size(); };
			ReactionClass * getDORrxn(int dorRxnIndex) const { return reactions.at(indexOfDORrxns.at(dorRxnIndex)); };
			int getDORrxnIndex(int dorRxnIndex) const { return indexOfDORrxns.at(dorRxnIndex); };
			int getDORrxnPosition(int dorRxnIndex) const { return reactionPositions.at(indexOfDORrxns.at(dorRxnIndex)); };

			void setUpLocalFunctionListForMolecules();

		protected:

			void init(
				string name,
				vector <string> &compName,
				vector <string> &defaultCompState,
				vector < vector<string> > &possibleCompStates,
				vector <bool> isIntegerComponent,
				System *system);


			//basic info
			System *system;
			string name;
			int type_id;

			//keeps track of the key information about a MoleculeType - the component
			int numOfComponents;
			string *compName;
			vector < vector < string > > possibleCompStates;
			int *defaultCompState;
			bool *isIntegerCompState;
			const bool population_type;


			//set of variables to keep track of equivalent (aka symmetric) components
			int n_eqComp;
			string *eqCompOriginalName;
			int * eqCompSizes;
			string **eqCompName;
			int **eqCompIndex;


			//Lists and vectors of everything we need to know
			MoleculeList * mList;

			vector <ReactionClass *> reactions; /* List of reactions that this type can be involved with */
			vector <int> reactionPositions;   /* the position in the reaction for this type of molecule */

			vector <int> indexOfDORrxns;


			vector <MoleculesObservable *> molObs;  /* list of things to keep track of */

			vector <TemplateMolecule *> allTemplates; /* keep track of all templates that exist of this type
															so that they are easy to delete from memory at the end */

			ReactionClass *rxn; /*used so we don't need to redeclare this at every call to updateRxnMembership */



		private:
			//Some iterators so we don't have to instantiate a new iterator every time
			vector<Molecule *>::iterator molIter;  /* to iterate over mInstances */
			vector<MoleculesObservable *>::iterator molObsIter; /* to iterate over observables */
			vector <ReactionClass *>::iterator rxnIter; /* to iterate over reactions */
	};



	//!  Each molecule in the system is represented by an instance of this.
	/*!
      The base unit of the NFsim program, this class maintains instances of
	  individual objects that exist in the simulation.  Molecules are able
	  to interact with their MoleculeType so that they also 'know' about the
	  reactions and observables they need to keep updated as they react and change
	  states.  There are also two static functions which are used to create
	  and delete bonds (bind and unbind).
	    @author Michael Sneddon
	 */
	class Molecule
	{
		public:

			/* constructors / deconstuctors */
			Molecule(MoleculeType * parentMoleculeType, int listId);
			~Molecule();

			/* basic get functions for name, type, complex, and IDs*/
			int getMolListId() const { return listId; };
			string getMoleculeTypeName() const { return parentMoleculeType->getName(); };
			MoleculeType * getMoleculeType() const { return parentMoleculeType; };
			int getUniqueID() const { return ID_unique; };
			bool isAlive() const { return isAliveInSim; };
			void setAlive(bool isAlive) { isAliveInSim = isAlive; };

			void setComplexID(int currentComplex) { this->ID_complex=currentComplex; }

			int getComplexID() const { return ID_complex; };
			Complex * getComplex() const { return (parentMoleculeType->getSystem()->getAllComplexes()).getComplex(ID_complex); };
			int getDegree();

			// get (non-unqiue) label for this molecule (cIndex==-1) or one of it's components (cIndex>=0)
			string getLabel(int cIndex) const;

			////////////////////////////////////////////////////////////////////////
			bool isPopulationType() const { return parentMoleculeType->isPopulationType(); } ;
			bool setPopulation( int count );
			int  getPopulation() const;
			bool incrementPopulation();
			bool decrementPopulation();

			///////////////////////////////////////////////////////////////////////
			int getComponentState(int cIndex) const { return component[cIndex]; };
			int getComponentIndexOfBond(int cIndex) const { return indexOfBond[cIndex]; };
			void setComponentState(int cIndex, int newValue);
			void setComponentState(string cName, int newValue);

			///////////// local function methods...
			void setLocalFunctionValue(double newValue,int localFunctionIndex);
			double getLocalFunctionValue(int localFunctionIndex);
			LocalFunction * getLocalFunction(int localFunctionIndex);
			void setUpLocalFunctionList();


			////////////////////////////////////////////////////////////////////

			/* accessor functions for checking binding sites */
			bool isBindingSiteOpen(int bIndex) const;
			bool isBindingSiteBonded(int bIndex) const;
			Molecule * getBondedMolecule(int bSiteIndex) const;
			int getBondedMoleculeBindingSiteIndex(int cIndex) const;

			int getRxnListMappingId(int rxnIndex) const { return rxnListMappingId[rxnIndex]; };
			void setRxnListMappingId(int rxnIndex, int rxnListMappingId) {
					this->rxnListMappingId[rxnIndex] = rxnListMappingId;
			};

			/* set functions for states, bonds, and complexes */
			//void setState(const char * stateName, int value);
			//void setState(int stateIndex, int value);
			void setBondTo(Molecule * m2, int bindingSiteIndex);
			void moveToNewComplex(int newComplexID) { ID_complex = newComplexID; };


			/* static functions which bind and unbind two molecules */
			static void bind(Molecule *m1, int cIndex1, Molecule *m2, int cIndex2);
			static void bind(Molecule *m1, string compName1, Molecule *m2, string compName2);
			static void unbind(Molecule *m1, int bSiteIndex);
			static void unbind(Molecule *m1, char * bSiteName);


			/* functions needed to traverse a complex and get all components
			 * which is important when we want to update reactions and complexes */
			void traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit);
			static void breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth);
			void depthFirstSearch(list <Molecule *> &members);

			/* when we are ready to begin simulations, moleculeType calls this function
			 * so that this molecule can add itself to all the necessary lists */
			void prepareForSimulation();

			/* function that tells this molecule that it changed states or bonds
			 * and it should update its reaction membership */
			void updateRxnMembership();
			void removeFromObservables();
			void addToObservables();
			//void updateDORs();

			//double getDORvalueFromGroup(string groupName, int valueIndex);



			/* DOR Functions*/
			// update local functions that may depend on this molecule,
			// i.e. this molecule's type appears in an observable that some local fcn references.
			// (this version only updateds lfcns on connected to this molecule) --Justin
			void updateTypeIIFunctions();
			// update local functions that may depend on this molecule,
			// i.e. this molecule's type appears in an observable that some local fcn references.
			// (this version updates lfcns on all product complexes) --Justin
			void updateTypeIIFunctions( vector <Complex *> & productComplexes );
			// update DOR reaction propensity that depends on this molecule (type I),
			// i.e. this molecule is an argument passed to a local function.
			void updateDORRxnValues();



			/* print functions for debugging */
			void printDetails();
			void printDetails(ostream &o);
			static void printMoleculeList(list <Molecule *> &members);

			static int getUniqueIdCount() { return uniqueIdCount; };
			static const int NOT_IN_RXN = -1;


			int isObs(int oIndex) const { return isObservable[oIndex]; };
			void setIsObs(int oIndex, int isObs) { isObservable[oIndex]=isObs; };


			/* used for traversing a molecule complex */
			bool hasVisitedMolecule;
			bool * hasVisitedBond;
			TemplateMolecule *isMatchedTo;

			/* used when reevaluating local functions */
			bool hasEvaluatedMolecule;


			static const int NOSTATE = -1;
			static const int NOBOND = 0;
			static const int NOINDEX = -1;


			//void addDependentUpdateMolecule(Molecule *m);
			//void removeDependentUpdateMolecule(Molecule *m);
			//void traverseBondedNeighborhoodForUpdate(list <Molecule *> &members, int traversalLimit);


			//Unnecessary functions now that were once used to
			//mark a molecule that is not in the simulation
			//bool isMolAlive() const { return !isDead; };
			//bool isMolDead() const { return isDead; };
			//void kill() { isDead = true; };
			//void create() { isDead = false; };


		protected:


			bool isPrepared;
			bool isAliveInSim;

			/* Set of IDs which identifies uniquely this molecule */
			int ID_complex;
			int ID_type;
			int ID_unique;
			int listId;


			static int uniqueIdCount;

			/* The type of this molecule */
			MoleculeType *parentMoleculeType;
			bool useComplex;

			/* track population properties */
			int population_count;

			/* store the states and bonds in arrays */


			///////////////////////////////////////////////////////////////////
			/* list of components */
			int *component;
			int numOfComponents;
			Molecule **bond;
			int *indexOfBond; /* gives the index of the component that is bonded to this molecule */


			//////////// keep track of local function values
			double *localFunctionValues;


			//keep track of which observables I match (and how many times I match it)
			int *isObservable;


			//Used to keep track of which reactions this molecule is in...
			int * rxnListMappingId;
			int nReactions;


			// dependent update molecule list, so that when this molecule updates,
			// it necessarily updates whatever is on this list.  This list will capture non
			// local interactions that need to be maintained that result from the connected-to syntax
			// list <Molecule *> dependentUpdateMolecules
			//list <Molecule *> dependentUpdateMolecules;


		private:

			static queue <Molecule *> q;
			static queue <int> d;
			static list <Molecule *>::iterator molIter;
			//static list <Molecule *>::iterator molIter2;

	};










	//!  Abstract Base Class that defines the interface for all reaction rules.
	/*!
	    A ReactionClass represents the set of reactions implied by a single reaction
	    rule along with the rate law and propensity of that reaction rule.  ReactionClasses
	    also store information needed to transform reactants.  A ReactionClass keeps a
	    reference to every reactant Molecule in the system that might participate in the
	    rule.  To determine if a particular set of connected Molecules is a possible reactant,
	    ReactionClasses use TemplateMolecules.  ReactionClasses come in several forms based
	    on how reactants are chosen and the rate law.  For typical ReactionClasses, each
	    reactant Molecule has an equal probability of reacting and therefore reactants
	    are selected to participate at random.  However, for ReactionClasses with rates
	    that depend functionally on local context, each reactant will have a different
	    probability and rate of reacting.  Such ReactionClasses have a distribution of rates
	    and so we refer to them as Distribution of Rates (DOR) Reactions.  DOR Reactions require
	    special handling because the participating reactants must be chosen based on its
	    weighted probability of reacting.  There are also reactionClasses that can support
	    functionally defined rate laws and michaelis-menton style reactions.  All of these
	    implementing classes are declared in the file reactions.hh.
	    @author Michael Sneddon
	*/
	class ReactionClass
	{
		// _NETGEN_
		friend class MatchSetIter;
		friend class Netgen;

		public:
			static const int NO_LIMIT = -3;

			static const int BASIC_RXN = 0;
			static const int DOR_RXN = 1;
			static const int OBS_DEPENDENT_RXN = 2;
			static const int POP_RXN = 3;  // deprecated
			static const int DOR2_RXN = 4;



			ReactionClass(string name, double rate, string baseRateParameterName, TransformationSet *transformationSet, System *s);
			virtual ~ReactionClass();

			int getNumOfReactants() const { return n_reactants; };

			string getName() const { return name; };
			double getBaseRate() const { return baseRate; };
			int getRxnType() const { return reactionType; };

			void setBaseRate(double newBaseRate,string newBaseRateName);
			void resetBaseRateFromSystemParamter();

			void setTraversalLimit(int limit) { this->traversalLimit = limit; };

			double get_a() const { return a; };
			virtual void printDetails() const;
			void fire(double random_A_number);

			//For DOR reactions
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) = 0;
			virtual int getDORreactantPosition() const { cerr<<"Trying to get DOR reactant Position from a reaction that is not of type DOR!"<<endl;
															cerr<<"this is an internal error, and so I will quit."<<endl; exit(1); return -1; };
			virtual int getDORreactantPosition2() const { cerr<<"Trying to get DOR reactant Position2 from a reaction that is not of type DOR2!"<<endl;
															cerr<<"this is an internal error, and so I will quit."<<endl; exit(1); return -1; };


			//The main virtual functions that must be implemented in all implementing classes
			virtual void init() = 0; //called when the reaction is added to the system
			virtual void prepareForSimulation() = 0; //called once everything is added to the system
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos) = 0;
			virtual void remove(Molecule *m, unsigned int reactantPos) = 0;

			virtual double update_a() = 0;



			/* turn the tag of this guy on */
			void tag() { tagged = true; };


			virtual int getReactantCount(unsigned int reactantIndex) const = 0;
			virtual int getCorrectedReactantCount(unsigned int reactantIndex) const = 0;
			virtual void printFullDetails() const = 0;


			void setRxnId(int rxnId) { this->rxnId = rxnId; };
			int getRxnId() const { return rxnId; };


			void turnOff_OnTheFlyObs() { onTheFlyObservables=false; };


			void setTotalRateFlag(bool totalRate) { totalRateFlag = totalRate; };


			// _NETGEN_
			void set_match( vector <MappingSet *> & match_set );
			void apply( vector <Molecule *> & product_molecules );


		protected:
			virtual void pickMappingSets(double randNumber) const=0;

			int rxnId;

			/* if this reaction is tagged, it outputs a message everytime it is fired */
			bool tagged;

			string name;
			int reactionType;
			unsigned int n_reactants;
			unsigned int n_mappingsets;

			System * system;

			double baseRate;
			string baseRateParameterName;
			double a;
			unsigned int fireCounter;

			unsigned int traversalLimit;

			TemplateMolecule **reactantTemplates;
			TransformationSet * transformationSet;
			MappingSet **mappingSet;

			bool onTheFlyObservables;
			bool isDimerStyle;

			list <Molecule *> products;
			list <Molecule *>::iterator molIter;

			// remember the molecule type of each product molecule a with typeII dependencies
			list <MoleculeType *> typeII_products;
			list <MoleculeType *>::iterator typeII_iter;

			//Used by the reaction class to make sure that it only updates
			//each complex once (for observables, and matchOnce reactants)
			vector <int> updatedComplexes;
			vector <Complex*> productComplexes;  // Justin 24Jun12
			vector <Complex*>::iterator complexIter; // Justin 24Jun12


			/* flag to identify if the macroscopic vs. microscopic rate is to
			 * be used (TotalRate = macroscopic)
			 */
			bool totalRateFlag;

			/* flag population reactants */
			bool *isPopulationType;

			/* if population reactants are identical, this is the discrete
			 * count correction for calculating the ratelaw
			 */
			int *identicalPopCountCorrection;
	};













	//!  Container to dynamically keep track of all system complexes.
	/*!
	    @author Michael Sneddon
	*/
	class Complex
	{
		public:
			Complex(System * s, int ID_complex, Molecule * m);
			~Complex();


			bool isAlive();
			int getComplexID() const { return ID_complex; };
			int getComplexSize() const {return complexMembers.size();};
			int getMoleculeCountOfType(MoleculeType *m);
			Molecule * getFirstMolecule() { return complexMembers.front(); };

			void mergeWithList(Complex * c);


			void updateComplexMembership(Molecule * m);


			void refactorToNewComplex(int new_ID_complex);

			void emptyComplexForever() {};

			static const int UNIFORM = 0;
			static const int FIXED_POINT = 1;
			static const int DIFFUSE_3D = 2;


			void printDegreeDistribution();
			void getDegreeDistribution(vector <int> &degreeDist);
			void printDetails();
			void printDetailsLong();

			//Diffusion functions
			double getDistance(Complex * c) {return 1000.0; };
			double getXpos() { return 0; };
			double getYpos() { return 0; };
			double getZpos() { return 0; };

			// get canonical label
			string getCanonicalLabel ( );
			// check if this is canonical
			bool isCanonical ( ) const { return is_canonical; };
			// unset canonical flag
			void unsetCanonical ( ) { is_canonical = false; };

			//This is public so that anybody can access the molecules quickly
			list <Molecule *> complexMembers;
			list <Molecule *>::iterator molIter;



		protected:
			// generate a canonical label using Nauty
			void   generateCanonicalLabel ( );

			System * system;
			int ID_complex;

			bool    is_canonical;
			string  canonical_label;

		private:

	};

	//!  Node object for creating labeled graphs from NFsim complexes.
	/*!
	    @author JustinHogg
	*/
	class Node
	{
	    public:
	        // constructor
	        Node ( Molecule * mol, int comp )
				: molecule(mol), component(comp), label( mol->getLabel(comp) ), index(-1) { };

	        // get methods
	        Molecule *  getMolecule ( ) const { return molecule; };
	        int         getComponent ( ) const { return component; };
	        string      getLabel ( ) const { return label; };
	        int         getIndex ( ) const { return index; };

	        // set methods
	        void  setIndex ( int val ) { index = val; };

	        // check if this is a molecule
	        bool  isMolecule ( ) const { return component == IS_MOLECULE; };

	        // compare methods
	        static bool less_by_label ( Node * node1, Node * node2 )
	        {   return node1->label < node2->label;   };

	        static bool less_by_index ( Node * node1, Node * node2 )
	        {   return node1->index < node2->index;   };

	        // definition
	        static const int IS_MOLECULE;

	    protected:

	        Molecule *  molecule;
	        int         component;
	        string      label;
	        int         index;
	};

	typedef  pair < Molecule *, int >  node_t;
	typedef  pair < node_t, Node * >   node_index_t;

}

#endif /*NFCORE_HH_*/
