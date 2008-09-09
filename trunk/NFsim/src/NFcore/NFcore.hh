//////////////////////////////////////////////////////////
// NFcore.hh
//
// Header file that contains the class declarations and definitions
// for the primary classes needed by the NFsim program.
//
//
// Michael Sneddon (michael.sneddon@yale.edu)
//
//////////////////////////////////////////////////////////
#ifndef NFCORE_HH_
#define NFCORE_HH_

#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <map>
#include <fstream>
#include <string>

#include "../NFutil/NFutil.hh"
#include "../NFoutput/NFoutput.hh"
#include "../NFreactions/NFreactions.hh"
#include "moleculeLists/moleculeList.hh"
#include "../NFfunction/NFfunction.hh"

#define DEBUG 0   			// Set to 1 to display all debug messages
#define BASIC_MESSAGE 1		// Set to 1 to display basic messages (eg runtime)

using namespace std;


//!  Contains the primary classes and functions of NFsim
/*!
  This namespace contains all the basic components needed for a typical simulation
  in NFsim.  The key classes are all defined in this header file for convenience, but
  the actual definitions are in the appropriately named cpp files.
  @author Michael Sneddon
*/
namespace NFcore 
{
	
	
	//Forward declarations to deal with cyclic dependencies 
//	class GroupOutputter;
	class MapGenerator;
	class MappingSet;
	class ReactantList;
	class TransformationSet;
	class MoleculeList;
	
	class GlobalFunction;
//	class TemplateMapping;
//	class Transformation;
	
	/*****************************************
	 * Class declarations
	 *    Here is the list of classes that are defined in this
	 *    header file and are used throughout the NFsim program
	 */
	
	class System;  /* the system that contains necessary information to set up, 
	                  run, and output results of the simulation */
	
	class MoleculeType;  /* indicates the "type" of molecule */
	class Molecule;  /* the actual class that represents instances of a molecule */
	class TemplateMolecule; /* used to determine if a molecule matches a particular pattern
	                           think of these as regular expressions for molecule comparison */
	
	class ReactionClass; /* defines a reaction class, (in other words, a rxn rule) */
	
	class Observable;  /* object that moniters counts of things we want to keep track of */
	
	class Complex;  /* collection of molecules that are bonded to each 
						other and whose membership dynamically can change*/					
//	class Group; /* collection of molecules that are not necessarily bonded together, but
//					can share certain properties */
//	class StateChangeListener;  /*used by molecules to let its group know its state changed */
	
	class ReactantList;
	
	


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
		public:
		
			static int NULL_EVENT_COUNTER;
			
			// Basic constructor and desconstructor methods
			System(string name);  /* creates a system with the given name that does not
				                         keep track dynamically of complex formation */
			System(string name, bool useComplex); /* creates a system with the given name
				                         that keeps track of complex formation if the useComplex
				                         parameter is set to true */
			~System(); /* destroys the system and cleans up all memory associated with it */
		
			// Basic functions to get the properties and objects of the system
			string getName() const { return name; };
			bool isUsingComplex() { return useComplex; };
			double getCurrentTime() const { return current_time; };
		
//			Group * getGroup(int ID_group) const { return allGroups.at(ID_group); };
//			int getGroupCount() const { return allGroups.size(); }
			int getObservableCount(int moleculeTypeIndex, int observableIndex) const;
			Observable * getObservableByName(string obsName);
			double getAverageGroupValue(string groupName, int valIndex);
		
			ReactionClass *getReaction(int rIndex) { return allReactions.at(rIndex); };
		
			MoleculeType * getMoleculeType(int mtIndex) { return allMoleculeTypes.at(mtIndex); };
			MoleculeType * getMoleculeTypeByName(string& name);
			int getNumOfMoleculeTypes() { return allMoleculeTypes.size(); };
			Molecule * getMoleculeByUid(int uid);
		    int getNumOfMolecules();
		
			//Functions used when setting up the system.  Note that most of these methods
			//are called automatically when you create an object (notable exception are
			//reactions), so use these methods with caution!  You probably don't need to
			//use them!
			int addMoleculeType(MoleculeType *moleculeType);
			void addReaction(ReactionClass *reaction);
			int createComplex(Molecule * m);
//			int addGroup(Group * g);
		
			void addGlobalFunction(GlobalFunction *gf);
			
			/* once all elements are added, we need to prepare and initialize for simulations */
			void prepareForSimulation();
		
		
			// A set of methods that allows us to update group properties for a subset of groups, 
			// or for all groups externally during the execution of a simulation */
			void updateGroupProperty(char * groupName, double *value, int n_values);
			void updateAllGroupProperty(double *value, int n_values);
		
		
			/* tell the system where to ouptut results*/
			void registerOutputFileLocation(string filename);
			void changeOutputFileLocation(string newFilename);
			
			void turnOnGlobalFuncOut() { this->outputGlobalFunctionValues=true; };
			void turnOffGlobalFuncOut() { this->outputGlobalFunctionValues=false; };
		
//			void addGroupOutputter(GroupOutputter * go);
		
			/* functions to output simulation properties to the registered file*/
			void outputAllObservableNames();
			void outputAllObservableCounts();
			void outputAllObservableCounts(double cSampleTime);
		
			/* functions that print out other information to the console */
			void printAllComplexes();
			void printAllReactions();
			void printAllGroups();
			void printIndexAndNames();
			void printAllMoleculeTypes();
			void printAllObservableCounts(double cSampleTime);
			void purgeAndPrintAvailableComplexList(); /* ONLY USE FOR DEBUG PURPOSES, AS THIS
															 DELETES ALL COMPLEX BOOKKEEPING */
			void outputComplexSizes(double cSampleTime);
			void outputMoleculeTypeCountPerComplex(MoleculeType *m);
			double outputMeanCount(MoleculeType *m);
			double calculateMeanCount(MoleculeType *m);
		
		
//			void outputGroupData() { outputGroupData(this->current_time); };
//			void outputGroupData(double cSampleTime);
		
		
		
		
		
		
			/* functions needed while running the simulation */
			Complex * getComplex(int ID_complex) const { return allComplexes.at(ID_complex); };
			Complex * getNextAvailableComplex();
			void notifyThatComplexIsAvailable(int ID_complex);
		
			/* run the simulation for a given length of time, output all results to the registered
			 * file.  The number of sample times is the number of times the function will output to
			 * a file divided equally throughout the elapsed time  */
			double sim(double time, long int sampleTimes);
		
			/* run the simulation up until the stopping time (but not exceding the stopping time. This
			 * will not output anything to file (so must be done manually) and returns the current time
			 * of the simulation, which will always be less than the stopping time */
			double stepTo(double stoppingTime);
			
			void singleStep();
		
			/* runs the simulation without output for the given time.  After, the clock is reset to
			 * zero so that additional simulations can take place */
			void equilibriate(double duration);
			void equilibriate(double duration, int statusReports);
			
			
			//For moleculeTypes to do a quick lookup and see where a particular reaction
			//is mapped to.  This is an optimization....
			void registerRxnIndex(int rxnId, int rxnPos, int rxnIndex) {rxnIndexMap[rxnId][rxnPos]=rxnIndex;};
			int getRxnIndex(int rxnId, int rxnPos) const { return rxnIndexMap[rxnId][rxnPos]; };
		
		protected:
		
			// The invariant system properties, created when the system is created
			string name;         /* arbitrary name of the system  */
			bool useComplex;     /* parameter that knows if we should be dynamically tracking complexes */
		    
		
			// The container objects that maintain the system configuration
			vector <MoleculeType *> allMoleculeTypes;  /* container of all MoleculeTypes in the simulation */
			vector <ReactionClass *> allReactions;    /* container of all Reactions in the simulation */
//			vector <Group *> allGroups;               /* container of all groups in the simulation */
			vector <Complex * > allComplexes;         /* container of all complexes in the simulation */
			queue <int> nextAvailableComplex;         /* queue tells us which complexes can be used next */
		
			
			vector <GlobalFunction *> globalFunctions; /* container of all functions available to the system */
		
			// Properties of the system that update in time
			double a_tot;        /* the sum of all a's (propensities) of all reactions */
			double current_time; /* keep track of the simulation time */
			ReactionClass * nextReaction;  /* keep track of the next reaction to fire */
		
		
			// functions needed only by the system while running a simulation
			double get_A_tot() const { return a_tot; };
			double recompute_A_tot();
			double getNextRxn();
		
		
			// Neccessary variables and methods for outputting
			ofstream outputFileStream; /* the stream to a file to write out the results */
			void outputGroupDataHeader();
//			GroupOutputter * go;
		    bool outputGlobalFunctionValues;
			
			//For moleculeTypes to do a quick lookup and see where a particular reaction
			//is mapped to.  This is an optimization....
			int **rxnIndexMap;
			
		
		private:
		
			//Iterators that allow fast traversal of the object containers
			vector<MoleculeType *>::iterator molTypeIter;  /* to iterate over allMoleculeType */
			vector <ReactionClass *>::iterator rxnIter;    /* to iterate over allReactions */
			vector <Complex *>::iterator complexIter;      /* to iterate over allComplexes */
//			vector <Group *>::iterator groupIter;          /* to iterate over allGroups */
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
					vector <string> &defaultCompState,
					vector < vector<string> > &possibleCompStates,
					System *system);
			
			~MoleculeType();
		
			string getName() const { return name; };
			System * getSystem() const { return system; };
			
			int getNumOfComponents() const { return numOfComponents; };
			string getComponentName(int cIndex) const { return compName[cIndex]; };
			void getPossibleComponentStates(int cIndex, list <string> &nameList);
			int getDefaultComponentState(int cIndex) const { return defaultCompState[cIndex]; };
			
			int getCompIndexFromName(string cName) const;
			string getComponentStateName(int cIndex, int cValue);
			
			
			
			
			void addEquivalentComponents(vector <vector <string> > &identicalComponents);
			bool isEquivalentComponent(string cName) const;
			bool isEquivalentComponent(int cIndex) const;
			void getEquivalencyClass(int *&components, int &n_components, string cName) const;
			
			
			
			
			
			
			/* get functions */
			
			int getNumOfObservables() const { return observables.size(); };
			string getObservableAlias(int obsIndex) const;
			Observable * getObservable(int obsIndex) const { return observables.at(obsIndex); };
			unsigned long int getObservableCount(int obsIndex) const;
		
			int getReactionCount() const { return reactions.size(); };
			int getTypeID() const { return type_id; };
		
			Molecule * getMolecule(int ID_molecule) const;
			int getMoleculeCount() const;
		
			int getRxnIndex(ReactionClass * rxn, int rxnPosition);
		
			
			//Functions for adding and removing molecules from the system during a simulation
			//Molecule *createMolecule();
			//void prepMolecule();
			//void deleteMolecule(Molecule *m);
			//void deleteJustOneMolecule(Molecule *m);
			
		
			/* set functions */
			Molecule *genDefaultMolecule();
			
			void addMoleculeToRunningSystem(Molecule *&mol);
			void removeMoleculeFromRunningSystem(Molecule *&m);
			
			void removeFromRxns(Molecule * m);
			//int addMolecule(Molecule *m);
			
			
			
			
			
			
			void addReactionClass(ReactionClass * r, int rPosition);
			void addObservable(Observable * o) { observables.push_back(o); }; //could add check here to make sure observable is of this type
			int createComplex(Molecule *m) { return system->createComplex(m); };
			void addTemplateMolecule(TemplateMolecule *t);
		
			/* handle DOR reactions */
			void addDORrxnClass(ReactionClass * r, int rPosition);
			int getDORrxnCount() const { return indexOfDORrxns.size(); };
			ReactionClass * getDORrxn(int index) const { return reactions.at(indexOfDORrxns.at(index)); };
			int getDORreactantPosition(int index) const { return reactionPositions.at(indexOfDORrxns.at(index)); };
			int getDORreactantIndex(int index) const { return indexOfDORrxns.at(index); };
		
			/* updates a molecules membership (assumes molecule is of type this) */
			void updateRxnMembership(Molecule * m);
		
			/* functions which handles observables */
			void removeFromObservables(Molecule * m);
			void addToObservables(Molecule * m);
			void outputObservableNames(ofstream &fout);
			void outputObservableCounts(ofstream &fout);
			void printObservableNames();
			void printObservableCounts();
		
			
			void printDetails() const;
		
			/* auto populate with default molecules */
			void populateWithDefaultMolecules(int moleculeCount);
		
			/* this method assumes all molecules in the simulation
			 * have been defined, and all reaction classes and observables
			 * have been added.  Then, this function will add those
			 * molecules to rxn lists and instantiate the counts of the observables. 
			 * In general, you do not need to worry about this function because
			 * it automatically gets called by the System when you prepare the System*/
			void prepareForSimulation();
		
		protected:
		
			System *system;
								
						string name;
						int type_id;
								
									
						//keeps track
						int numOfComponents;
						string *compName;
						vector < vector < string > > possibleCompStates;
						int *defaultCompState;
						
						////////////////////////////////////////////////////////////////////////////////
						int n_eqComp;
						int * eqCompSizes;
						string **eqCompName;
						int **eqCompIndex;
		
			MoleculeList * mList;
			//vector <Molecule *> mInstances;  /* List of all molecules that exist */
			vector <ReactionClass *> reactions; /* List of reactions that this type can be involved with */
			vector <int> reactionPositions;   /* the position in the reaction for this type of molecule */
		
			vector <int> indexOfDORrxns;
		
			vector <Observable *> observables;  /* list of things to keep track of */
		
			vector <TemplateMolecule *> allTemplates; /* keep track of all templates that exist of this type
															so that they are easy to delete from memory at the end */
			
			
			
			
			
			
			
			
			//Maintains the set of equivalent sites used only when creating reaction
			//rules
//			unsigned int n_eqSites;
//			unsigned int * eqSiteSizes;
//			string **eqSiteName;
//			int **eqSiteIndex;
//			
//			int *isSiteEq;
//			
//			
//			unsigned int n_eqStates;
//			unsigned int * eqStateSizes;
//			string **eqStateName;
//			int **eqStateIndex;
//					
//			int *isStateEq;
		
		private:
			vector<Molecule *>::iterator molIter;  /* to iterate over mInstances */
			vector<Observable *>::iterator obsIter; /* to iterate over observables */
			vector <ReactionClass *>::iterator rxnIter; /* to iterate over reactions */
	};
	
	
	
	//!  Each molecule in the system is represented by an instance of this.
	/*!
      The base unit of the NFsim program, this class maintains instances of
	  individual objects that exist in the simulation.  All molecules that exist
	  must be generated before the simulation starts.  Molecules are then able
	  to interact with thier MoleculeType so that they also 'know' about the
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
			int getMoleculeID() const { return ID_number; };
			int getMolListId() const { return listId; };
			string getMoleculeTypeName() const { return parentMoleculeType->getName(); };
			MoleculeType * getMoleculeType() const { return parentMoleculeType; };
			int getUniqueID() const { return ID_unique; };
		
			int getComplexID() const { return ID_complex; };
			Complex * getComplex() const { return parentMoleculeType->getSystem()->getComplex(ID_complex); };
			int getDegree();
		
			
			///////////////////////////////////////////////////////////////////////
			int getComponentState(int cIndex) const { return component[cIndex]; };
			int getComponentIndexOfBond(int cIndex) const { return indexOfBond[cIndex]; };
			void setComponentState(int cIndex, int newValue);
			void setComponentState(string cName, int newValue);
			
			
			
			////////////////////////////////////////////////////////////////////
			
			
			//int getState(char * stateName) const { return 0; }
			//int getState(int stateIndex) const { return states[stateIndex]; };
		
			/* accessor functions for checking binding sites */
			bool isBindingSiteOpen(int bIndex) const;
			bool isBindingSiteBonded(int bIndex) const;
			Molecule * getBondedMolecule(int bSiteIndex) const;
			//int getBsiteIndexOfBond(int bSiteIndex) const { return bSiteIndexOfBond[bSiteIndex]; };
		
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
		//	void clear();
		//	void setHasVisited(int bSiteIndex);
			void traverseBondedNeighborhood(list <Molecule *> &members, int traversalLimit);
			static void breadthFirstSearch(list <Molecule *> &members, Molecule *m, int depth);
		
			/* when we are ready to begin simulations, moleculeType calls this function
			 * so that this molecule can add itself to all the necessary lists */
			void prepareForSimulation();
		
			/* function that tells this molecule that it changed states or bonds
			 * and it should update its reaction membership */
			void updateRxnMembership();
			void removeFromObservables();
			void addToObservables();
			void updateDORs();
		
			double getDORvalueFromGroup(string groupName, int valueIndex);
		
		
		
			/* group functions */
//			void addListener(StateChangeListener * l) { listeners.push_back(l); };
		
			void notifyGroupsThatRateMayChange();
		
		
			/* print functions for debugging */
			void printDetails() const;
			static void printMoleculeList(list <Molecule *> &members);
		
			static int getUniqueIdCount() { return uniqueIdCount; };
			static const int NOT_IN_RXN = -1;
		
			
			/* used for traversing a molecule complex */
			//bool * hasVisitedBond;
			bool hasVisitedMolecule;
			
			
			static const int NOSTATE = -1;
			static const int NOBOND = 0;
			static const int NOINDEX = -1;
			
		protected:
		

			bool isPrepared;
			
			
			/* Set of IDs which identifies uniquely this molecule */
			int ID_number;
			int ID_complex;
			int ID_type;
			int ID_unique;
			int listId;
		
//			list <StateChangeListener *> listeners;
		
			static int uniqueIdCount;
		
			/* The type of this molecule */
			MoleculeType *parentMoleculeType;
			bool useComplex;
		
		
			/* store the states and bonds in arrays */
			
			
			///////////////////////////////////////////////////////////////////
			/* list of components */
			int *component;
			int numOfComponents;
			Molecule **bond;
			int *indexOfBond; /* gives the index of the component that is bonded to this molecule */
			
			
			
			//int *states;
			//Molecule ** bonds;
			//int * bSiteIndexOfBond;  // index of this molecule in bonded 
		
			
			/* used for traversing a molecule complex */
			int * rxnListMappingId;
			int nReactions;
			
			
		private:
			
			static queue <Molecule *> q;
			static queue <int> d;
			static list <Molecule *>::iterator molIter;
			
//			list <StateChangeListener *>::iterator listenerIter; /* to iterate over groups this molecule is in */
	};
	
	
	

	
	
	
	
	
	
	//!  Abstract Base Class that defines the interface for all reaction rules.
	/*!
	    @author Michael Sneddon
	*/
	class ReactionClass
	{
		public:
			static const int NO_LIMIT = -3;
			
			static const int BASIC_RXN = 0;
			static const int DOR_RXN = 1;
			static const int OBS_DEPENDENT_RXN = 2;
			
			

			ReactionClass(string name, double rate, TransformationSet *transformationSet);
			virtual ~ReactionClass();
			
			int getNumOfReactants() const { return n_reactants; };
			
			string getName() const { return name; };
			double getBaseRate() const { return baseRate; };
			int getRxnType() const { return reactionType; };

			void setBaseRate(double newBaseRate) { this->baseRate=newBaseRate; update_a(); };
			void setTraversalLimit(int limit) { this->traversalLimit = limit; };
			
			double get_a() const { return a; };
			void printDetails() const;
			void fire(double random_A_number);
			
			
			
			//The main virtual functions that must be implemented in all implementing classes
			virtual void init() = 0; //called when the reaction is added to the system
			virtual void prepareForSimulation() = 0; //called once everything is added to the system
			virtual bool tryToAdd(Molecule *m, unsigned int reactantPos) = 0;
			virtual void remove(Molecule *m, unsigned int reactantPos) = 0;
			
			virtual double update_a() = 0;
			
			
			virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) = 0;
			
			
			
			virtual unsigned int getReactantCount(unsigned int reactantIndex) const = 0;
			virtual void printFullDetails() const = 0;
			
			
			void setRxnId(int rxnId) { this->rxnId = rxnId; };
			int getRxnId() const { return rxnId; };
			
		protected:
			virtual void pickMappingSets(double randNumber) const=0;
			
			int rxnId;
			
			string name;
			int reactionType;
			unsigned int n_reactants;
			
			double baseRate;
			double a;
			unsigned int fireCounter;
			
			unsigned int traversalLimit;
			
			TemplateMolecule **reactantTemplates;
			TransformationSet * transformationSet;
			MappingSet **mappingSet;
			
			
			
			list <Molecule *> products;
			list <Molecule *>::iterator molIter;
	};
	
	
	
	
	
	
	
	
	//!  Used for matching Molecule objects to a given pattern
	/*!
        @author Michael Sneddon
	 */
	class TemplateMolecule {
	public:
		TemplateMolecule(MoleculeType * parentMoleculeType);
		~TemplateMolecule();
	
		/* get functions */
		MoleculeType * getMoleculeType() const { return parentMoleculeType; };
		string getMoleculeTypeName() const { return parentMoleculeType->getName(); };
	
		unsigned int getNumCompareStates() const { return stateIndex.size(); };
		unsigned int getStateIndex(int state) const { return stateIndex.at(state); };
		unsigned int getStateValue(int state) const { return stateValue.at(state); };
	
		unsigned int getNumBindingSites() const { return bonds.size(); };
		bool isBindingSiteOpen(int bIndex) const;// { return bonds.at(bIndex)->isOpen(); };
		bool isBindingSiteBonded(int bIndex) const; //{ return bonds.at(bIndex)->isBonded(); };
		TemplateMolecule * getBondedTemplateMolecule(int bIndex) const;
		unsigned int getTemplateBsiteIndexFromMoleculeBsiteIndex(int molBsiteIndex);
	
		int getBindingSiteIndex(int bIndex) const;
		int getBindingSiteIndexOfBondedTemplate(int bIndex) const ;
		
		
		/* set functions */
		void setHasVisited(int bSite);
		unsigned int addEmptyBindingSite(int bSiteIndex);
		unsigned int addEmptyBindingSite(const char * bSiteName);
		void addOccupiedBindingSite(const char * bSiteName);
		static void bind(TemplateMolecule *t1, int bSiteIndex1, TemplateMolecule *t2, int bSiteIndex2);
		static void bind(TemplateMolecule *t1, const char * bSiteName1, TemplateMolecule *t2, const char * bSiteName2);
	
		void addStateValue(int stateIndex, int stateValue);
		void addStateValue(string stateName, int stateValue);
		void addNotStateValue(char * stateName, int notStateValue);
		void clear() { this->matchMolecule = 0; for(unsigned int i=0; i<hasVisitedBond.size(); i++) hasVisitedBond.at(i) = false; hasVisited=false; };
	
	
		/* the primary function and purpose of a template molecule 
			   is to compare itself to an instance of a molecule */
		bool compare(Molecule * m);
		static bool compareBreadthFirst(TemplateMolecule *tm, Molecule *m);
	
		/*
		 * Used to check if a particular state value matches - used when parsing an xml file and
		 * generating reaction transformations.
		 */
		bool isStateValue(const char * stateName, int stateValue);
		bool isBonded(const char * bSiteName);
	
	
	//	void addTemplateMapping(TemplateMapping *tm);
	//	bool compare(Molecule * m, MappingSet *mappingSet);
	
	
	
		void printDetails() const;
	
	
		
		
		///////////////////////////////////////////////////////////////////
		void addMapGenerator(MapGenerator *mg);
		bool compare(Molecule *m, MappingSet *ms);
		bool contains(TemplateMolecule *tempMol);
	
	
		Molecule * matchMolecule;
		vector <bool> hasVisitedBond;
				
				
		/////////////////////////////////////////////////////////
		bool hasVisited;
		
	protected:
	
	
	
		MoleculeType *parentMoleculeType;  // ptr to indicate type of molecule
		vector <int> stateIndex; // saves index into state array of MoleculeType
		vector <int> stateValue; // saves the value we need to have in the state array
	
		vector <int> notStateIndex;
		vector <int> notStateValue;
	
		vector <int> bSiteIndex; // saves index into bSite array in MoleculeType
		vector <TemplateMolecule *> bonds; // tells us a binding site
		vector <int> bSiteIndexOfBond;
	
		vector <int> sitesThatMustBeOccupied; // 
	

		
		vector <MapGenerator *> mapGenerators;
		vector <MapGenerator *>::iterator mgIter;
		vector <int>::iterator intVecIter;
		
		
		static queue <TemplateMolecule *> tmq;
		static queue <Molecule *> mq;
		static list <TemplateMolecule *> tml;
		static queue <int> d;
		static list <TemplateMolecule *>::iterator tmIter;
		
	};
	
	
	

	//!  Tracks the counts of predefined observables in the simulation.
	/*!
        @author Michael Sneddon
	 */
	class Observable
	{
		public:
			Observable(string aliasName, TemplateMolecule * templateMolecule);
			~Observable();
		
			/* methods used to keep the observable count up to date */
			bool isObservable(Molecule * m) const;
			void add() { count++;};
			void subtract() { if(count==0){ cout<<"Error in observable count!!"<<endl; exit(1); } count--;};
		
			/* methods used to get observable information */
			unsigned long int getCount() const {return (unsigned long int) count;};
			string getAliasName() const { return aliasName; };
			
			void addReferenceToMyself(mu::Parser * p);
		
		protected:
			string aliasName;   /* The name that will be output for this observable */
			TemplateMolecule * templateMolecule; /* The template molecule which represents what we want to observe */
			double count; /* the number of molecules that match this observable, its a double so that functions can use it (as a reference) */
		
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
		
			int getComplexID() const { return ID_complex; };
			int getComplexSize() const {return complexMembers.size();};
			int getMoleculeCountOfType(MoleculeType *m);
		
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
		
		protected:
			list <Molecule *> complexMembers;
			System * system;
			int ID_complex;
		
			int diffusionType;
			double x_pos;
			double y_pos;
			double z_pos;
		
		private:
			list <Molecule *>::iterator molIter;
	
	};
	
	
	
	//!  A defined set of Molecule objects that have some joint property
	/*!
        @author Michael Sneddon
	 */
//	class Group
//	{
//		public:
//			Group(string groupName, System * s, int stateIndex);
//			virtual ~Group();
//		
//			string getName() const { return groupName; };
//			double getValue(unsigned int valIndex);
//			int getNumberInGroup() const { return groupMembers.size(); };
//			Molecule * getMolecule(int mIndex) { return groupMembers.at(mIndex); };
//		
//			virtual void addToGroup(Molecule * m);
//			virtual void addToGroupWithoutListener(Molecule * m);
//		
//			virtual void notify(Molecule *changedMolecule, int oldStateValue, int newStateValue);
//		
//			/* an function that lets classes that inherit from group update
//			 * some of its values that are not dependent on state, but rather
//			 * depend on some outside value at arbitrary points in the simulation */
//			virtual void updateGroupProperty(double * values, int n_values);
//		
//			virtual void updateReactionRates();
//		
//			virtual void printDetails();
//		
//		protected:
//			vector <Molecule *> groupMembers;
//			System * system;
//		
//			string groupName;
//			int Group_ID;
//		
//			vector <char *> valueNames;
//			vector <double> value;  
//			/* keeps track of the sum of the state values - 
//				                  override notify to change what this does */
//			int stateIndex;
//		
//			bool areReactionsUpToDate;
//		
//			vector <Molecule *>::iterator molIter;	
//	};
	
	
	
	//!  Listens to state changes in a Molecule and reports back to a Group
	/*!
        @author Michael Sneddon
	 */
//	class StateChangeListener
//	{
//		public:
//			StateChangeListener(Molecule * m, Group * g, int stateIndex);
//			~StateChangeListener();
//		
//			/* called by molecule to alert that state has changed */
//			void notify(Molecule *changedMolecule, int stateIndexOfChange);
//		
//			/* called by molecule when the molecule is ready to have
//			 * all the group's reactions updated */
//			void updateGroupReactions();
//		
//			/* simple functions for molecule to talk to this group */
//			string getGroupName() const { return group->getName(); };
//			double getValue(int valueIndex) const { return group->getValue(valueIndex); };
//		
//		
//		protected:
//			Molecule * molecule;
//			Group * group;
//			int oldStateValue;
//			int newStateValue;
//		
//			int stateIndex;
//	
//	};
	
	
}

#endif /*NFCORE_HH_*/
