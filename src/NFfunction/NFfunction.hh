#ifndef NFFUNCTION_HH_
#define NFFUNCTION_HH_


#include "muParser/muParser.h"
#include "../NFcore/NFcore.hh"


using namespace std;


namespace NFcore {

	class System;
	class ReactionClass;
	class Observable;
	class Complex; // this seems superfluous, but won't compile without it --Justin

	//! Parses mathmatical functions that can be easily used anywhere
	/*!
	    Built from the muParser freeware package, this factory creates Parser objects
	    that can be used to evaluate arbitrary functions.  The function to evaluate
	    is given as a string (which can contain the given functions or constants below,
	    case sensitive) along with pointers to particular varables that exist.  The Parser,
	    when initialized, translates the string into bytecode that can be evaluated fast.
	    At runtime, the bytecode that was created can be evaluated and every time it is
	    evaluated, it uses the current value of the given pointers to variables.  Thus,
	    NFsim can use functions of arbitrary complexity for defining reaction rates
	    and parameters.

	    <b>Supported functions:</b>
	    <table border=0>
	    	<tr><td>name&nbsp&nbsp</td><td>n_args&nbsp&nbsp</td><td>description</td></tr>
	        <tr><td>sin</td><td>1</td><td>sine function (arg is in radians)</td></tr>
    		<tr><td>cos</td><td>1</td><td>cosine function</td></tr>
			<tr><td>tan</td><td>1</td><td>tangens function</td></tr>
			<tr><td>asin</td><td>1</td><td>arcus sine function</td></tr>
			<tr><td>acos</td><td>1</td><td>arcus cosine function</td></tr>
			<tr><td>atan</td><td>1</td><td>arcus tangens function</td></tr>
			<tr><td>sinh</td><td>1</td><td>hyperbolic sine function</td></tr>
			<tr><td>cosh</td><td>1</td><td>hyperbolic cosine</td></tr>
			<tr><td>tanh</td><td>1</td><td>hyperbolic tangens function</td></tr>
			<tr><td>asinh</td><td>1</td><td>hyperbolic arcus sine function</td></tr>
			<tr><td>acosh</td><td>1</td><td>hyperbolic arcus tangens function</td></tr>
			<tr><td>atanh</td><td>1</td><td>hyperbolic arcur tangens function</td></tr>
			<tr><td>log2</td><td>1</td><td>logarithm to the base 2</td></tr>
			<tr><td>log10</td><td>1</td><td>logarithm to the base 10</td></tr>
			<tr><td>log</td><td>1</td><td>logarithm to the base 10</td></tr>
			<tr><td>ln</td><td>1</td><td>logarithm to base e (2.71828...)</td></tr>
			<tr><td>exp</td><td>1</td><td>e raised to the power of x</td></tr>
			<tr><td>sqrt</td><td>1</td><td>square root of a value</td></tr>
			<tr><td>sign</td><td>1</td><td>sign function -1 if x<0; 1 if x>0</td></tr>
			<tr><td>rint</td><td>1</td><td>round to nearest integer</td></tr>
			<tr><td>abs</td><td>1</td><td>absolute value</td></tr>
			<tr><td>if</td><td>3</td><td>if ... then ... else ...</td></tr>
			<tr><td>min</td><td>var.</td><td>min of all arguments</td></tr>
	 		<tr><td>max</td><td>var.</td><td>max of all arguments</td></tr>
			<tr><td>sum</td><td>var.</td><td>sum of all arguments</td></tr>
			<tr><td>avg</td><td>var.</td><td>mean value of all arguments</td></tr>
		</table>

	    <b>Supported constants</b> (to a precision of at least 10^-8):
	    <table border=0>
	    	<tr><td>name&nbsp&nbsp</td><td>value&nbsp&nbsp</td><td>description</td></tr>
	        <tr><td>_PI</td><td>3.141...</td><td>you all know pi, right?</td></tr>
	        <tr><td>_e</td><td>2.718...</td><td>the base of natural logs</td></tr>
	        <tr><td>_Na</td><td>6.022...e23</td><td>Avogadro's number</td></tr>
	    </table>

	    For muParser documentation, see:
	    http://muparser.sourceforge.net/

	*/
	class FuncFactory {
		public:

			/*!
			    Use this function to create a new function parser that can operate on the given
			    variables.  The vectors should contain (in the same order) the variable names
			    and pointers to the actual variable values.  Do not use local variables or variables
			    that will be forgotten before the Parser is destroyed.  To evaluate functions, you can
			    either just call the Eval() function of the Parser object (as in p->Eval()), or, call
			    the Eval function of FuncFactory to do error checking (as in FuncFactory::Eval(p)).
			*/
			static mu::Parser * create(string function, vector <string> & variableNames, vector <double *> & variablePtrs);


			/*!
				Creates a function without any variables and without any functions defined. You will have
				to do that yourself.  But, this does add in any predefined constants we want, and that is
				why this function exists.
			*/
			static mu::Parser * create();

			/*!
				Evaluates the given Parser object safely, meaning exceptions and errors are caught
				and program execution is terminated - if you have to do better than simple program
				termination, well than you can catch your own errors (see the muParser documenation).
			*/
			static double Eval(mu::Parser *p);

			/*!
				Runs a couple tests of the FuncFactory and the muParser functions to make sure everything
				is working as promised.  (To run this, use the arguements '-test mathFuncParser')
			*/
			static void test();

	};



	//! Defines functions to be used globally in a simulation.
	/*!
	    This small class is a small wrapper for the mu parser that allows the System
	    to keep track of functions that will be used throughout the course of the simulation.
	    These functions are only evaluated as needed to either output the value or
	    when recomputing the rate of some reaction.
	 */
	class GlobalFunction {

		public:
			/*!
				Creates a GlobalFunction with the given variables (which should be Observable objects
				that the System has) as well as a set of parameter values.  Note that creating a GlobalFunction
				does not initialize its parser.  The initialize, you have to call the prepareForSimulation() function
				which is currently handled by the System.
			*/
			GlobalFunction(string name,
					string funcExpression,
					vector <string> &varRefNames,
					vector <string> &varRefTypes,
					vector <string> &paramNames,
					System *s);

			/*!
				Deletes the GlobalFunction along with all of its associated variable and constant parameter information.
			*/
			~GlobalFunction();

			/*!
				This actually initializes the Function Parser so that it can be used.
			*/
			void prepareForSimulation(System *s);

			void updateParameters(System *s);



			/*!
				Simply gives the name of the function nicely (meaning something like func1()) for debugging / outputing.
			*/
			string getNiceName() const {return name+"()";};

			/*!
				Simply gives the name of the function only (without the open and close parentheses) .
			*/
			string getName() const {return name;};

			/*!
				For Debugging, prints out the details of the function including the defined variables and constant
				parameters as well as what the function currently evaluates to.
			*/
			void printDetails();
			void printDetails(System *s);


			/*!
				When a reaction uses this function, we have to keep track of it here so that
				we can notify it to update is propensity when this function changes...
			*/
			void attatchRxn(ReactionClass *r);


			int getNumOfVarRefs() const { return (int) n_varRefs; };
			string getVarRefName(int varRefIndex) const {
				return varRefNames[varRefIndex];
			}
			string getVarRefType(int varRefIndex) const {
				return varRefTypes[varRefIndex];
			}

			/*!
				This is the actual Parser object that keeps track of the function and has references to all of its
				arguments.  It is publicly visible, so be careful with it!  Use this variable to evaluate the function.
				See the mu parser documentation and the FuncFactory class for details on how to use this.
			*/
			mu::Parser *p;

		protected:

			string name;
			string funcExpression;

			unsigned int n_varRefs;
			string *varRefNames;
			string *varRefTypes;

			unsigned int n_params;
			string *paramNames;
	};



	class FunctionReference
	{

		public:
			FunctionReference(string name, string expression, string referencedFuncName) {
				this->name = name;
				this->expression=expression;
				this->referencedFuncName=referencedFuncName;
			};
			~FunctionReference() {};

			string name;
			string expression;
			string referencedFuncName;

	};




	class StateCounter {
		public:
			StateCounter(string name, MoleculeType *mt, string stateName);
			~StateCounter();

			void add(Molecule *m);
			void reset() { value=0; };
			int getValue() const { return (int)value; };


			string name;
			MoleculeType *mt;
			int stateIndex;
			double value;
	};









	class LocalFunction {

		public:

			LocalFunction(System *s,
								string name,
								string originalExpression,
								string parsedExpression,
								vector <string> &args,
								vector <string> &varRefNames,
								vector <string> &varObservableNames,
								vector <Observable *> & varObservables,
								vector <int> &varRefScope,
								vector <string> paramNames);
			~LocalFunction();


			string getName() const;
			string getNiceName() const;
			string getExpression() const;
			string getParsedExpression() const;

			// set/get whether this evaluates on complex complex
			bool getEvaluateComplexScope() const;
			void setEvaluateComplexScope( bool val );

			void printDetails(System *s);

			void prepareForSimulation(System *s);


			double getValue(Molecule *m, int scope);
			double evaluateOn(Molecule *m, int scope);
			// this version evaluates local fcn on a complex with SPECIES scope
			double evaluateOn(Complex *c);


			void addTypeIMoleculeDependency(MoleculeType *mt);
			void updateParameters(System *s);

			static const int SPECIES = 0;
			static const int MOLECULE = 1;

			mu::Parser *p;
		protected:

			// this variable defaults to false. If we detect that a rule requires
			// this function to evaluate on a species scope at any time, this is set
			// to true.  Otherwise, this allows us to just say zero if a typeII local
			// function molecule requests to be evaluated over the species scope, because,
			// by golly, that is never needed.  Is that clear?
			bool isEverEvaluatedOnSpeciesScope;

			string name;
			string nicename;
			string originalExpression;
			string parsedExpression;

			System * system;

			unsigned int n_args;
			string *argNames;


			unsigned int n_params;
			string *paramNames;

			unsigned int n_varRefs;
			string *varRefNames;
			string *varObservableNames;
			int *varRefScope;
			Observable **varLocalObservables;


			static list <Molecule *> molList;
			static list <Molecule *>::iterator molIter;

			//Here we store back pointers into both type I and type II molecules
			//Remember that type I molecules must store the value of this function
			//locally so that it can be used in DOR reactions.  Type II molecules
			//do not have the local value explicitly, but local functions should still
			//know 'of' them in case of future speedups that might use this information

			//@todo : change these to arrays from vectors!!!

			vector <MoleculeType *> typeI_mol;
			vector <int> typeI_localFunctionIndex;
			int n_typeIImolecules;
			MoleculeType ** typeII_mol;
			//vector <MoleculeType *> typeII_mol;
			vector <int> typeII_localFunctionIndex;


	};



	class CompositeFunction {
			public:
				CompositeFunction(System *s,
						string name,
						string expression,
						vector <string> &functions,
						vector <string> & argNames,
						vector <string> &paramNames);
				~CompositeFunction();


				string getName() const {return name;};

				void updateParameters(System *s);

				void finalizeInitialization(System *s);

				void prepareForSimulation(System *s);


				void setGlobalObservableDependency(ReactionClass *r, System *s);

				double evaluateOn(Molecule **molList, int *scope, int *curReactantCounts, int n_reactants);

				void printDetails(System *s);

				int getNumOfArgs() const;
				string getArgName(int aIndex) const;

				void addTypeIMoleculeDependency(MoleculeType *mt);


			protected:

				string name;
				string originalExpression;
				string parsedExpression;


				unsigned int n_allFuncs;
				string * allFuncNames;


				unsigned int n_args;
				string * argNames;

				unsigned int n_params;
				string *paramNames;

				//here are the referenced global functions
				int n_gfs;
				string * gfNames;
				GlobalFunction ** gfs;
				double * gfValues;

				//stores list of all local functions
				int n_lfs;
				string * lfNames;
				LocalFunction ** lfs;


				int n_reactantCounts;
				double * reactantCount;

				//stores list of all local functions and how they are referenced
				int n_refLfs;
				int *refLfInds;
				string *refLfRefNames;
				int *refLfScopes;
				double * refLfValues;

				mu::Parser *p;
		};





}














#endif /*NFFUNCTION_HH_*/
