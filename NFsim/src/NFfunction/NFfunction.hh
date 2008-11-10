#ifndef NFFUNCTION_HH_
#define NFFUNCTION_HH_


#include "muParser/muParser.h"
#include "../NFcore/NFcore.hh"


using namespace std;


namespace NFcore {

	class System;
	class ReactionClass;
	class Observable;

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

	//double FuncFactory::PI = 3.14159265358979323846;
	//double FuncFactory::NA ;
	//double FuncFactory::E ;








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
					string funcString,
					vector <string> &argNames,
					vector <string> &argTypes,
					vector <string> &paramConstNames,
					vector <double> &paramConstValues);

			/*!
				Deletes the GlobalFunction along with all of its associated variable and constant parameter information.
			*/
			~GlobalFunction();

			/*!
				This actually initializes the Function Parser so that it can be used.
			*/
			void prepareForSimulation(System *s);

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


			/*!
				When a reaction uses this function, we have to keep track of it here so that
				we can notify it to update is propensity when this function changes...
			*/
			void attatchRxn(ReactionClass *r);


			int getNumberOfArgs() const { return (int) n_args; };
			string getArgName(int argIndex) const {
				if((unsigned)argIndex<n_args && argIndex>=0) return argNames[argIndex];
				cerr<<"invalid argIndex given in GlobalFunction."<<endl; exit(1); };
			string getArgType(int argIndex) const {
				if((unsigned)argIndex<n_args && argIndex>=0) return argTypes[argIndex];
				cerr<<"invalid argIndex given in GlobalFunction."<<endl; exit(1); };

			/*!
				This is the actual Parser object that keeps track of the function and has references to all of its
				arguments.  It is publicly visible, so be careful with it!  Use this variable to evaluate the function.
				See the mu parser documentation and the FuncFactory class for details on how to use this.
			*/
			mu::Parser *p;

		protected:

			string name;
			string funcString;

			unsigned int n_args;
			string *argNames;
			string *argTypes;

			unsigned int n_paramConst;
			string *paramNames;
			double *paramValues;
	};



	class StateCounter {
		public:
			StateCounter(string name, MoleculeType *mt, string stateName);
			~StateCounter();

			void add(Molecule *m);
			void reset() { value=0; };
			int getValue() const { return value; };


			string name;
			MoleculeType *mt;
			int stateIndex;
			double value;
	};


	class LocalFunction {

		public:
			LocalFunction(string name,
					string funcString,
					vector <Observable *> &observables,
					vector <StateCounter *> &stateCounters,
					vector <string> &paramConstNames,
					vector <double> &paramConstValues);
			~LocalFunction();

		//	void prepareForSimulation(System *s);

			string getNiceName() const {return name+"()";};

			string getName() const {return name;};

			void printDetails();

		//	void attatchRxn(ReactionClass *r);

			double evaluateOn(Molecule *m);


//			int getNumberOfArgs() const { return (int) n_args; };
//			string getArgName(int argIndex) const {
//				if(argIndex<n_args && argIndex>=0) return argNames[argIndex];
//				cerr<<"invalid argIndex given in GlobalFunction."<<endl; exit(1); };
//			string getArgType(int argIndex) const {
//				if(argIndex<n_args && argIndex>=0) return argTypes[argIndex];
//				cerr<<"invalid argIndex given in GlobalFunction."<<endl; exit(1); };

			mu::Parser *p;
		protected:
			//List of observables that this local function depends on
			Observable ** obs;
			unsigned int n_obs;
			int * obsVal;

			StateCounter ** sc;
			unsigned int n_sc;

			string name;
			string funcString;

			unsigned int n_paramConst;
			string *paramNames;
			double *paramValues;


	};




}














#endif /*NFFUNCTION_HH_*/
