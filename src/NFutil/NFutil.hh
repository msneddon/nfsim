//////////////////////////////////////////////////////////
// NFutil.h
//
// Header file that contains the class declarations and definitions
// for any math utility functions that are needed for running the
// NFsim program.
//
// Michael Sneddon (michael.sneddon@yale.edu)
//
//////////////////////////////////////////////////////////
#ifndef NFUTIL_H_
#define NFUTIL_H_

// JUSTIN -- added stdlib.h to resolve problem finding "exit"
#include <stdlib.h>
#include <string>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>



using namespace std;

//!  General utility functions for NFsim, including a random number generator and a function parser.
/*!
    @author Michael Sneddon
 */
namespace NFutil {


	//!  Seeds the random number generator used in all simulations
	/*!
	   Seed the random number generator with a positive 32bit integer
	   If you don't call this function, the current time will be used
	   as a seed so each run will be different.
    	@author Michael Sneddon
	 */
	void SEED_RANDOM( unsigned long  seed );


	//!  Uniform random number on the interval (0,max]
	/*!
		This function returns a random double on the half open interval
		(0,max].  This is important for selecting the next reaction
		to fire because we want the value to always be greater than
		zero, but the value can equal the max value.
	    @author Michael Sneddon
	*/
	double RANDOM( double max );


	//!  Uniform random number on the interval (0,1)
	/*!
		This function returns a random double on the open interval
		(0,1).  This is important for choosing the next dt in the simulation
		because we will take the log of this number and we don't want
		dt=0 or dt=infinity.
		    @author Michael Sneddon
	*/
	double RANDOM_OPEN();

	//!  Uniform random number on the closed interval [0,1]
	/*!
		This function returns a random double on the closed interval [0,1].
		    @author Michael Sneddon
	*/
	double RANDOM_CLOSED();


	//!  Uniform random number on the interval [min,max]
	/*!
		This is a multipurpose function that selects a random number
		on the half open interval [min, max).  This allows you to quickly
		select a random element from a list of known size.  The function
		call would look something like: RANDOM_INT(0,sizeOfList).  This
		is used to select the next molecule to be included in a reaction.
		    @author Michael Sneddon
	*/
	int RANDOM_INT(unsigned long min, unsigned long max);

	//!  Normally distributed random number with mean 0 and variance 1.
	/*!
		Returns a normally distributed random number.
		Based on the java implementation of Random.nextGaussian, this method uses the
		polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia, as described by
		Donald E. Knuth in The Art of Computer Programming, Volume 2: Seminumerical
		Algorithms, section 3.4.1, subsection C, algorithm P.
			    @author Michael Sneddon
	*/
	double RANDOM_GAUSSIAN();



	//!  Parses and converts std::string objects to double values.
	/*!
		Just uses the standard library operations for converting strings to
		double values.  Throws a run time exception if the parse failed.
			@author Michael Sneddon
	*/
	double convertToDouble(const std::string& s);

	//!  Parses and converts std::string objects to int values.
	/*!
		Just uses the standard library operations for converting strings to
		integers.  Throws a run time exception if the parse failed.
			@author Michael Sneddon
	*/
	int convertToInt(const std::string& s);


	//!  Converts double values to their string representations
	/*!
			@author Michael Sneddon
	*/
	string toString(double x);

	//!  Converts integer values to their string representations
	/*!
			@author Michael Sneddon
	*/
	string toString(int x);


	//! Removes leading and trailing whitespace (including tabs) from the string
	/*!
	 */
	void trim(string& str);





	const double PI = 3.14159265358979323846;
	const double NA = 6.02214179e23;

}




#endif /*NFUTIL_H_*/
