//////////////////////////////////////////////////////////
// NFutil.h
//
// Header file that contains the class declarations and definitions
// for any math utility functions that are needed for running the
// NFsim program.  Currently, this just includes code to generate
// random numbers.
//
// Michael Sneddon (michael.sneddon@yale.edu)
//
//////////////////////////////////////////////////////////
#ifndef NFUTIL_H_
#define NFUTIL_H_


namespace NFutil {

	/* Seed the random number generator with a positive 32bit integer
	 * If you don't call this function, the current time will be used
	 * as a seed so each run will be different */
	void SEED_RANDOM( unsigned long  seed );
	
	/* This function returns a random double on the half open interval
	 * (0,max].  This is important for selecting the next reaction
	 * to fire because we want the value to always be greater than
	 * zero, but the value can equal the max value. */
	double RANDOM( double max );
	
	/* This function returns a random double on the closed interval
	 * (0,1).  This is important for choosing the next dt in the simulation
	 * because we will take the log of this number and we don't want
	 * dt=0 or dt=infinity. */
	double RANDOM_CLOSED();
	
	/* This is a multipurpose function that selects a random number
	 * on the half open interval [min, max).  This allows you to quickly
	 * select a random element from a list of known size.  The function
	 * call would look something like: RANDOM_INT(0,sizeOfList).  This
	 * is used to select the next molecule to be included in a reaction. */
	int RANDOM_INT(unsigned long min, unsigned long max);
	
	
	/* Returns a normally distributed random number. 
	 * Based on the java implementation of Random.nextGaussian, this method uses the
	 * polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia, as described by 
	 * Donald E. Knuth in The Art of Computer Programming, Volume 2: Seminumerical 
	 * Algorithms, section 3.4.1, subsection C, algorithm P.
	 */
	double RANDOM_GAUSSIAN();
}
#endif /*NFUTIL_H_*/
