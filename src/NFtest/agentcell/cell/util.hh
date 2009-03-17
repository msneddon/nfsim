#ifndef CHEMOTAXIS_UTIL_HH_
#define CHEMOTAXIS_UTIL_HH_

#include <math.h>
#include "../../../NFutil/NFutil.hh"




using namespace NFutil;


/** Contains some classes / functions that are useful for chemotaxis models.
 * 
 *  
 */
namespace ChemotaxisUtil {
	

/** Generates random numbers sampled from a gamma distribution.
 * 
 *  A nice implementation of an algorithm that samples from a gamma
 *  distributionwith code  modified from the book :
 * 
 *  Press, Teukolsky, Vetterling, and Flannery.  "Numerical Recipes, The art
 *  of Scientific Computing", Third Edition. 2007 Cambridge University Press.
 * 
 *  It additionally has a parameter that allows the entire distribution
 *  to be shifted by a given amount and a generate function that resamples
 *  until the chosen number is within a particular range.  In the chemotaxis
 *  model, this code is used primarily to choose a new random angle from the
 *  distribution of changes in direction (fig 3 from Berg & Brown, Chemotaxis
 *  in Escherichia coli analysed by Three-dimensional Tracking, Nature,
 *  volume 239, (1972) )
 */
class GammaSampler {
	
	public:
		/**
		 * Constructor to generate a basic gamma distribution.
		 * @param alpha the alpha parameter (note, alpha=k=shape).
		 * @param beta the beta parameter (note, beta=1/theta=1/scale).
		 */
		GammaSampler(double alpha, double beta);
		
		/**
		 * Constructor to generate a gamma distribution that is shifted by the given amount.
		 * @param alpha the alpha parameter (note, alpha=k=shape).
		 * @param beta the beta parameter (note, beta=1/theta=1/scale).
		 * @param meanShift shifts the entire distribution by this amount.
		 */
		GammaSampler(double alpha, double beta, double meanShift);
		
		
		/**
		 * Generates a new random number (using the NFutil Mersenne Twister random 
		 * number generator) sampled from the initialized gamma distribution.
		 */
		double gen();
		
		/**
		 * Generates a new random number (using the NFutil Mersenne Twister random 
		 * number generator) sampled from the initialized gamma distribution that
		 * lies within (inclusive) of the range given.  Carefull here!  Your range
		 * should be such that it doesn't take forever to find a number in the range!
		 * @param minAcceptableValue the minimum random value this function will return.
		 * @param maxAcceptableValue the maximum random value this function will return.
		 */
		double gen(double minAcceptableValue, double maxAcceptableValue);
	
	protected:
		double alpha;
		double oalpha;
		double beta;
		double a1;
		double a2;
		double meanShift;
	
	private:
		void init(double alpha, double beta, double meanShift);
};


// 2) random rotation matrix generator


void genUniformRandRotation3d(double rotMatrix[3][3]);
void genRot3dAboutAxis(double rotMatrix[3][3], double axis[3], double angleInRad);
void genRotFromAngles(double rotMatrix[3][3], double angleX, double angleY, double angleZ);

void applyRotation(double rotMatrix[3][3], double vec[3]);
void applyDisplacement(double pos[3], double dir[3], double distance);


}

#endif /*CHEMOTAXIS_UTIL_HH_*/
