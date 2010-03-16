#include "util.hh"
#include <iostream>

using namespace ChemotaxisUtil;
using namespace std;



/*
 * Class that generates random numbers drawn from a gamma distribution
 * with alpha and beta parameters.  The mean of such a distribution is alpha / beta
 * and the variance is alpha / beta^2.  Often, parameters for the gamma distribution
 * are given as the shape (k) and the scale (theta) in which case alpha = k and
 * beta = 1/theta.  
 * 
 */
GammaSampler::GammaSampler(double alpha, double beta)
{
	init(alpha, beta, 0.0);
}

GammaSampler::GammaSampler(double alpha, double beta, double meanShift)
{
	init(alpha, beta, meanShift);
}

void GammaSampler::init(double alpha, double beta, double meanShift)
{
	if (alpha <= 0.) {cout<<"in gamma sampler, alpha cannot be <= 0"<<endl; exit(1);}
	if (alpha < 1.) alpha += 1.;
	this->alpha = alpha;
	this->oalpha = alpha;
	this->beta = beta;
	this->a1 = alpha-1./3.;
	this->a2 = 1./sqrt(9.*a1);
	this->meanShift = meanShift;
}


//Algorithm for generating random numbers from a gamma distribution, this
//code is modified for c++ using random number generators from NFsim 
//from Numerical Recipies, 3rd edition, page 369-371.
double GammaSampler::gen()
{
	double u,v,x;
	do {
		do {
			x = RANDOM_GAUSSIAN();
			v = 1.+a2*x;
		} while (v<=0);
		v=v*v*v;
		u=RANDOM_CLOSED();
	} while ( (u > (1.-0.331*x*x*x*x)) 
			&& ((log(u)) > ((0.5*x*x)+a1*(1.-v+log(v)))) ); //this line is rarely executed
	if(alpha==oalpha)
		return (a1*v/beta) +meanShift;
	else {
		do {
			u=RANDOM_CLOSED();
		} while (u==0);
		return (pow(u,1./oalpha)*a1*v/beta) +meanShift;
	}
}

double GammaSampler::gen(double minAcceptableValue, double maxAcceptableValue)
{
	double g;
	do {
		 g = this->gen();
	} while(g<minAcceptableValue || g>maxAcceptableValue);
	return g;
}




// To pick a random rotation in 2D, one is able to just pick a random angle with
// which to rotate (on the interval [0,2pi]) and just rotate.  For three dimensions,
// that is not valid and produces a non-uniform distribution.  The general method
// for generating a truly random rotation in n dimensions is slow, but there is
// fortunately a fast method for 3d.  See Shoemake, K, 1992, "Uniform random
// rotations" in Graphics Gems III, Kirk, D., ed. Cambridge, MA: Academic Press,
// pp. 124-132.  A brief explanation of the method can also be found in:
//  Press, Teukolsky, Vetterling, and Flannery.  "Numerical Recipes, The art
//  of Scientific Computing", Third Edition. 2007 Cambridge University Press. (pp. 1130-1)
void ChemotaxisUtil::genUniformRandRotation3d(double rotMatrix[3][3])
{
	//First generate two sets of points ((u0,u1) and (u2,u3) that lie inside the unit circle
	float u0,u1,u2,u3;
	do {
		u0 = 2.*RANDOM_CLOSED()-1.;
		u1 = 2.*RANDOM_CLOSED()-1.;
	} while((u0*u0+u1*u1)>1.);
	
	do {
		u2 = 2.*RANDOM_CLOSED()-1.;
		u3 = 2.*RANDOM_CLOSED()-1.;
	} while((u2*u2+u3*u3)>1.);
	
	
	//Use those points to generate a single random point on a 4d sphere
	float x0 = u0;
	float x1 = u1;
	float temp = sqrt( (1.-u0*u0-u1*u1) / (u2*u2+u3*u3) );
	float x2 = u2 * temp;
	float x3 = u3 * temp;
	
	//Based on that sphere, we can do some simple multiplication to get
	//our rotation matrix.  As a bonus: no need for calling those slow trig functions!
	//rotMatrix[row][col]
	rotMatrix[0][0] = 1.-2.*(x1*x1+x2*x2);
	rotMatrix[0][1] = 2.*(x0*x1-x3*x2);
	rotMatrix[0][2] = 2.*(x0*x2+x3*x1);
	
	rotMatrix[1][0] = 2.*(x0*x1+x3*x2);
	rotMatrix[1][1] = 1.-2.*(x0*x0+x2*x2);
	rotMatrix[1][2] = 2.*(x1*x2-x3*x0);
	
	rotMatrix[2][0] = 2.*(x0*x2-x3*x1);
	rotMatrix[2][1] = 2.*(x1*x2+x3*x0);
	rotMatrix[2][2] = 1.-2.*(x0*x0+x1*x1);
}



void ChemotaxisUtil::genRot3dAboutAxis(double rotMatrix[3][3], double axis[3], double angleInRad)
{
	//We absolutely must normalize first!  The slightest error here throws everything off!
	double mag = sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
	double u1 = axis[0]/=mag; double u2 = axis[1]/=mag; double u3 = axis[2]/=mag;
	double c=cos(angleInRad); double s=sin(angleInRad);
	
	rotMatrix[0][0] = (1.-c)*u1*u1+c;
	rotMatrix[0][1] = (1.-c)*u1*u2-s*u3;
	rotMatrix[0][2] = (1.-c)*u1*u3+s*u2;
	
	rotMatrix[1][0] = (1.-c)*u1*u2+s*u3;
	rotMatrix[1][1] = (1.-c)*u2*u2+c;
	rotMatrix[1][2] = (1.-c)*u2*u3-s*u1;
	
	rotMatrix[2][0] = (1.-c)*u1*u3-s*u2;
	rotMatrix[2][1] = (1.-c)*u2*u3+s*u1;
	rotMatrix[2][2] = (1.-c)*u3*u3+c;
}




void ChemotaxisUtil::genRotFromAngles(double rotMatrix[3][3], double angleX, double angleY, double angleZ)
{
	//Precompute the needed trig functions
	double cosX = cos(angleX), sinX = sin(angleX);
	double cosY = cos(angleY), sinY = sin(angleY);
	double cosZ = cos(angleZ), sinZ = sin(angleZ);
	
	//Generate the rotation matrix
	rotMatrix[0][0] = cosZ*cosY;;
	rotMatrix[0][1] = -sinZ*cosX + cosZ*sinY*sinX;
	rotMatrix[0][2] = sinZ*sinX+cosZ*sinY*cosX;
	rotMatrix[1][0] = sinZ*cosY;
	rotMatrix[1][1] = cosZ*cosX+sinZ*sinY*sinX;
	rotMatrix[1][2] = -cosZ*sinX+sinZ*sinY*cosX;
	rotMatrix[2][0] = -sinY;
	rotMatrix[2][1] = cosY*sinX;
	rotMatrix[2][2] = cosY*cosX;
}






void ChemotaxisUtil::applyRotation(double rotMatrix[3][3], double vec[3])
{
	double newX = rotMatrix[0][0]*vec[0] + rotMatrix[0][1]*vec[1] + rotMatrix[0][2]*vec[2];
	double newY = rotMatrix[1][0]*vec[0] + rotMatrix[1][1]*vec[1] + rotMatrix[1][2]*vec[2];
	double newZ = rotMatrix[2][0]*vec[0] + rotMatrix[2][1]*vec[1] + rotMatrix[2][2]*vec[2];
	vec[0] = newX;
	vec[1] = newY;
	vec[2] = newZ;
}



void ChemotaxisUtil::applyDisplacement(double pos[3], double dir[3], double distance)
{
	pos[0] = pos[0]+distance*dir[0];
	pos[1] = pos[1]+distance*dir[1];
	pos[2] = pos[2]+distance*dir[2];
}



