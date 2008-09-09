//#include "environment.hh"
//
//
//
//
//
//using namespace NG_swim;
//
//
//Environment::Environment() 
//{
//}
//
//double Environment::getLigConc(double xPos, double yPos, double zPos, double time)
//{
////	if(xPos<=0) return 0;
//	return 0;//xPos*10e-8;
//}
//
//void Environment::getStartPosition(int cellNumber, double pos[3])
//{
//	pos[Cell::X] = 0;
//	pos[Cell::Y] = 0;
//	pos[Cell::Z] = 0;
//}
//
//
//
//void Environment::tryToMove(
//					double p[3],
//					double p2[3],
//					double u[3],
//					double up[3]) {
//	//For the empty environment, no constraints on where we can move, so do nothing
//}
//
//
//
//ConstantEnvironment::ConstantEnvironment(double ligandConc) 
//{
//	this->ligandConc=ligandConc;
//}
//
//void ConstantEnvironment::tryToMove(
//		double p[3],
//		double p2[3],
//		double u[3],
//		double up[3]) 
//{
//}
//
//double ConstantEnvironment::getLigConc(double xPos, double yPos, double zPos, double time)
//{
//	return ligandConc;
//}
//
//
//
//
//
//
//
//
//CapillaryEnvironment::CapillaryEnvironment()
//{
////	this->mouthZcoordinate = 1000; //Start 1000um (or 1mm) from base
//	
//	
//	//The following values comes from Futrelle and Berg, Nature 1972	
////	this->initialCapillaryConc = 2.0e-3; //Molar
////	this->diffusionConstant = 890; //um^2 /s (of aspartate)
////	this->radiusOfCapillary = 100; //um, or 0.01cm, 0.1mm
//	
//	
//	//Parameters from Bainer, Park, Cluzel, J Microbio Methods 2003
//	//this->mouthZcoordinate = 5000; //Start 5000um from base (so that -2000um is bottom of the jar
//	//this->initialCapillaryConc = 1e-2;
//	//this->diffusionConstant = 890;
//	//this->radiusOfCapillary = 250; //diameter = 0.5mm, so radius = 0.25mm or 0.025cm or 250um
//	
//	
////	this->initialCapillaryConc = 2e-3;
////	this->diffusionConstant = 890;
////	
////	Zt = 1500;	//Z coordinate of the top of the tube
////	Zb = -200;		//Z coordinate of the base of the tube
////	Zc = 1000;   //Z coordinate of the mouth of the capillary
////	Re = 1000;	//Radius of the entire environment
////	Rc = 120;	//Radius of the capillary
////	
//	Nt[0]=0; Nt[1]=0; Nt[2]=-1;
//	Nb[0]=0; Nb[1]=0; Nb[2]=1;
//	
//	//Hueng-won's capillary...
//	this->initialCapillaryConc = 1e-2;
//	this->diffusionConstant = 890;
//	Zt = 10000;
//	Zc = 7000;
//	Zb = 0;
//	Re = 4000; //2525;
//	Rc = 250;
//	
////	Zt = 500;	//%Z coordinate of the top of the tube
////	Zb = -10;	//%Z coordinate of the base of the tube
////	Zc = 350;   //%Z coordinate of the mouth of the capillary
////	Re = 200;	//%Radius of the entire environment
////	Rc = 80;	//%Radius of the capillary
//}
//
//
//
//
//
//
//CapillaryEnvironment::CapillaryEnvironment(
//					double mouthZcoordinate,
//					double initialCapillaryConc,
//					double diffusionConstant,
//					double radiusOfCapillary
//					)
//{
//	cout<<"wrong capillary constructor called"<<endl;
//	exit(0);
//	
////	this->mouthZcoordinate = mouthZcoordinate;
////	this->initialCapillaryConc = initialCapillaryConc;
////	this->diffusionConstant = diffusionConstant;
////	this->radiusOfCapillary = radiusOfCapillary;
//}
//			
//double CapillaryEnvironment::getLigConc(double xPos, double yPos, double zPos, double time)
//{
//	if(time<=0) return 0;
//	
//	double r0 = sqrt( xPos*xPos + yPos*yPos);
//	if(r0<Rc && zPos>Zc) return initialCapillaryConc;
//	
//	
//	double diffX = xPos;
//	double diffY = yPos;
//	double diffZ = Zc-zPos;
//
//	double r2 = diffX*diffX+diffY*diffY+diffZ*diffZ;
//	double r = sqrt(r2);
//	
//	double c = ( (initialCapillaryConc*Rc*Rc) /
//			(2.*r*(sqrt(NFutil::PI*diffusionConstant*time))) )  *  
//			(exp(-r2/(4.*diffusionConstant*time)))  /  
//			(1.+((3.*Rc*r) /
//			(4.*diffusionConstant*time)));
//	return c;
//}
//
//void CapillaryEnvironment::getStartPosition(int cellNumber, double pos[3])
//{
//	
//	//Choose a random point in the medium, so, choose a random
//	//point in the sphere
//	
//	double x=0,y=0,z=5000;
////	do {
////		
////		z = NFutil::RANDOM_OPEN()*(Zt-Zb)+Zb;
////		x = NFutil::RANDOM_OPEN()*(2*Re)-Re;
////		y = NFutil::RANDOM_OPEN()*(2*Re)-Re;
////		double r = sqrt(x*x+y*y);
////		
////		if(r>Re) continue;
////		if(r<Rc && z>Zc) continue;
////		break;
////	} while(true);
////	
//	
//	pos[Cell::X] = x;
//	pos[Cell::Y] = y;
//	pos[Cell::Z] = z;
//}
//
//
//void CapillaryEnvironment::tryToMove(
//					double p[3],
//					double p2[3],
//					double u[3],
//					double up[3]) {
//	
//	while(!updatePosition(p,p2,u,up));
//}
//
//
//
//
//void norm(double v[3])
//{
//	double mag = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
//	v[0]/=mag; v[1]/=mag; v[2]/=mag;
//}
//
//double dist(double p1[3], double p2[3])
//{
//	double dx = p1[0]-p2[0];
//	double dy = p1[1]-p2[1];
//	double dz = p1[2]-p2[2];
////	cout<<"dx: " <<dx<<" dy: "<<dy<<" dz: "<<dz<<endl;
//	return sqrt(dx*dx+dy*dy+dz*dz);
//}
//
//void translate(double p[3], double v[3], double d)
//{
//	p[0] = p[0]+d*v[0];
//	p[1] = p[1]+d*v[1];
//	p[2] = p[2]+d*v[2];
//	//cout<<"here"<<endl;
//}
//
//void cout_v(const char *name, double v[3])
//{
//	cout<<name<<":\t"<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<endl;
//}
//void cout_v(double v[3])
//{
//	cout<<v[0]<<","<<v[1]<<","<<v[2]<<";"<<endl;
//}
//bool intersectPlane(double p[3], double u[3], double n[3], double d, double q[3])
//{
//	const int x=0,y=1,z=2;
//	double c = u[x]*n[x]+u[y]*n[y]+u[z]*n[z];
//	//cout<<"c:"<<c<<endl;
//	if(fabs(c)<=0.00001)
//	{
//		cout<<"Near parallel"<<endl;
//		return false;
//	}
//	double alpha = (d-(p[x]*n[x]+p[y]*n[y]+p[z]*n[z])) / c;
//	if(alpha<0)
//	{
//		cout<<"No intersection"<<endl;
//		return false;
//	}
//	//cout<<"alpha"<< alpha<<endl;
//	//cout<<u[x]<<" "<<alpha<<endl;
//	q[x] = p[x]+alpha*u[x];
//	q[y] = p[y]+alpha*u[y];
//	q[z] = p[z]+alpha*u[z];
//	//cout_v("intersection:",q);
//	return true;
//}
//
//
//
//bool intersectTube(double p[3], double u[3], double Re, double q[3], double tubeNormal[3])
//{
//	const int x=0,y=1,z=2;
//	
//	//Solution to the quadradic equation intersecting with an infinite
//	//sphere of radius Rc.  The equation:
//	// (p[x]+t*u[x])^2 + (p[y]+t*u[y])^2 = r^2
//	
//	double a = (u[x]*u[x]+u[y]*u[y]);
//	double b = 2.*p[x]*u[x]+2.*p[y]*u[y];
//	double c = p[x]*p[x]+p[y]*p[y]-Re*Re;
//	double rootTerm = b*b-4*a*c;
//	
//	//Make sure the root term is positive, otherwise we don't intersect
//	if(rootTerm<0) return false;
//	
//	//We use this better alternate and equivalent formula for the quadradic equation
//	// [ instead of: t1=(-b+sqrt(b*b-4*a*c))/(2*a);]
//	//that does not make the round off errors like the usual equation may sometimes do
//	//(see any version of Numerical Recipies for the explanation)
//	double sign = 0; if(b>=0) sign=1; else sign=-1;
//	double temp = -0.5*(b+sign*sqrt(rootTerm));
//	if(fabs(a)<=0.0000001 || fabs(temp)<=0.0000001) { cout<<"here\n";return false;  }//return false if we find only one intersection
//	double t1 = temp/a;
//	double t2 = c/temp;
//	
//	//cout<<temp<<endl;
//	//cout<<t1<<"\t---"<<t2<<endl;
//	//cout<<t1a<<"\t---"<<t2a<<endl;
//	
//	//Based on the sign of our parameters t, we can determine if the intersection is in front, behind,
//	//and based on its value, we can determine which came first....
//	if(fabs(t1)<0.00001 && t2>=0) {
//		q[x] = p[x]+t2*u[x];
//		q[y] = p[y]+t2*u[y];
//		q[z] = p[z]+t2*u[z];
//	} else if (fabs(t2)<0.00001 && t1>=0) {
//		q[x] = p[x]+t1*u[x];
//		q[y] = p[y]+t1*u[y];
//		q[z] = p[z]+t1*u[z];
//	} else if(t1<0 && t2<0) {
//		//we should never get here, but this is the case that the tube is completely
//		//behind us  - again, which should never happen as we are in the tube
//		cout<<"\n\n!!! Error in intersectTube: the tube should never be behind us!\n\n";
//		return false;
//	} else if(t1>0) {
//		q[x] = p[x]+t1*u[x];
//		q[y] = p[y]+t1*u[y];
//		q[z] = p[z]+t1*u[z];
//	} else if(t2>0) {
//		q[x] = p[x]+t2*u[x];
//		q[y] = p[y]+t2*u[y];
//		q[z] = p[z]+t2*u[z];
//	} else {
//		cout<<"\n\nSomething unexpected happened when trying to find the intersection with the tube!!\n";
//		cout<<"We should have collided with it, but didn't find the proper intersections!!\n\n";
//		return false;
//	}
//	
//	tubeNormal[x]=-q[x];
//	tubeNormal[y]=-q[y];
//	tubeNormal[z]=0;  //should always be zero, by orientation of the tube
//	norm(tubeNormal);
//	
//	return true;
//}
//
//
//bool intersectCap(double p[3], double u[3], double p2[3], double Rc, double Zc, double q[3], double capNormal[3])
//{
//	const int x=0,y=1,z=2;
//	
//	//Solution to the quadradic equation intersecting with an infinite
//	//sphere of radius Rc.  The equation:
//	// (p[x]+t*u[x])^2 + (p[y]+t*u[y])^2 = r^2
//	
//	double a = (u[x]*u[x]+u[y]*u[y]);
//	double b = 2.*p[x]*u[x]+2.*p[y]*u[y];
//	double c = p[x]*p[x]+p[y]*p[y]-Rc*Rc;
//	double rootTerm = b*b-4*a*c;
//	
//	//Make sure the root term is positive, otherwise we don't intersect
//	if(rootTerm<0) return false;
//	
//	//We use this better alternate and equivalent formula for the quadradic equation
//	//that does not make the round off errors like the usual equation may sometimes do
//	//(see any version of Numerical Recipies for the explanation)
//	double sign = 0; if(b>=0) sign=1; else sign=-1;
//	double temp = -0.5*(b+sign*sqrt(rootTerm));
//	if(fabs(a)<=0.0000001 || fabs(temp)<=0.0000001) return false;  //return false if we find only one intersection
//	double t1 = temp/a;
//	double t2 = c/temp;
//	
////	cout<<t1<<"\t---"<<t2<<endl;
//	
//	//Right away, let us calculate the two intersection points;
//	double q1[3];
//	q1[x] = p[x]+t1*u[x];
//	q1[y] = p[y]+t1*u[y];
//	q1[z] = p[z]+t1*u[z];
//	double q2[3];
//	q2[x] = p[x]+t2*u[x];
//	q2[y] = p[y]+t2*u[y];
//	q2[z] = p[z]+t2*u[z];
//	
//	//if both intersections are further than how far we traveled, then we clearly
//	//are not intersecting on this segment
//	double d1 = dist(p,q1);
//	double d2 = dist(p,q2);
//	double maxD = dist(p,p2);
//	if(d1>maxD && d2>maxD) return false;
//	
//	//if both intersections are behind us, then we must be outside the capillary facing
//	//away from the capillary so any intersection does not matter
//	if(t1<0.0000001 && t2<0.0000001) return false;
//	
//	
//	//if both intersections are in front of us, then we must be outside of the
//	//capillary, so handle that case...
//	if(t1>0.0000001 && t2>0.0000001)
//	{
//		
//		if(q1[z]<Zc && q2[z]<Zc)
//		{
//			// if both intersections lie below the plane, then we didn't hit a thing!
//			return false;
//		}
//		
//		if(q1[z]<Zc && q2[z]>=Zc && p2[z]>Zc)
//		{
//			//Here there are two intersections, both in front of us, so we must be outside the cap
//			//Also, one of the intersections crosses below the plane of the mouth (so it touches nothing)
//			//and the second point where we are moving to lies above the plane of the mouth.  Therefore,
//			//we are outside the cap and are going at an upwards angle into the cap.  That is fine!
//			//We are in!!
//			if(d2<maxD)
//			{
//				q[x]=q2[x]; q[y]=q2[y]; q[z]=q2[z];
//				capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//				return true;
//			}
//				
//			
//		}
//		if(q2[z]<Zc && q1[z]>=Zc && p2[z]>Zc)
//				{
//					//Here there are two intersections, both in front of us, so we must be outside the cap
//					//Also, one of the intersections crosses below the plane of the mouth (so it touches nothing)
//					//and the second point where we are moving to lies above the plane of the mouth.  Therefore,
//					//we are outside the cap and are going at an upwards angle into the cap.  That is fine!
//					//We are in!!
//					if(d1<maxD)
//					{
//						q[x]=q1[x]; q[y]=q1[y]; q[z]=q1[z];
//						capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//						return true;
//					}
//						
//					
//				}
//		
//		//Now, whichever intersection is closer is the one we have to find.  Since
//		//we are outside the cap, we assume the normal of the cap is facing out
//		if(d1<d2 && d1<maxD) {
//			q[x]=q1[x]; q[y]=q1[y]; q[z]=q1[z];
//			capNormal[x]=q[x]; capNormal[y]=q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else if(d2<=d1 && d2<maxD) {
//			q[x]=q2[x]; q[y]=q2[y]; q[z]=q2[z];
//			capNormal[x]=q[x]; capNormal[y]=q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else { 
//			//This is the case that the intersection was further than our next point, so
//			//we don't actually intersect yet
//			return false;
//		}
//	}
//	
//	
//	//If one intersection is clearly in front, and the other behind, then we must be in the
//	//capillary, so handle that case (note that normals are now pointing in!
//	else if(t1>0.0000001 && t2<0.0000001) {
//		if(q1[z]>Zc && d1<maxD) {  //remember to check that the intersection was in the capillary!
//			q[x]=q1[x]; q[y]=q1[y]; q[z]=q1[z];
//			capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else { //Getting here means we are swimming out of the cap... oh my!
//			return false; 
//		}
//	} else if(t2>0.0000001 && t1<0.0000001) {
//		if(q2[z]>Zc && d2<maxD) {  //remember to check that the intersection was in the capillary!
//			q[x]=q2[x]; q[y]=q2[y]; q[z]=q2[z];
//			capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else { //Getting here means we are swimming out of the cap... oh my!
//			return false; 
//		}
//	}
//	
//	
//	//This is where things get interesting...  Imagine that we are a point sitting on the
//	//wall of a capillary (meaning one of the t values ~0).  We are pointing in the direction
//	//of another side of the capillary (the other t is >0).  That must mean we are inside the
//	//capillary (otherwise we would be pointing away from the capillary) and we want to make sure
//	//we intersect with the right one.  This shouldn't really ever happen if our steps are small,
//	//but in general we must account for this.  So let's check this case out...
//	else if(fabs(t1)<0.00001 && t2>0 && d2<maxD)
//	{
//		if(q2[z]>Zc) {  //remember to check that the intersection was in the capillary!
//			q[x]=q2[x]; q[y]=q2[y]; q[z]=q2[z];
//			capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else { //Getting here means we are swimming out of the cap again... oh my!
//			return false; 
//		}
//	}
//	else if(fabs(t2)<0.00001 && t1>0 && d1<maxD)
//	{
//		if(q1[z]>Zc) {  //remember to check that the intersection was in the capillary!
//			q[x]=q1[x]; q[y]=q1[y]; q[z]=q1[z];
//			capNormal[x]=-q[x]; capNormal[y]=-q[y]; capNormal[z]=0; norm(capNormal);
//			return true;
//		} else { //Getting here means we are swimming out of the cap again... oh my!
//			return false; 
//		}
//	}
//	
//	
//	else
//	{
//		cout<<"Didn't catch something!!"<<endl;
//		cout<<t1<<"\t"<<t2<<endl;
//		cout<<d1<<"\t"<<d2<<"\t"<<maxD<<endl;
//		exit(0);
//	}
//	
//	q[0]=0; q[1]=0; q[2]=0;
//	return false;
//}
//
//
//// p: orignal point
//// u: original direction that got us to p2
//// p2: location of next point beyond the reflecting plane
//// q: intersection point of line segment p2-p with plane
//// up: some up vector perpindicular to u, but must be updated
//// n: normal of the plane
//
//
//
///*
// *    ^                        / 
// *    |                       /
// *   (up)                    /
// *    |                     /
// *   [p]----(u)-->--------[q]-------[p2] 
// *                        /
// *                       /
// *              <--(n)--/ 
// *                     / 
// */
//
//void reflect(double p[3], double u[3], double p2[3], double up[3], double q[3], const double n[3], double rv[3])
//{
//	const int x=0,y=1,z=2;
//	
//	//create vector from intersection point to original point
//	double v[3];
//	v[x]=p[x]-q[x]; v[y]=p[y]-q[y]; v[z]=p[z]-q[z];
//	norm(v);  //normalize v;
//	
//	//use equation 2*(v dot n)*n - v to get reflection vector
//	double v_dot_n = v[x]*n[x]+v[y]*n[y]+v[z]*n[z];
//	rv[x] = 2.*(v_dot_n)*n[x]-v[x];
//	rv[y] = 2.*(v_dot_n)*n[y]-v[y];
//	rv[z] = 2.*(v_dot_n)*n[z]-v[z];
//	
//	//set the direction to be the new direction
//	u[x]=rv[x];
//	u[y]=rv[y];
//	u[z]=rv[z];
//	
//	//set the last point to be the intersection point
//	cout_v(p);
//	p[x]=q[x];
//	p[y]=q[y];
//	p[z]=q[z];
//	norm(u);      //normalize the new direction
//	
//	
//	//We have to choose the new up vector (which doesn't matter in this case, it
//	//just has to be orthogonal to the direction.  That is why we just change
//	//one of the coordinates so that we are now in the correct direction
//	if(u[z]>u[x] && u[z]>u[y])
//		up[z]=-(up[x]*u[x]+up[y]*u[y])/u[z];
//	else if(u[y]>u[x])
//		up[y]=-(up[x]*u[x]+up[z]*u[z])/u[y];
//	else
//		up[x]=-(up[y]*u[y]+up[z]*u[z])/u[x];
//	norm(up);
//	
//	
//	
//	//Finally, we have to set the final point (p2).  We must, however, conserve
//	//the displacement, so that we don't go on forever and travel the correct distance
//	
//	// To do that, we first find distance between intersection and p2, then
//	// translate the intersection by the new direction to get the new p2
//	double dis = dist(q,p2);
//	p2[0]=q[0]; p2[1]=q[1]; p2[2]=q[2];
//	translate(p2, rv, dis);
//	
//	cout_v(q);
//}
//
//
//
//
//
//void CapillaryEnvironment::checkPos()
//{
//	//Needed parameters
////	double Zt = 1500;	//Z coordinate of the top of the tube
////	double Nt[3]; Nt[0]=0; Nt[1]=0; Nt[2]=-1;
////	double Zb = 0;		//Z coordinate of the base of the tube
////	double Nb[3]; Nb[0]=0; Nb[1]=0; Nb[2]=1;
////	double Zc = 1000;   //Z coordinate of the mouth of the capillary
////	double Re = 500;	//Radius of the entire environment
////	double Rc = 100;	//Radius of the capillary
//	
//	
//	
//	//Input will be p (last position, that was valid), u(direction of movement from
//	//last position), p2(location that we are trying to move to)
//	
//	double p[3];
//	p[0]=50;
//	p[1]=50;
//	p[2]=50;
//	
//	cout<<"here"<<endl;
//
//	
//	
//	double u[3];
//	u[0]=0.5; u[1]=0.1; u[2]=0.5;
//	norm(u);
//	double up[3];
//	up[0]=0.5;
//	up[1]=0.5;
//	up[2] = -(u[0]*up[0]+u[1]*up[1])/u[2];
//	norm(up);
//	
//	double p2[3];
//	p2[0]=p[0];p2[1]=p[1];p2[2]=p[2];
//	translate(p2,u,40000);
//	
//	const int x=0,y=1,z=2;
//	
//	
//	
//	
//	cout_v("direction",u);
//	cout_v("p1",p);
//	cout_v("p2",p2);
//	
//	
//	//The strategy is to test if the point we are trying to move to is outside the
//	//possible enclosed space as defined by the particular environment.  
//	
//	
//	
//	
//	
//	
//	
//	
//	cout<<endl<<endl<<"points:"<<endl;
//	while(!updatePosition(p,p2,u, up));
//	cout_v("p2",p2);
//	
//
//	
////	//calculate r (distance from center axis of the tube)
////	double r = sqrt( p[x]*p[x] + p[y]*p[y] + p[z]*p[z] );
////	
////	if(r>Re)
////	{
////		
////		
////		
////	}
//	
//	
//	
//	
//	
//	
//	
//	//return c;
//}
//
//
//
//bool CapillaryEnvironment::updatePosition(
//		double p[3],
//		double p2[3],
//		double u[3],
//		double up[3])
//{
//	
//	double nRefPlane[3]; nRefPlane[x]=0; nRefPlane[y]=0; nRefPlane[z]=0;
//	double qFinal[3]; qFinal[x]=0; qFinal[y]=0; qFinal[z]=0;
//	bool intersects = false; double minDist = 0;
//	
////	cout<<"------------------------"<<endl;
////	cout_v("p",p);
////	cout_v("p2",p2);
////	cout_v("u",u);
//	
//	//Calculate the radius from the center of the 
//	double r = sqrt( p2[x]*p2[x] + p2[y]*p2[y]);
//	
//	//This block will always keep us below the top of the chamber
//	if(p2[z]>=Zt && r<Rc)
//	{
//		double q[3]; q[0]=-1; q[1]=-1; q[2]=-1;
//		if(!intersectPlane(p,u,Nt,-Zt,q)) {
//			cout<<"Error!  There should be an intersection with the top!!"<<endl;
//			cout_v("p:",p);cout_v("u:",u);cout_v("p2:",p2);cout_v("up:",up);cout_v("q:",q);
//			exit(1);
//		}
//		if(!intersects) {
//			intersects = true;
//			qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//			nRefPlane[x]=Nt[x]; nRefPlane[y]=Nt[y]; nRefPlane[z]=Nt[z];
//			minDist = dist(p,q);
//		} else {
//			double cDist= dist(p,q);
//			if(cDist<minDist) {
//				qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//				nRefPlane[x]=Nt[x]; nRefPlane[y]=Nt[y]; nRefPlane[z]=Nt[z];
//				minDist = cDist;
//			}
//		}
//	}
//	
//	//This block always keeps us above the bottom of the chamber
//	if(p2[z]<=Zb) {
//		double q[3]; q[0]=-1; q[1]=-1; q[2]=-1;
//		if(!intersectPlane(p,u,Nb,Zb,q)) {
//			cout<<"Error!  There should be an intersection with the base!!"<<endl;
//			cout_v("p:",p);cout_v("u:",u);cout_v("p2:",p2);cout_v("up:",up);cout_v("q:",q);
//			exit(1);
//		}
//		if(!intersects) {
//			intersects = true;
//			qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//			nRefPlane[x]=Nb[x]; nRefPlane[y]=Nb[y]; nRefPlane[z]=Nb[z];
//			minDist = dist(p,q);
//		} else {
//			double cDist= dist(p,q);
//			if(cDist<minDist) {
//				qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//				nRefPlane[x]=Nb[x]; nRefPlane[y]=Nb[y]; nRefPlane[z]=Nb[z];
//				minDist = cDist;
//			}
//		}
//	}
//	
//	//This block always keeps us inside the sides of the cylindrical chamber
//
//	if(r>Re) {
//		double q[3]; q[0]=-1; q[1]=-1; q[2]=-1;
//		double tubeNormal[3];tubeNormal[0]=-1; tubeNormal[1]=-1; tubeNormal[2]=-1;
//		if(!intersectTube(p,u,Re,q,tubeNormal))
//		{
//			cout<<"Error!  There should be an intersection with the tube!!"<<endl;
//			cout_v("p:",p);cout_v("u:",u);cout_v("p2:",p2);cout_v("up:",up);cout_v("q:",q);
//			exit(1);
//		}
//		if(!intersects) 
//		{
//			intersects = true;
//			qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//			nRefPlane[x]=tubeNormal[x]; nRefPlane[y]=tubeNormal[y]; nRefPlane[z]=tubeNormal[z];
//			minDist = dist(p,q);
//		}
//		else
//		{
//			double cDist= dist(p,q);
//			if(cDist<minDist)
//			{
//				qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//				nRefPlane[x]=tubeNormal[x]; nRefPlane[y]=tubeNormal[y]; nRefPlane[z]=tubeNormal[z];
//				minDist = cDist;
//			}
//		}
//	}
//	
//	
//	double r0 = sqrt( p[x]*p[x] + p[y]*p[y]);
//	
//	
//	
//	
//	//We only check for intersection with the capillary if it
//	// a) is not completely contained inside the capillary radius
//	// b) is not completely below the plane of the mouth of the capillary
//	if(!(r0<Rc && r<Rc) && !(p2[z]<Zc && p[z]<Zc))
//	{
//		double q[3]; q[0]=-1; q[1]=-1; q[2]=-1;
//		double capNormal[3];capNormal[0]=-1; capNormal[1]=-1; capNormal[2]=-1;
//		if(intersectCap(p,u,p2,Rc,Zc,q,capNormal))
//		{
//			if(!intersects) 
//			{
//				intersects = true;
//				qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//				nRefPlane[x]=capNormal[x]; nRefPlane[y]=capNormal[y]; nRefPlane[z]=capNormal[z];
//				minDist = dist(p,q);
//			}
//			else
//			{
//				double cDist= dist(p,q);
//				if(cDist<minDist)
//				{
//					qFinal[x]=q[x]; qFinal[y]=q[y]; qFinal[z]=q[z];
//					nRefPlane[x]=capNormal[x]; nRefPlane[y]=capNormal[y]; nRefPlane[z]=capNormal[z];
//					minDist = cDist;
//				}
//			}
//		}
//	}
//
//	if(intersects)
//	{
//		//
//		double rv[3]; rv[0]=0; rv[1]=0;rv[2]=0;
//		reflect(p,u,p2,up,qFinal,nRefPlane,rv);
//		return false;
//	}
//	return true;
//}
//
//
//
//
//
