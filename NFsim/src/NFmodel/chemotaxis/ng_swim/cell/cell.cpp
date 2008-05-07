#include "cell.hh"

#include <string>
#include <sstream>





template <class T>
inline std::string toString (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}







//motor index is 5, CW observable index is 0
#define MOTOR_INDEX 5
#define TUMBLE_INDEX 0

using namespace std;
using namespace NG_swim;

unsigned int Cell::cellIdCounter = 0;

Cell::Cell(	
		Environment *e,
		double speed,
		double rotDiffusionConstant,
		double startTime,
		const char * outputDirectoryPath,
		NGparam &p
		)
{
	this->cellId = ++Cell::cellIdCounter;
	
	
	////////////////////////////////// Setup output file locations....
	
	//Use some strings to manipulate the path into the filenames we want
	string strPath(outputDirectoryPath);
	
	string cellIdNumber = toString(cellId);
	strPath.append("/c");
	strPath.append(cellIdNumber);
	strPath.append("/");
	
	cout<<strPath<<endl;
	
	string strCellTraj(strPath);
	strCellTraj.append("cellTrajectory.txt");
	const char *cellTrajFileName = strCellTraj.c_str();
	
	string strMolTraj(strPath);
	strMolTraj.append("moleculeTrajectory.txt");
	const char *molTrajFileName = strMolTraj.c_str();
	
	string strClusterKey(strPath);
	strClusterKey.append("clusterKey.txt");
	const char *clusterKeyFileName = strClusterKey.c_str();
	
	string strClusterValues(strPath);
	strClusterValues.append("clusterValues.txt");
	const char *clusterValuesFilename = strClusterValues.c_str();
	
	string strMotorSwitching(strPath);
	strMotorSwitching.append("motorSwitching.txt");
	const char *motorSwitchingFilename = strMotorSwitching.c_str();
	
	
	
	// With shape = 4 scale = 18.32045567939674 location = -4.606176605609249 it
	// fits the distribution plotted in Fig.3, p.501 of Berg & Brown (1972)
	gs = new ChemotaxisUtil::GammaSampler(4.,1./18.32045567939674,-4.606176605609249);
	
	
	////////////////////////////////// Create the system....
	
	// New System
	
	
	

	p.setOutputFileNames(molTrajFileName,clusterKeyFileName,clusterValuesFilename,motorSwitchingFilename);
	
	this->chemotaxisSystem = create("liteSys",p);
		
	/*  Old system... 
	this->chemotaxisSystem = NG::init_NG_system(molTrajFileName, clusterKeyFileName, clusterValuesFilename, motorSwitchingFilename);
	chemotaxisSystem->outputAllObservableNames();
	*/
	
	
	

	//Set output streams for the cell movement
	outputFileStream.open(cellTrajFileName);
	outputFileStream.setf(ios::scientific);

	
	//Set the variables passed in...
	this->env = e;
	this->speed = speed;
	this->rotDiffusionConstant = rotDiffusionConstant;
	currentTime = startTime;//startTime;
	shouldOutputDirection=true;
	
	//Init the cell's position, direction, and orientation
	env->getStartPosition(cellId, pos);
	currentLigConc = env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
	dir[X]=1; dir[Y]=0; dir[Z]=0;
	up[X]=0; up[Y]=1; up[Z]=0;

	//Init the rotation matrix
	rotMat[0][0]=1;   rotMat[0][1]=0;   rotMat[0][2]=0;
	rotMat[1][0]=0;   rotMat[1][1]=1;   rotMat[1][2]=0;
	rotMat[2][0]=0;   rotMat[2][1]=0;   rotMat[2][2]=1;
	
	//Cells should always start with a random orientation
	changeDirRandom();
	
	
	
	
	
//	for(int i=0; i<1000000;i++){
//	ChemotaxisUtil::genUniformRandRotation3d(rotMat);
//	ChemotaxisUtil::applyRotation(rotMat,dir);
//	ChemotaxisUtil::applyRotation(rotMat,up);
//	
//	double mag = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//	dir[X]/=mag; dir[Y]/=mag; dir[Z]/=mag;
//	
//	mag = sqrt(up[X]*up[X]+up[Y]*up[Y]+up[Z]*up[Z]);
//	up[X]/=mag; up[Y]/=mag; up[Z]/=mag;
//	if(!rotMat[0][0]>0 && !rotMat[0][0]>=0 ){ cout<<"last i:"<<i<<endl; break; }
//	}
//	
//	
//	
//	
//	//dir[X] = NFutil::RANDOM_CLOSED()-0.5; dir[Y] = NFutil::RANDOM_CLOSED()-0.5; dir[Z] = NFutil::RANDOM_CLOSED()-0.5;
//	//double mag = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//	//dir[X]/=mag; dir[Y]/=mag; dir[Z]/=mag;
//	
//	
//	
//	
//	cout<<rotMat[0][0]<<"\t"<<rotMat[0][1]<<"\t"<<rotMat[0][2]<<endl;
//	cout<<rotMat[1][0]<<"\t"<<rotMat[1][1]<<"\t"<<rotMat[1][2]<<endl;
//	cout<<rotMat[2][0]<<"\t"<<rotMat[2][1]<<"\t"<<rotMat[2][2]<<endl;
//	
//	cout<<endl;
//	cout<<"dir:"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//	cout<<"up:"<<up[X]<<"\t"<<up[Y]<<"\t"<<up[Z]<<endl;
//	cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//	rotate(rotMat);
//	cout<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//	cout<<endl;
//	
//	
//	
//	cout<<rotMat[0][0]<<"\t"<<rotMat[0][1]<<"\t"<<rotMat[0][2]<<endl;
//	cout<<rotMat[1][0]<<"\t"<<rotMat[1][1]<<"\t"<<rotMat[1][2]<<endl;
//	cout<<rotMat[2][0]<<"\t"<<rotMat[2][1]<<"\t"<<rotMat[2][2]<<endl;
//	
//	
//	cout<<"dir:"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//	cout<<"up:"<<up[X]<<"\t"<<up[Y]<<"\t"<<up[Z]<<endl;
//	cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//	
//	cout<<endl;
//	cout<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//	double mag = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//	cout<<"mag:"<<mag<<endl;
//	rotate(rotMat);
//	cout<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//	mag = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//	cout<<"mag:"<<mag<<endl;
//	cout<<endl;
//	
//	
//	
//	
//	cout<<"------"<<endl;
//	cout<<"dir:"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//		cout<<"up:"<<up[X]<<"\t"<<up[Y]<<"\t"<<up[Z]<<endl;
//		cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//		
//	ChemotaxisUtil::genRot3dAboutAxis(rotMat,dir,3.14);
//	rotate(rotMat);
//	
//	cout<<"dir:"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//			cout<<"up:"<<up[X]<<"\t"<<up[Y]<<"\t"<<up[Z]<<endl;
//			cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//			cout<<"------"<<endl;
//	
//	
//	
//	
//	
//	
//
//	double magd = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//	double magu = sqrt(up[X]*up[X]+up[Y]*up[Y]+up[Z]*up[Z]);
//	double angle = acos((dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z])/(magd*magu));
//	cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//	cout<<"angle:" << angle*180/NFutil::PI<<endl;
//	
//	
	
	
//	double lastDir[3];
//	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
//	for(int i=0;i<10000;i++)
//	{
////		cout<<endl;
//		//changeDirRandom();
//		changeDirDistribution();
//		
//		double magd = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//		double magu = sqrt(lastDir[X]*lastDir[X]+lastDir[Y]*lastDir[Y]+lastDir[Z]*lastDir[Z]);
//		double angle = acos((dir[X]*lastDir[X]+dir[Y]*lastDir[Y]+dir[Z]*lastDir[Z])/(magd*magu));
//
//		//cout<<"lastdir:"<<lastDir[X]<<"\t"<<lastDir[Y]<<"\t"<<lastDir[Z]<<endl;
//		//cout<<"up:"<<up[X]<<"\t"<<up[Y]<<"\t"<<up[Z]<<endl;
//		//cout<<i<<" dir:"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z]<<endl;
//		cout<<angle*180/NFutil::PI<<endl;
//				
//		
////		double magd = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
////		cout<<"mag: "<<magd<<endl;
////		double magu = sqrt(up[X]*up[X]+up[Y]*up[Y]+up[Z]*up[Z]);
////		double angle = acos((dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z])/(magd*magu));
////		cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
////		cout<<"angle:" << angle*180/NFutil::PI<<endl;
//		
//		
//		
//		lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
//	}
//	
//	
//			double magd = sqrt(dir[X]*dir[X]+dir[Y]*dir[Y]+dir[Z]*dir[Z]);
//			double magu = sqrt(up[X]*up[X]+up[Y]*up[Y]+up[Z]*up[Z]);
//			double angle = acos((dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z])/(magd*magu));
//			cout<<"mad d:"<<magd<<" magu:"<<magu<<endl;
//			cout<<"dot product:" << dir[X]*up[X]+dir[Y]*up[Y]+dir[Z]*up[Z]<<endl;
//			cout<<"angle:" << angle*180/NFutil::PI<<endl;
	
	
	

	outputCellHeader();
	
	
//	for(int i=0; i<3000; i++)
//	{
//		this->outputCellValues();
//		swimToNewPosition(0.01);
		//cout<<pos[X]<<","<<pos[Y]<<","<<pos[Z]<<";"<<endl;
//	}
	
	
	
//	equilibriate(5);
//	this->stepTo(20,0.05);
	
	
	
	//Set cell parameters

	
//	rotationalDiffusion = 0.06205; // rad/s
//	// 
//	//  Berg, Random Walks in Biology, pg 84::   0.062 radians^2/sec
//	//  Emonet, (bioinformatics) agent cell, table 1: 0.06205 radians^2/sec
//	//  Andrews Iglesias, optimal noise filtering: 0.16 radians^2 / sec
//	//From berg? 0.16^2//From agent cell 0.06205 rad / sec;;
//	speed = 20; // um/s
//	currentLigandConcentration = this->environment->getLigandConcentrationInMolar(pos[X],pos[Y],pos[Z]);
//	
//	int numOfMotorsTumbling = chemotaxisSystem->getObservableCount(MOTOR_INDEX,TUMBLE_INDEX); 
//	if(numOfMotorsTumbling>0) { currentMovement=TUMBLE; lastMovement=TUMBLE; }
//	else { currentMovement = SWIM; lastMovement=SWIM; }
//	
//	//Init stats variabls
//	stepDisplacement = 0; //um
//	totalDisplacement = 0; //um
//	lastPos[X]=pos[X]; lastPos[Y]=pos[Y]; lastPos[Z]=pos[Z];
//	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
//	startDir[X]=dir[X]; startDir[Y]=dir[Y]; startDir[Z]=dir[Z];
//	angle = 0; //radians
//	angleFromStart=0;
//	
//	
//	runTimeSpent = 0;
}

Cell::~Cell()
{
	outputFileStream.close();
	delete chemotaxisSystem;
	delete gs;
}





//Changes the cell direction uniformly random over all possible new directions
void Cell::changeDirRandom()
{
	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
	lastUp[X]=up[X]; lastUp[Y]=up[Y]; lastUp[Z]=up[Z];
	
	ChemotaxisUtil::genUniformRandRotation3d(rotMat);
	ChemotaxisUtil::applyRotation(rotMat,dir);
	ChemotaxisUtil::applyRotation(rotMat,up);
}


//Changes the cell direction choosen from the given gamma distribution
void Cell::changeDirDistribution()
{
	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
	lastUp[X]=up[X]; lastUp[Y]=up[Y]; lastUp[Z]=up[Z];
	
	//First, rotate the up vector randomly between 0 and 2pi along the direction axis
	//This will not change the direction
	double angleInRad=NFutil::RANDOM_CLOSED()*2.*NFutil::PI;
	ChemotaxisUtil::genRot3dAboutAxis(rotMat,dir,angleInRad);
	ChemotaxisUtil::applyRotation(rotMat,up);
	
	//Second, rotate the direction by the given degree by the set gamma distribution
	angleInRad = gs->gen(0.,180.)*NFutil::PI/180.;
	ChemotaxisUtil::genRot3dAboutAxis(rotMat,up,angleInRad);
	ChemotaxisUtil::applyRotation(rotMat,dir);
}


//The diffusionConstant must be in units rads^2/s
void Cell::changeDirRotDiffusion(double elapsedTime)
{
	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
	lastUp[X]=up[X]; lastUp[Y]=up[Y]; lastUp[Z]=up[Z];
	
	double rotDifTerm = sqrt(elapsedTime * 2. * rotDiffusionConstant);
	double angleX = rotDifTerm * NFutil::RANDOM_GAUSSIAN();
	double angleY = rotDifTerm * NFutil::RANDOM_GAUSSIAN();
	double angleZ = rotDifTerm * NFutil::RANDOM_GAUSSIAN();
	
	ChemotaxisUtil::genRotFromAngles(rotMat,angleX,angleY,angleZ);
	ChemotaxisUtil::applyRotation(rotMat,dir);
	ChemotaxisUtil::applyRotation(rotMat,up);
}

void Cell::moveToNewPosition(double elapsedTime)
{
	double distanceTraveled = speed*elapsedTime;
	lastPos[X]=pos[X]; lastPos[Y]=pos[Y]; lastPos[Z]=pos[Z];
	ChemotaxisUtil::applyDisplacement(pos,dir,distanceTraveled);
}




void Cell::equilibriate(double duration)
{
	cout<<"Cell "<<cellId<<" is equilibriating...  please wait."<<endl;
	
	int n_values = 1;
	double * ligandConc = new double[n_values];
	ligandConc[0] = env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
	
	chemotaxisSystem->updateAllGroupProperty(ligandConc, n_values);
	chemotaxisSystem->equilibriate(duration);
	int numOfMotorsTumbling = chemotaxisSystem->getObservableCount(MOTOR_INDEX,TUMBLE_INDEX); 
	if(numOfMotorsTumbling>0) { currentMovement=TUMBLE; lastMovement=TUMBLE; }
	else { currentMovement = SWIM; lastMovement=SWIM; }
	
	delete [] ligandConc;
}



double Cell::stepTo(double endTime, double dt)
{
//	cout<<"Cell "<<cellId<<" is swimming for "<<endTime<<" seconds ...";  
//	clock_t start,finish;
//	double time;
//	start = clock();
	double simTime = chemotaxisSystem->getCurrentTime();
	double lastTime = simTime;
	double eTime = 0;
	double nextStoppingTime = simTime + dt;
	while(currentTime+dt<endTime)
	{
		//Advance the stochastic simulation, and update time variables
		//Note that there are two time measurements, the simulation time (simTime, lastTime) and the 
		//external system time (currentTime, endTime).  the simulation time always starts at zero when
		//the cell is created.  But a cell can be created in the middle of a simulation
		//in which case its currentTime matches the external time
		simTime = chemotaxisSystem->stepTo(nextStoppingTime);  //advance rxn sim
		chemotaxisSystem->outputAllObservableCounts(); //output rxn stats
		outputCellValues(); //output cell stats
		eTime = simTime-lastTime; //calculate elapsed time since last step
		lastTime = simTime; //remember for next time
		nextStoppingTime += dt; //calculate the next stopping time
		currentTime += eTime; //advance the current system time

		//check to see if the cell is now swimming or tumbling
		int numOfMotorsTumbling = chemotaxisSystem->getObservableCount(MOTOR_INDEX,TUMBLE_INDEX); 
		//cout<<"Num Tumble: "<<numOfMotorsTumbling<<endl;
		if(numOfMotorsTumbling>0) currentMovement=TUMBLE;
		else currentMovement = SWIM;

		//If we start tumbling, but were swimming, we need to choose a new direction
		if(currentMovement==TUMBLE && lastMovement==SWIM)
		{
			
			//This lets us change directions based on the gamma distribution
			this->changeDirDistribution();
			//This lets us change directions completely randomly
			this->changeDirRandom();

			//Update our last position, because we did not move
			lastMovement=TUMBLE;
		}
		//If we are swimming, it doesn't matter what we did last time, just swim already!
		else if(currentMovement==SWIM)
		{
			swimToNewPosition(eTime);
			lastMovement=SWIM;
		}
		
		//Finally, we have to update ligand concentration
		double * value = new double[1];
		currentLigConc= env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
		value[0] = currentLigConc;
		chemotaxisSystem->updateAllGroupProperty(value, 1);
		delete [] value;
	}

//	// Finish and check the run time;
//    finish = clock();
//    time = (double(finish)-double(start))/CLOCKS_PER_SEC;
//    runTimeSpent+=time;
//	cout<<"   done...  Runtime now spent on this cell: "<<runTimeSpent<<" s"<<endl;  
	return currentTime;
}

void Cell::swimToNewPosition(double elapsedTime)
{
	moveToNewPosition(elapsedTime);
	env->tryToMove(this->lastPos,this->pos,this->dir,this->up);
	changeDirRotDiffusion(elapsedTime);
}


void Cell::outputCellHeader()
{
	outputFileStream<<"#\tTIME(s)\tX_POS(um)\tY_POS(um)\tZ_POS(um)\tLigand(M)\tAction(0=swim)";
	if(shouldOutputDirection) 
	{
		outputFileStream<<"\tX_DIR\tY_DIR\tZ_DIR";
	}
	outputFileStream<<endl;
}


void Cell::outputCellValues()
{
	outputFileStream<<"\t"<<currentTime<<"\t"<<pos[X]<<"\t"<<pos[Y]<<"\t"<<pos[Z]<<"\t"<<currentLigConc;
	if(currentMovement==SWIM)outputFileStream<<"\t0";
	else outputFileStream<<"\t1";
	
	if(shouldOutputDirection) 
	{
		outputFileStream<<"\t"<<dir[X]<<"\t"<<dir[Y]<<"\t"<<dir[Z];
	}
	outputFileStream<<endl;
}
