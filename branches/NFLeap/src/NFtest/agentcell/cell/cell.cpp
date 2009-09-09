#include "cell.hh"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>






using namespace std;

unsigned int AgentCell::cellIdCounter = 0;


const double AgentCell::DefaultGammaAlpha = 4.;
const double AgentCell::DefaultGammaBeta = 1./18.32045567939674;
const double AgentCell::DefaultGammaOffset = -4.606176605609249;

const double AgentCell::DefaultSpeed = 20; //uM per second
const double AgentCell::DefaultRotDifConst = 0.06205; // rad^2/s


AgentCell::AgentCell(
		System *s,
		Environment *e,
		double speed,
		double rotDiffusionConstant,
		double cheYpThreshold,
		double startTime,
		string outputDirectoryPath
		)
{
	this->system = s;
	this->cellId = ++AgentCell::cellIdCounter;


	//Set the variables passed in...
	this->env = e;
	this->speed = speed;
	this->rotDiffusionConstant = rotDiffusionConstant;
	currentTime = startTime;//startTime;

	//Init the cell's position, direction, and orientation
	env->getStartPosition(cellId, pos);
	currentLigConc = env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
	dir[X]=1; dir[Y]=0; dir[Z]=0;
	up[X]=0; up[Y]=1; up[Z]=0;

	//Init the rotation matrix
	rotMat[0][0]=1;   rotMat[0][1]=0;   rotMat[0][2]=0;
	rotMat[1][0]=0;   rotMat[1][1]=1;   rotMat[1][2]=0;
	rotMat[2][0]=0;   rotMat[2][1]=0;   rotMat[2][2]=1;


	////////////////////////////////// Setup output file locations....

	string strCellTraj(outputDirectoryPath);
	strCellTraj.append("cellTraj");
	const char *cellTrajFileName = strCellTraj.c_str();


	// With shape = 4 scale = 18.32045567939674 location = -4.606176605609249 it
	// fits the distribution plotted in Fig.3, p.501 of Berg & Brown (1972)
	gs = new ChemotaxisUtil::GammaSampler(4.,1./18.32045567939674,-4.606176605609249);



	//Set output streams for the cell movement (Now in binary!!! )
	fileName = cellTrajFileName;
	outputFileStream.open((fileName+".dat").c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
	if(!outputFileStream.is_open()) {
		cout<<"could not open stream to output file."<<endl;
		cout<<(fileName+".dat").c_str()<<endl;
		exit(1);
	}
	outputCellHeader();

	//outputFileStream.open(cellTrajFileName);
	//outputFileStream.setf(ios::scientific);


	this->cheYhistory=0;
	this->cheYhisPos=0;
	this->meanCheYp=0;
	this->cheYhistorySum=0;

	this->cheYpThreshold = (int)round(cheYpThreshold);// usually is about ~ 1795;
	this->motorState = this->CCW;
	this->flagellaState = this->APART;
	this->lastFlagellaState = this->APART;
	this->apartDuration = 0;
	this->droppingTumble = false;


	this->boxcarTimeWidth = 0.3;
	this->cheYhistorySize = 0;// = (int)round(boxcarTimeWidth/dt);


	//Cells should always start with a random orientation
	changeDirRandom();
	binaryFileOutputCounter = 0;
}

AgentCell::~AgentCell()
{
	outputCellHeader();
	outputFileStream.close();
	motorFileStream.close();
	delete gs;
}





//Changes the cell direction uniformly random over all possible new directions
void AgentCell::changeDirRandom()
{
	lastDir[X]=dir[X]; lastDir[Y]=dir[Y]; lastDir[Z]=dir[Z];
	lastUp[X]=up[X]; lastUp[Y]=up[Y]; lastUp[Z]=up[Z];

	ChemotaxisUtil::genUniformRandRotation3d(rotMat);
	ChemotaxisUtil::applyRotation(rotMat,dir);
	ChemotaxisUtil::applyRotation(rotMat,up);
}


//Changes the cell direction choosen from the given gamma distribution
void AgentCell::changeDirDistribution()
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
void AgentCell::changeDirRotDiffusion(double elapsedTime)
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

void AgentCell::moveToNewPosition(double elapsedTime)
{
	double distanceTraveled = speed*elapsedTime;
	lastPos[X]=pos[X]; lastPos[Y]=pos[Y]; lastPos[Z]=pos[Z];
	ChemotaxisUtil::applyDisplacement(pos,dir,distanceTraveled);
}




void AgentCell::equilibriate(double duration, double dt)
{
	int progressStep = round(((duration)/dt)/40);
	int currentStep = 0;
	cout<<"start -------------------------------------- end"<<endl;
	cout<<"      ";


	//boxcarTimeWidth = 0.3;
	this->cheYhistorySize = (int)round(boxcarTimeWidth/dt);

	//since we are equilibrating, we have to reinitialize the
	//cheYp history.
	this->cheYhistory =new double [cheYhistorySize];
	for(int i=0; i<cheYhistorySize; i++) {
		cheYhistory[i]=(double)this->system->getObservableByName("Yp")->getCount();
		cheYhistorySum += cheYhistory[i];
	}
	this->cheYhisPos = 0;


	//cout<<"Current Yp: "<<this->system->getObservableByName("Yp")->getCount() <<endl;
	//cout<<"Current Yp Sum: "<<cheYhistorySum<<endl;
	//cout<<"Mean Yp: "<<cheYhistorySum/(double)cheYhistorySize;


	//Get and set the ligand concentration at the starting position...
	this->currentLigConc = env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
	double L = currentLigConc;
	this->system->setParameter("L",L);
	this->system->updateSystemWithNewParameters();


	for(double time=0; time<duration; time+=dt) {

		//Run the system for the given time
		system->equilibrate(dt);

		//  Update the cheY history array
		cheYhistorySum-=cheYhistory[cheYhisPos];
		cheYhistory[cheYhisPos]=system->getObservableByName("Yp")->getCount();
		cheYhistorySum+=cheYhistory[cheYhisPos];
		//cout<<"this cheYp: "<<cheYhistory[cheYhisPos]<<endl;
		cheYhisPos++;
		if(cheYhisPos>=cheYhistorySize) cheYhisPos=0;
		meanCheYp = cheYhistorySum/(double)cheYhistorySize;
		//cout<<cheYhisPos<<" : "<<meanCheYp<<endl;


		//Use the history to determine motor state
		if(meanCheYp>cheYpThreshold) {
			motorState = this->CW;
		} else {
			motorState = this->CCW;
		}

		this->lastFlagellaState=flagellaState;
		//Use motor state to determine flagella state, and drop 20% of the tumbles
		if(motorState == this->CW) {
			flagellaState = this->APART;
		} else {
			flagellaState = this->BUNDLED;
		}
		if(flagellaState==this->APART) {
			if(apartDuration<(dt/2)) {
				if(NFutil::RANDOM_CLOSED()<0.2) droppingTumble = true;
			}
			if(droppingTumble) {
				flagellaState=this->BUNDLED;
			}
			apartDuration+=dt;
		} else { apartDuration=0; };



		//cout<<"Flagella State: "<< flagellaState<<endl;

		currentStep++;
		if(currentStep%progressStep==0) {
			cout<<"*"; cout.flush();
		}
	}

	cout<<endl;
}



double AgentCell::stepTo(double endTime, double dt)
{
	int progressStep = round(((endTime - currentTime)/dt)/38);
	int currentStep = 0;
	cout<<"start -------------------------------------- end"<<endl;
	cout<<"      ";

	int CWrot = 0;
	int CCWrot = 0;
	int swim = 0;
	int tumble = 0;
	while(currentTime<=endTime)
	{
		//First we output the current information
		system->outputAllObservableCounts(currentTime);
		this->outputCellValues();

		//Next we simulate for the time step
		system->stepTo(currentTime+dt);

		//  Update the cheY history array
		cheYhistorySum-=cheYhistory[cheYhisPos];
		cheYhistory[cheYhisPos]=system->getObservableByName("Yp")->getCount();
		cheYhistorySum+=cheYhistory[cheYhisPos];
		cheYhisPos++;
		if(cheYhisPos>=cheYhistorySize) cheYhisPos=0;
		meanCheYp = cheYhistorySum/(double)cheYhistorySize;


		//Use the history to determine motor state
		if(meanCheYp>cheYpThreshold) {
			motorState = this->CW;
			CWrot++;
		} else {
			motorState = this->CCW;
			CCWrot++;
		}


		//Use motor state to determine flagella state, and drop 20% of the tumbles
		this->lastFlagellaState=flagellaState;
		if(motorState == this->CW) {
			flagellaState = this->APART;
			tumble++;
		} else {
			flagellaState = this->BUNDLED;
			swim++;
		}
		if(flagellaState==this->APART) {
			if(apartDuration<(dt/2)) {
				if(NFutil::RANDOM_CLOSED()<0.2) droppingTumble = true;
			}
			if(droppingTumble) {
				flagellaState=this->BUNDLED;
			}
			apartDuration+=dt;
		} else { apartDuration=0; droppingTumble = false; };


		///////// Move the cell
		//If we start tumbling, but were swimming, we need to choose a new direction
		if(flagellaState==APART && lastFlagellaState==BUNDLED)
		{
		  //This lets us change directions based on the gamma distribution
		  this->changeDirDistribution();
		  //This lets us change directions completely randomly
		  //this->changeDirRandom();
		}
		//If we are swimming, it doesn't matter what we did last time, just swim already!
		else if(flagellaState==BUNDLED)
		{
			swimToNewPosition(dt);
		}

		//and update the time
		currentTime+=dt;


		//Now update the current ligand concentration at the new location
		double L = env->getLigConc(pos[X],pos[Y],pos[Z],currentTime);
		if(L!=currentLigConc)
		{
			currentLigConc=L;
			this->system->setParameter("L",currentLigConc);
			this->system->updateSystemWithNewParameters();
		}

		currentStep++;
		if(currentStep%progressStep==0) {
			cout<<"*"; cout.flush();
		}
	}

	cout<<endl;
	cout<<"Final CW bias: "<<((float)CWrot/((float)CCWrot+(float)CWrot))<<endl;
	cout<<"Swimming bias: "<<(float)swim/((float)swim+(float)tumble)<<endl;
	return currentTime;
}

void AgentCell::swimToNewPosition(double elapsedTime)
{
	moveToNewPosition(elapsedTime);
	env->tryToMove(this->lastPos,this->pos,this->dir,this->up);
	changeDirRotDiffusion(elapsedTime);
}


void AgentCell::outputCellHeader()
{
	//new, binary technique
	ofstream o((fileName+".hd").c_str());
	o<<"#\trows\tTIME(s)\tX_POS(um)\tY_POS(um)\tZ_POS(um)\tLigand(M)\tMotor("<<this->CCW<<"=CCW)\tFlagella("<<this->BUNDLED<<"=BUNDLED)\tMeanCheYp";
	o<<"\n";
	o<<"\t"<<binaryFileOutputCounter<<"\t0\t0\t0\t0\t0\t0\t0\t0\t0";
	o.close();
}


void AgentCell::outputCellValues()
{
	binaryFileOutputCounter++;

	outputFileStream.write((char *)&currentTime, sizeof(double));
	outputFileStream.write((char *)&pos[X], sizeof(double));
	outputFileStream.write((char *)&pos[Y], sizeof(double));
	outputFileStream.write((char *)&pos[Z], sizeof(double));
	outputFileStream.write((char *)&currentLigConc, sizeof(double));

	double v = this->motorState;
	outputFileStream.write((char *)&v, sizeof(double));

	v = (double)this->flagellaState;
	outputFileStream.write((char *)&v, sizeof(double));

	v = (double)this->meanCheYp;
	outputFileStream.write((char *)&v, sizeof(double));

}


