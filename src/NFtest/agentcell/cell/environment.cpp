#include "environment.hh"






Environment::Environment()
{
}

double Environment::getLigConc(double xPos, double yPos, double zPos, double time)
{
//	if(xPos<=0) return 0;
	return 0;//xPos*10e-8;
}

void Environment::getStartPosition(int cellNumber, double pos[3])
{
	pos[AgentCell::X] = 0;
	pos[AgentCell::Y] = 0;
	pos[AgentCell::Z] = 0;
}



void Environment::tryToMove(
					double p[3],
					double p2[3],
					double u[3],
					double up[3]) {
	//For the empty environment, no constraints on where we can move, so do nothing
}



ConstantEnvironment::ConstantEnvironment(double ligandConc)
{
	this->ligandConc=ligandConc;
}

void ConstantEnvironment::tryToMove(
		double p[3],
		double p2[3],
		double u[3],
		double up[3])
{
}

double ConstantEnvironment::getLigConc(double xPos, double yPos, double zPos, double time)
{
	return ligandConc;
}



void LinearEnvironment::getStartPosition(int cellNumber, double pos[3]){
	               pos[AgentCell::X] = 0;
	               pos[AgentCell::Y] = 0;
	               pos[AgentCell::Z] = 500;
               }


