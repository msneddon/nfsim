////////////////////////////////////////////////////////////////////////////////
//
//    NFsim: The Network Free Stochastic Simulator
//    A software platform for efficient simulation of biochemical reaction
//    systems with a large or infinite state space.
//
//    Copyright (C) 2009,2010,2011,2012
//    Michael W. Sneddon, James R. Faeder, Thierry Emonet
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
//    For more information on NFsim, see http://emonet.biology.yale.edu/nfsim
//
////////////////////////////////////////////////////////////////////////////////



#ifndef NFSIM_HH_
#define NFSIM_HH_

//Include "mpi.h" in Scheduler.h first
#include "NFscheduler/Scheduler.h"

//Include the core files needed to run the simulation
#include "NFcore/NFcore.hh"
#include "NFutil/NFutil.hh"
#include "NFinput/NFinput.hh"
#include "NFreactions/NFreactions.hh"

//Include the specific tests
#include "NFfunction/NFfunction.hh"
//#include "NFtest/compare/compare.hh"
//#include  "NFtest/transformations/transformations.hh"
#include  "NFtest/simple_system/simple_system.hh"
#include  "NFtest/transcription/transcription.hh"
#include  "NFtest/tlbr/tlbr.hh"
#include  "NFtest/agentcell/agentcell.hh"



//! Runs a given System with the specified arguments
/*!
  @author Michael Sneddon
*/
bool runFromArgs(System *s, map<string,string> argMap, bool verbose);


//! Initialize a system from command line flags
/*!
  @author Michael Sneddon
*/
System *initSystemFromFlags(map<string,string> argMap, bool verbose);








#endif /*NFSIM_HH_*/
