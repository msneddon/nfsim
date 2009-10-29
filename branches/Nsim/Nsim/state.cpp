/*
 * state.cpp
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#include <iostream>
#include "state.hh"

using namespace std;
using namespace Nsim;

State::State(){
	cout << "State constructor called.\n";
}

State::State(string state){
	cout << "State constructor called.\n";
	this->state = state;
}
