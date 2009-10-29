/*
 * bindingSite.cpp
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#include <iostream>
#include "bindingSite.hh"

using namespace std;
using namespace Nsim;

BindingSite::BindingSite(){
	cout << "BindingSite constructor called.\n";
}

BindingSite::BindingSite(string name){
	cout << "BindingSite constructor called.\n";
	this->name = name;
}
