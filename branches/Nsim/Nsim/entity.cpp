/*
 * entity.cpp
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#include <iostream>
#include <vector>
#include "entity.hh"

using namespace std;
using namespace Nsim;

Entity::Entity(){
	cout << "Entity constructor called.\n";
}

Entity::Entity(vector<Entity> entity, vector<BindingSite> bindingSite, vector<State> state){
	cout << "Entity constructor called.\n";
	this->entity = entity;
	this->bindingSite = bindingSite;
	this->state = state;
}
