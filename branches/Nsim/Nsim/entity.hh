/*
 * entity.hh
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#ifndef ENTITY_HH_
#define ENTITY_HH_

#include <vector>

using namespace std;

namespace Nsim{
	class Entity{
	public:

		vector<Entity> entity;
		vector<BindingSite> bindingSite;
		vector<State> state;

		Entity();
		Entity(vector<Entity> entity, vector<BindingSite> bindingSite, vector<State> state);
		~Entity(){};
	};
}

#endif /* ENTITY_HH_ */
