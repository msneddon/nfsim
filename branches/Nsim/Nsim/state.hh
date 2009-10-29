/*
 * state.hh
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#ifndef STATE_HH_
#define STATE_HH_

#include <string>

using namespace std;

namespace Nsim{
	class State{
	public:

		string state;

		State();
		State(string state);
		~State(){}
	};
}

#endif /* STATE_HH_ */
