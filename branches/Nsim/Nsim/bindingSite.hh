/*
 * bindingSite.hh
 *
 *  Created on: Oct 28, 2009
 *      Author: Len
 */

#ifndef BINDINGSITE_HH_
#define BINDINGSITE_HH_

#include <string>

using namespace std;

namespace Nsim{
	class BindingSite{
	public:

		string name;

		BindingSite();
		BindingSite(string name);
		~BindingSite(){};
	};
}

#endif /* BINDINGSITE_HH_ */
