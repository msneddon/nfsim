#ifndef NFINPUT_HH_
#define NFINPUT_HH_

#include <iostream>
#include <map>
#include <stdio.h> 
#include <exception>


#include "../NFcore/NFcore.hh"
#include "TinyXML/tinyxml.h"

using namespace NFcore;



//!  Functionality to handle input from XML files or command line arguments
/*!
    @author Michael Sneddon
 */
namespace NFinput {

	// class for comparing strings needed for remembering parameters in a map
//	class strCmp {
//	    public:
//	    	bool operator()( const char* s1, const char* s2 ) const {
//	      		return strcmp( s1, s2 ) < 0;
//	    	}
//	};
	

	//! Maintains information about a component of a TemplateMolecule.
	/*!
    	@author Michael Sneddon
	 */
	class component {
		public:
			component(TemplateMolecule *t, int type, string name);
			~component();
			
			TemplateMolecule * t;
			unsigned int type;
			string name;
			
			const static unsigned int BSITE = 0;
			const static unsigned int STATE = 1;
	};

	
	
	//! Maintains information about a component of a TemplateMolecule.
	/*!
    	@author Michael Sneddon
	 */
	System * initializeFromXML(
			string filename,
			bool verbose);
	
	//! Reads the parameter XML block and puts them in the parameter map.
	/*!
    	@author Michael Sneddon
	 */
	bool initParameters(
			TiXmlElement *pListOfParameters, 
			map <string,double> &parameter, 
			bool verbose);
	
	//! Reads the MoleculeType XML block and adds the MoleculeTypes to the system.
	/*!
    	@author Michael Sneddon
	 */
	bool initMoleculeTypes(
			TiXmlElement * pListOfMoleculeTypes, 
			System * system,
			map<string,int> &allowedStates, 
			bool verbose);
	
	//! Reads a Species XML block, creates the molecules and adds them to the system.
	/*!
    	@author Michael Sneddon
	 */
	bool initStartSpecies(
			TiXmlElement * pListOfSpecies, 
			System * system, 
			map <string,double> &parameter, 
			map<string,int> &allowedStates, 
			bool verbose);
	
	//! Reads a reactionRule XML block and adds the rules to the system.
	/*!
    	@author Michael Sneddon
	 */
	bool initReactionRules(
			TiXmlElement * pListOfReactionRules, 
			System * system, 
			map <string,double> &parameter, 
			map<string,int> &allowedStates, 
			bool verbose);
	
	//! Reads an observable XML block and adds the new observables to the system.
	/*!
    	@author Michael Sneddon
	 */
	bool initObservables(
			TiXmlElement * pListOfObservables, 
			System * system, 
			map <string,double> &parameter, 
			map<string,int> &allowedStates, 
			bool verbose);
	
	
	//! Reads a pattern XML block and returns the set of new TemplateMolecule objects.
	/*!
    	@author Michael Sneddon
	 */
	TemplateMolecule *readPattern(
			TiXmlElement * pListOfMol, 
			System * s, map <string,double> &parameter,
			map<string,int> &allowedStates,
			string patternName,
			map <string, TemplateMolecule *> &templates, 
			map <string, component> &comps,
			bool verbose);
	
	
	
	
	//! Parses command line arguments from the console nicely.
	/*!
    	@author Michael Sneddon
	 */
	bool parseArguments(int argc, const char *argv[], map<string,string> &argMap);

}



#endif /*NFINPUT_HH_*/
