#ifndef NFINPUT_HH_
#define NFINPUT_HH_

#include <iostream>
#include <map>
#include <stdio.h> 
#include <exception>


#include "../NFcore/NFcore.hh"
#include "TinyXML/tinyxml.h"

using namespace NFcore;

namespace NFinput {

	// class for comparing strings needed for remembering parameters in a map
	class strCmp {
	    public:
	    	bool operator()( const char* s1, const char* s2 ) const {
	      		return strcmp( s1, s2 ) < 0;
	    	}
	};
	
	
	System * initializeFromXML(char * filename);
	
	bool initParameters(TiXmlElement *pListOfParameters, map <const char*,double, strCmp> &parameter);
	bool initMoleculeTypes(TiXmlElement * pListOfMoleculeTypes, System * system);
	bool initStartSpecies(TiXmlElement * pListOfSpecies, System * system, map <const char*,double, strCmp> &parameter);
	
	bool initReactionRules(TiXmlElement * pListOfReactionRules, System * system, map <const char*,double, strCmp> &parameter);
	bool initObservables(TiXmlElement * pListOfObservables, System * system, map <const char*,double, strCmp> &parameter);
	
	TemplateMolecule *readPattern(TiXmlElement * pListOfMol, System * s, map <const char*,double, strCmp> &parameter, const char *patternName);
	
}



#endif /*NFINPUT_HH_*/
