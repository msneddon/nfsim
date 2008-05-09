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
//	class strCmp {
//	    public:
//	    	bool operator()( const char* s1, const char* s2 ) const {
//	      		return strcmp( s1, s2 ) < 0;
//	    	}
//	};
	
	
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

	
	
	
	System * initializeFromXML(string filename);
	
	bool initParameters(
			TiXmlElement *pListOfParameters, 
			map <string,double> &parameter, 
			bool verbose);
	
	
	bool initMoleculeTypes(
			TiXmlElement * pListOfMoleculeTypes, 
			System * system,
			map<string,int> &allowedStates, 
			bool verbose);
	
	
	bool initStartSpecies(TiXmlElement * pListOfSpecies, System * system, map <string,double> &parameter, map<string,int> &allowedStates, bool verbose);
	
	bool initReactionRules(TiXmlElement * pListOfReactionRules, System * system, map <string,double> &parameter, map<string,int> &allowedStates, bool verbose);
	bool initObservables(TiXmlElement * pListOfObservables, System * system, map <string,double> &parameter, map<string,int> &allowedStates, bool verbose);
	
	TemplateMolecule *readPattern(
			TiXmlElement * pListOfMol, 
			System * s, map <string,double> &parameter,
			map<string,int> &allowedStates,
			string patternName,
			map <string, TemplateMolecule *> &templates, 
			map <string, component> &comps,
			bool verbose);
	
	
	bool addTransformations(TiXmlElement * pListOfProducts, 
			System * s, 
			map <string,double> &parameter, 
			const char *patternName,
			map <const char*, TemplateMolecule *> &reactants,
			ReactionClass *r);
	
	
}



#endif /*NFINPUT_HH_*/
