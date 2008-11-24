#ifndef NFINPUT_HH_
#define NFINPUT_HH_

#include <iostream>
#include <map>
#include <stdio.h>
#include <exception>


#include "../NFcore/NFcore.hh"
#include "../NFoutput/NFoutput.hh"
#include "../NFfunction/NFfunction.hh"
#include "../NFreactions/reactions/reaction.hh"
#include "TinyXML/tinyxml.h"

using namespace NFcore;



//!  Functionality to handle input from XML files or command line arguments
/*!
	This is a very straightforward set of functions and simple classes that
	together parse and handle all the major input into NFsim.  The raw XML
	parsing code comes from the TinyXML project and can be found in its
	unedited form in the TinyXML directory.  Thanks TinyXML!  You perform great
	and you are definitely tiny!
    @author Michael Sneddon
 */
namespace NFinput {

	//! Maintains information about a component of a TemplateMolecule.
	/*!
    	@author Michael Sneddon
	 */
	class component {
		public:
			component(TemplateMolecule *t, string name);
			component(MoleculeType *mt, string name);
			~component();

			TemplateMolecule * t;
			MoleculeType *mt;
			string name;

			string symPermutationName;
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


	//! Reads the Function XML block and adds the Functions to the system.
	/*!
	   	@author Michael Sneddon
	*/
	bool initGlobalFunctions(
		TiXmlElement * pListOfFunctions,
		System * system,
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





	bool FindReactionRuleSymmetry(
			TiXmlElement * pRxnRule,
			System * s,
			map <string,double> &parameter,
			map<string,int> &allowedStates,
			map <string, component> &symComps,
			map <string, component> &symRxnCenter,
			bool verbose);

	bool readPatternForSymmetry(
			TiXmlElement * pListOfMol,
			System * s,
			string patternName,
			map <string, component> &comps,
			map <string, component> &symComps,
			bool verbose);

	bool generateRxnPermutations(vector<map<string,component> > &permutations,
			map<string,component> &symComps,
			map<string,component> &symRxnCenter);



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
			map <string, component> &symMap,
			bool verbose);

	//! Reads a pattern XML block and returns the set of new TemplateMolecule objects.
	/*!
    	@author Michael Sneddon
	 */

		bool readProductPattern(
			TiXmlElement * pListOfMol,
			System * s, map <string,double> &parameter,
			map<string,int> &allowedStates,
			string patternName,
			vector <MoleculeType *> &productMoleculeTypes,
			vector < vector <int> > &stateInformation,
			vector < vector <int> > &bindingSiteInformation,
			bool verbose);



	//! Parses command line arguments from the console nicely.
	/*!
    	@author Michael Sneddon
	 */
	bool parseArguments(int argc, const char *argv[], map<string,string> &argMap);


	int parseAsInt(map<string,string> &argMap, string argName, int defaultValue);
	double parseAsDouble(map<string,string> &argMap, string argName, double defaultValue);




	void walk(System *s);




	bool createComplexOutputDumper(string paramStr, System *s, bool verbose);
	bool createSystemDumper(string paramStr, System *s, bool verbose);

	void parseComplexDump( vector <Outputter> &outputter, string arg );
}



#endif /*NFINPUT_HH_*/
