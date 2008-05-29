#ifndef NFOUTPUT_HH_
#define NFOUTPUT_HH_

#include "../NFcore/NFcore.hh"

#include <vector>




using namespace std;

//! Manages some of the outputting functionality.
/*!
	@author Michael Sneddon
 */
namespace NFcore
{
	//Forward declarations of classes from NFcore.hh

	class System;
	class Molecule;		
	class TemplateMolecule; 
	class Group;
	
	class GroupOutputter
	{
		public:
			GroupOutputter(System * s, 
					string groupName, 
					string groupKeyFileName, 
					vector <TemplateMolecule *> &keyTemplates, 
					vector <char *> &templateNames, 
					vector <char *> &filenames, 
					vector <unsigned int> &values);
			~GroupOutputter();
			
			void writeGroupKeyFile();
			void writeOutputFileHeader();
			void writeStateToOutputFile(double cSampleTime);
			
			
		protected:
			
			System *s;
			string groupKeyFileName;
			string groupName;
			
			vector <TemplateMolecule *> keyTemplates;
			vector <char *> templateNames;
			vector <ofstream *> outputStreams;
			vector <int> values;
			
			//For this to all work, the number of groups must be static (they can be empty) but
			//must be static.  We will enforce this with this variable
			int groupCount;
			
		private:
			vector <Molecule *>::iterator molIter;
			vector <ofstream *>::iterator streamIter;
			vector <TemplateMolecule *>::iterator tempIter;
			
			
	};
	
	

}









#endif /*NFOUTPUT_HH_*/
