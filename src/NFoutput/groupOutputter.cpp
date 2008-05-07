//////////////////////////////////////////////////////////
// groupOutputter.cpp
//
// Defines the methods of the groupOutputter.  The basic
// idea of the groupOutputter is that the user provides a list of
// filenames where each file 
//
//
// Michael Sneddon (michael.sneddon@yale.edu)
//
//////////////////////////////////////////////////////////
#include "NFoutput.hh"


using namespace NFcore;



/**************************************************
 * Create a GroupOutputter with the following properties:
 * 
 *  char * groupKeyFileName - the filename for the output of the group key file which
 *                            will list each group number, name, and the makeup of the
 *                            molecules in that group
 * 
 * vector <TemplateMolecule *> &keyTemplates - a list of template molecules that you would like
 *                                            to count up in each group
 * 
 * vector <char *> &templateNames - a list of the names that coorespond to the template molecules 
 *                                  already given (used in the header of the group key file)
 * 
 * vector <char *> &filenames - the list of files that you want to output at each time step.  These
 *                              files will contain, for each time, all the groups and thier associated
 *                              value from the values vector
 * 
 * vector <unsigned int> &value - the values with which to look up in each group for the given output file
 * 
 * 
 * Therefore, every GroupOutputter will output a key file that gives the groups and thier
 * contents (as specified by the keyTemplates and templateNames vectors.  At each time step,
 * all of the output files, each potentially outputting a different value, will be allowed
 * to dump thier information to the given file.
 * 
 * 
 * Format of output files
 * 
 *   group key file
 * 
 *   #    Group#     Name   Size    TemplateMolecule1   TemplateMolecule2 ....
 * 
 * 
 *   output file number i:
 * 
 *   #    Time       [GroupIndex0]              [GroupIndex1]              .....
 *        x          group0.value[i] @ time=x   group1.value[i] @ time=x   .....
 *        .
 *        .
 * 
 * ************************************************/
GroupOutputter::GroupOutputter(System * s,
	const char * groupName,
	const char * groupKeyFileName, 
	vector <TemplateMolecule *> &keyTemplates,
	vector <const char *> &templateNames,
	vector <const char *> &filenames, 
	vector <unsigned int> &values)
{
	this->groupKeyFileName = groupKeyFileName;
	this->groupName = groupName;
	this->s=s;
	
	//Save the list of template molecules (and thier names) into our vector
	vector <TemplateMolecule *>::iterator tempIter;
	for( tempIter = keyTemplates.begin(); tempIter != keyTemplates.end(); tempIter++ )
		this->keyTemplates.push_back((*tempIter));
	vector <const char *>::iterator tempNameIter;
	for( tempNameIter = templateNames.begin(); tempNameIter != templateNames.end(); tempNameIter++ )
		this->templateNames.push_back((*tempNameIter));
	
	
	//Make sure sizes match up!!
	unsigned int fileCount = filenames.size();
	if(fileCount!=values.size())
	{
		cerr<<"ERROR:  You created a group outputter where the number of filenames and values don't match!"<<endl;
		exit(1);	
	}
	
	
	//open up all of our streams...
	vector <const char *>::iterator filenameIter; int v=0;
	for( filenameIter = filenames.begin(); filenameIter != filenames.end(); filenameIter++, v++ )
	{
		ofstream *o = new ofstream();
		o->open((*filenameIter));
		o->setf(ios::scientific);
		
		outputStreams.push_back(o);
		this->values.push_back(values.at(v));
	}
		
	
	groupCount = s->getGroupCount();
}


// Deconstructor for GroupOutputter
GroupOutputter::~GroupOutputter()
{
	ofstream *o;
	while(outputStreams.size()>0)
	{
		o = outputStreams.back();
		o->flush();
		o->close();
		outputStreams.pop_back();
		delete o;
	}
}

/**
 * This method writes out the list of all groups and the members to a file
 * that we call the key file.  This will let you match up the group id numbers
 * to the memebers of the group.
 */
void GroupOutputter::writeGroupKeyFile()
{
	//Open up the output stream
	ofstream keyStream;
	keyStream.open(groupKeyFileName);
	keyStream.setf(ios::scientific);
	
	//Write the header
	keyStream<<"#\tGroup#\tName\tSize";
	vector <const char *>::iterator tempNameIter;
	for(tempNameIter = templateNames.begin(); tempNameIter != templateNames.end(); tempNameIter++)
		keyStream<<"\t"<<(*tempNameIter);
	keyStream<<endl;
	
	
	//Create a counter vector so we can count what is in each group
	vector <int> counter;
	for( tempIter = keyTemplates.begin(); tempIter != keyTemplates.end(); tempIter++ )
		counter.push_back(0);
		
		
	//We need an iterator to traverse the group
	vector <TemplateMolecule *>::iterator tempIter;
	
	
	//Loop through all the groups and get the scoop on each of them
	int n_groups = s->getGroupCount();
	for(int i=0; i<n_groups; i++)
	{
		//First output the basic stuff
		Group *g = s->getGroup(i);
		const char *gName = g->getName();
		if(strcmp(gName,groupName)!=0) continue;
		
		int groupSize = g->getNumberInGroup();
		keyStream<<"\t"<<i;
		keyStream<<"\t"<<g->getName();
		keyStream<<"\t"<< groupSize;
		
		//See which templates match each molecule in each group
		for(int k=0; k<groupSize; k++)
		{
			int c=0;
			for( tempIter = keyTemplates.begin(); tempIter != keyTemplates.end(); tempIter++, c++) 
				if((*tempIter)->compare(g->getMolecule(k)))
					counter.at(c)++;
		}
		
		
		//Output the template counts
		for(unsigned int c=0; c<counter.size(); c++)
		{
			keyStream<<"\t"<<counter.at(c);
			counter.at(c)=0;
		}
		keyStream<<endl;
	}
	
	// Say goodbye
	keyStream.flush();
	keyStream.close();
}


/**
 * This file outputs the actual data to the file, writing out the given sample time
 */
void GroupOutputter::writeStateToOutputFile(double cSampleTime)
{
	int v=0, n_groups=0;
	for( streamIter = outputStreams.begin(); streamIter != outputStreams.end(); streamIter++,v++)
	{
		(*(*streamIter))<<"\t"<<cSampleTime;
		n_groups = s->getGroupCount();
		if(n_groups!=this->groupCount)
		{
			cerr<<"ERROR!! The number of groups changed!! (reported from groupOutputter)."<<endl;
			exit(1);
		}
		for(int i=0; i<n_groups; i++)
		{
			const char *gName = s->getGroup(i)->getName();
			if(strcmp(gName,groupName)!=0) continue;
			(*(*streamIter))<<"\t"<< s->getGroup(i)->getValue(values.at(v));
		}
		(*(*streamIter))<<endl;
	}
}


/*
 * This method outputs the header information of the group output file
 */
void GroupOutputter::writeOutputFileHeader()
{
	int v=0, n_groups=0;
	for( streamIter = outputStreams.begin(); streamIter != outputStreams.end(); streamIter++,v++)
	{
		(*(*streamIter))<<"#\tTime";
		int n_groups = s->getGroupCount();
		if(n_groups!=this->groupCount)
		{
			cerr<<"ERROR!! The number of groups changed!! (reported from groupOutputter)."<<endl;
			exit(1);
		}
		for(int i=0; i<n_groups; i++)
		{
			const char *gName = s->getGroup(i)->getName();
			if(strcmp(gName,groupName)!=0) continue;
			(*(*streamIter))<<"\t"<<i;
		}
		(*(*streamIter))<<endl;
	}
}


