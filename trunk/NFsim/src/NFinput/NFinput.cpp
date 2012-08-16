#include "NFinput.hh"



#include <algorithm>


using namespace NFinput;
using namespace std;



component::component(TemplateMolecule *t, string name)
{
	mt=0;
	uniqueId="";
	symPermutationName="";
	numOfBondsLabel="";
	stateConstraintLabel="";


	this->t=t;
	this->name = name;
}

component::component(MoleculeType *mt, string name)
{
	t=0;
	uniqueId="";
	symPermutationName="";
	numOfBondsLabel="";
	stateConstraintLabel="";

	this->mt=mt;
	this->name = name;
}

component::~component()
{
	t=0;
	mt=0;
}



System * NFinput::initializeFromXML(
		string filename,
		bool blockSameComplexBinding,
		int globalMoleculeLimit,
		bool verbose,
		int &suggestedTraversalLimit,
		bool evaluateComplexScopedLocalFunctions )
{
	if(!verbose) cout<<"reading xml file ("+filename+")  \n\t[";
	if(verbose) cout<<"\tTrying to read xml model specification file: \t\n'"<<filename<<"'"<<endl;


	TiXmlDocument doc(filename.c_str());
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		if(verbose) cout<<"\t\tread was successful... beginning parse..."<<endl<<endl;

		//First declare our system
		System *s;

		//Read in the root node, which should give us the system's name
		TiXmlHandle hDoc(&doc);
		TiXmlElement *pModel = hDoc.FirstChildElement().Node()->FirstChildElement("model");
		if(!pModel) { cout<<"\tNo 'model' tag found.  Quitting."; return NULL; }

		//Make sure the basics are there
		string modelName;
		if(!pModel->Attribute("id"))  {
			if(!blockSameComplexBinding) s=new System("nameless",false,globalMoleculeLimit);
			else s=new System("nameless",true,globalMoleculeLimit);
			if(verbose) cout<<"\tNo System name given, so I'm calling your system: "<<s->getName()<<endl;
		}
		else  {
			modelName=pModel->Attribute("id");
			//We have to add complex bookkeeping if we are blocking same complex binding
			if(!blockSameComplexBinding) s=new System(modelName,false,globalMoleculeLimit);
			else s=new System(modelName,true,globalMoleculeLimit);
			if(verbose) cout<<"\tCreating system: "<<s->getName()<<endl;
		}

		// set evaluation of complex-scoped local functions (true or false)
		s->setEvaluateComplexScopedLocalFunctions(evaluateComplexScopedLocalFunctions);

		//Read the key lists needed for the simulation and make sure they exist...
		TiXmlElement *pListOfParameters = pModel->FirstChildElement("ListOfParameters");
		if(!pListOfParameters) { cout<<"\tNo 'ListOfParameters' tag found.  Quitting."; delete s; return NULL; }
		TiXmlElement *pListOfFunctions = pModel->FirstChildElement("ListOfFunctions");
		//(we do not enforce that functions must exist... yet)  if(!pListOfFunctions) { cout<<"\tNo 'ListOfParameters' tag found.  Quitting."; delete s; return NULL; }
		TiXmlElement *pListOfMoleculeTypes = pListOfParameters->NextSiblingElement("ListOfMoleculeTypes");
		if(!pListOfMoleculeTypes) { cout<<"\tNo 'ListOfMoleculeTypes' tag found.  Quitting."; delete s; return NULL; }
		TiXmlElement *pListOfSpecies = pListOfMoleculeTypes->NextSiblingElement("ListOfSpecies");
		if(!pListOfSpecies) { cout<<"\tNo 'ListOfSpecies' tag found.  Quitting."; delete s; return NULL; }
		TiXmlElement *pListOfReactionRules = pListOfSpecies->NextSiblingElement("ListOfReactionRules");
		if(!pListOfReactionRules) { cout<<"\tNo 'ListOfReactionRules' tag found.  Quitting."; delete s; return NULL; }
		TiXmlElement *pListOfObservables = pListOfReactionRules->NextSiblingElement("ListOfObservables");
		if(!pListOfObservables) { cout<<"\tNo 'ListOfObservables' tag found.  Quitting."; delete s; return NULL; }


		//Now retrieve the parameters, so they are easy to look up in the future
		//and save the parameters in a map we call parameter
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading parameter list..."<<endl;
		map<string, double> parameter;
		if(!initParameters(pListOfParameters, s, parameter, verbose))
		{
			cout<<"\n\nI failed at parsing your Parameters.  Check standard error for a report."<<endl;
			if(s!=NULL) delete s;
			return NULL;
		}

		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of MoleculeTypes..."<<endl;
		map<string,int> allowedStates;
		if(!initMoleculeTypes(pListOfMoleculeTypes, s, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your MoleculeTypes.  Check standard error for a report."<<endl;
			if(s!=NULL) delete s;
			return NULL;
		}


		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Species..."<<endl;
		if(!initStartSpecies(pListOfSpecies, s, parameter, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your species.  Check standard error for a report."<<endl;
			if(s!=NULL) delete s;
			return NULL;
		}


		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Observables..."<<endl;
		if(!initObservables(pListOfObservables, s, parameter, allowedStates, verbose, suggestedTraversalLimit))
		{
			cout<<"\n\nI failed at parsing your observables.  Check standard error for a report."<<endl;
			if(s!=NULL) delete s;
			return NULL;
		}



		if(!verbose) cout<<"-";
		else if(pListOfFunctions) cout<<"\n\tReading list of Functions..."<<endl;
		if(pListOfFunctions)
		{
			if(!initFunctions(pListOfFunctions, s, parameter, pListOfObservables,allowedStates,verbose)) {
				cout<<"\n\nI failed at parsing your Global Functions.  Check standard error for a report."<<endl;
				if(s!=NULL) delete s;
				return NULL;
			}
		}



		//We have to read reactionRules AFTER observables because sometimes reactions
		//might depend on some observable...
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Reaction Rules..."<<endl;

		if(!initReactionRules(pListOfReactionRules, s, parameter, allowedStates, blockSameComplexBinding, verbose, suggestedTraversalLimit))
		{
			cout<<"\n\nI failed at parsing your reaction rules.  Check standard error for a report."<<endl;
			if(s!=NULL) delete s;
			return NULL;
		}

		/////////////////////////////////////////
		// Parse is finally over!  Now we just have to take care of some final details.

		//Finish up the output message
		if(!verbose) cout<<"-]\n";

		//We no longer prepare the simulation here!  You have to do it yourself

		return s;
	}
	else
	{
		cout<<"\nError reading the file.  I could not find / open it, or it is not valid xml."<<endl;
	}


	return 0;
}




bool NFinput::initParameters(TiXmlElement *pListOfParameters, System *s, map <string,double> &parameter, bool verbose)
{
	try {
		//Loop through all the parameter elements
		TiXmlElement *pParamElement;
		for ( pParamElement = pListOfParameters->FirstChildElement("Parameter");
				pParamElement != 0; pParamElement = pParamElement->NextSiblingElement("Parameter"))
		{

			//Try to extract out the parameter name
			string paramName;
			if(!pParamElement->Attribute("id")) {
				cerr<<"\t\t!!A Parameter is undefined! It is missing the 'id' attribute!  Quitting.\n";
				return false;
			} else {
				paramName = pParamElement->Attribute("id");
			}

			//Try to parse the value of the parameter
			string paramValue;
			if(!pParamElement->Attribute("value")) {
				cerr<<"\t\t!!A Parameter '"<<paramName<<"' does not have the 'value' attribute defined! Quitting.\n";
				return false;
			} else {
				paramValue = pParamElement->Attribute("value");
			}

			//See if we can get the parameter value into a double
			double d = 0;
			try {
				d = NFutil::convertToDouble(paramValue);
			} catch (std::runtime_error e) {
				cerr<<"I couldn't understand your parameter: "<<paramName<<endl;
				cerr<<e.what()<<endl;
				return false;
			}

			// Make sure we haven't tried to add this parameter before
			if(parameter.find(paramName)!=parameter.end()) {
				cerr<<"You tried to define the parameter '"+paramName+"' more than once.  I think"<<endl;
				cerr<<"you made a mistake.  I'll forgive you, but I'm quitting now."<<endl;
				return false;
			}

			//Save the parameter into our parameter map
			parameter[paramName]=d;
			s->addParameter(paramName,d);
			if(verbose) cout<<"\t\t Identified parameter:\t"<<paramName<<"\tValue:"<<d<<endl;
		}
		return true;
	} catch (...) {
		cerr<<"Undefined exception thrown while parsing the parameters."<<endl;
		return false;
	}
	return false;
}













/**
 *
 * The strategy is to look at one MoleculeType at a time, make sure that moleculeType contains
 * valid information, then, and only then, create it.
 *
 *
 */
bool NFinput::initMoleculeTypes(
		TiXmlElement * pListOfMoleculeTypes,
		System * s,
		map<string,int> &allowedStates,
		bool verbose)
{
	try {
		vector <string> compLabels;
		vector <string> defaultCompState;
		vector <vector <string> > possibleComponentStates;
		vector <vector <string> > identicalComponents;
		vector <bool> isIntegerComponent;
		//Organized as each vector in this vector contains the set of names of components that are identical

		vector <int> firstSymSiteToAppend;  //If we find a symmetric component, we should go
		                                    //back and name append a '1' to the name.

		//Loop through the MoleculeType tags...
		TiXmlElement *pMoTypeEl;
		for ( pMoTypeEl = pListOfMoleculeTypes->FirstChildElement("MoleculeType"); pMoTypeEl != 0; pMoTypeEl = pMoTypeEl->NextSiblingElement("MoleculeType"))
		{
			//Check if MoleculeType tag has a name...
			if(!pMoTypeEl->Attribute("id")) {
				cerr<<"!!!Error:  MoleculeType tag must contain the id attribute.  Quitting."<<endl;
				return false;
			}

			//Read in the MoleculeType Name
			string typeName = pMoTypeEl->Attribute("id");

			//Make sure the name isn't null & output a message if needed
			string lowerCaseTypeName = "";
			lowerCaseTypeName.reserve(typeName.size());
			std::transform(typeName.begin(),typeName.end(),lowerCaseTypeName.begin(),::tolower);
			if(typeName.compare("null")==0 || typeName.compare("trash")==0){
				if(verbose) cout<<"\t\tSkipping Moleculetype of name: '" + typeName + "'"<<endl;
				continue;
			}
			if(verbose) cout<<"\t\tReading and Creating MoleculeType: "+typeName+"(";


			//Check for population flag  --justin
			bool isPopulation = false;
			if ( pMoTypeEl->Attribute("population") )
			{
				string popFlag = pMoTypeEl->Attribute("population");
				if( popFlag.compare("1") == 0 )
				{
					isPopulation = true;
					if(verbose) cout << "\t\tTreating molecule type '" << typeName << "' as a population" << endl;
				}
			}

			//Get the list of components in the moleculeType
			TiXmlElement *pListOfComp = pMoTypeEl->FirstChildElement("ListOfComponentTypes");
			if(pListOfComp)
			{
				if ( isPopulation )
				{   // Populations aren't allowed to have components!  --Justin
					cerr<<"!!!Error:  Population type " << typeName << " is not allowed to have components! Quitting." << endl;
					return false;
				}

				//Loop through the list of components
				TiXmlElement *pComp;
				for ( pComp = pListOfComp->FirstChildElement("ComponentType"); pComp != 0; pComp = pComp->NextSiblingElement("ComponentType"))
				{
					//Check again for errors by making sure the component has a name
					if(!pComp->Attribute("id")) {
						if(verbose) cout<<" ?? ...\n";
						cerr<<"!!!Error:  ComponentType tag in MoleculeType: '" + typeName + "' must contain the id attribute.  Quitting."<<endl;
						return false;
					}

					//Read in the component name and output its value
					string compName = pComp->Attribute("id");
					if(verbose) {
						if(compLabels.size()!=0) cout<<","+compName; else cout<<compName;
					}
					//cout<<"\n analyzing: "<<compName<<endl;

					//First check if the component Name already exists, if so, we gotta do more!
					//This means that one of the sites are symmetric, so we must handle it correctly
					int pos=0;
					for(vector<string>::iterator it = compLabels.begin(); it != compLabels.end(); it++,pos++ ) {

						//cout<<"comparing: "<<(*it)<<" to "<<compName<<endl;
						if((*it)==compName) {
							bool shouldAdd = true;
							for(unsigned int k=0; k<firstSymSiteToAppend.size(); k++) {
								if(firstSymSiteToAppend.at(k)==pos) { shouldAdd=false; break; }
							}
							if(shouldAdd) firstSymSiteToAppend.push_back(pos);


							string newCompName = compName;
							string num="0"; bool matchedSiteName=false;
							for(unsigned int is=0; is<identicalComponents.size(); is++)
							{
								if(identicalComponents.at(is).at(0)==(compName+"1"))
								{
									unsigned int lastIndex = identicalComponents.at(is).size();
									std::stringstream lastIndexStream; lastIndexStream << lastIndex+1;
									num = lastIndexStream.str();
									//cout<<"identified num = "<<num<<endl;
									newCompName = compName+num;
									identicalComponents.at(is).push_back(newCompName);
									matchedSiteName = true;
								}
							}
							if(!matchedSiteName) {  // We have to handle the first matched site differently
								num="2";
								newCompName = compName+num;
								vector <string> v;
								v.push_back(compName+"1");
								v.push_back(newCompName);
								identicalComponents.push_back(v);
							}
							compName = newCompName;
							//cout<<"giving name: "<<compName<<endl;
							break;
						}
					}


					compLabels.push_back(compName);

					bool isIntegerState=false;

					//Take a look at the list of allowed states - this is where we determine if the component
					//site should be an integer state or a string
					vector <string> possibleStateNames;
					TiXmlElement *pListOfAllowedStates = pComp->FirstChildElement("ListOfAllowedStates");
					if(pListOfAllowedStates)
					{

						// First, loop through all the allowed states to determine if we need to parse this as an
						// integer site or a string site, by counting up the number of sting or integer sites
						int stringComponentNameCount = 0;
						int intComponentNameCount = 0;  int maxStateValue = -1;
						bool throwStateValLessThanZeroError = false;
						int plusMinusComponentNameCount = 0;

						TiXmlElement *pAlStates=0;
						for ( pAlStates = pListOfAllowedStates->FirstChildElement("AllowedState"); pAlStates != 0; pAlStates = pAlStates->NextSiblingElement("AllowedState"))
						{
							if(!pAlStates->Attribute("id")) {
								cerr<<"!!!Error:  AllowedState tag in ComponentType '"+compName+"' of MoleculeType: '" +
									typeName + "' must contain the id attribute.  Quitting."<<endl;
								return false;
							}

							string aState = pAlStates->Attribute("id");
							try {
								int stateVal= NFutil::convertToInt(aState);
								//Also check now to make sure that stateVal is nonnegative - we don't allow that because
								//of the way states are indexed
								if(stateVal<0) { throwStateValLessThanZeroError = true; }
								if(stateVal>maxStateValue) { maxStateValue = stateVal; }
								intComponentNameCount++;
							} catch (std::runtime_error e) {
								//If we catch an error, than the state value was a string (not an integer) and
								//we should handle as such.
								if(aState!="PLUS" && aState!="MINUS") {
									stringComponentNameCount++;
								}
								else {
									plusMinusComponentNameCount++;
								}
							}
						}


						// If we should parse the states as strings...
						if(stringComponentNameCount>0) {
							if(plusMinusComponentNameCount>0) {
								cerr<<"Warning!! when creating Component '"+compName+"' of MoleculeType: '"+
									typeName +"'.\n  State values are being parsed as a string, and you " +
									"have included 'PLUS' or 'MINUS' as state labels.  In this context, " +
									"'PLUS' and 'MINUS' do NOT increment or decrement the state values.\n" +
									"If you want to increment or decrement an integer state, there can be no " +
									"strings used as component names!  All state values must be integers.  So be careful!" <<endl;
							}

							int allowedStateCount = 0;
							for ( pAlStates = pListOfAllowedStates->FirstChildElement("AllowedState"); pAlStates != 0; pAlStates = pAlStates->NextSiblingElement("AllowedState"))
							{
								string aState = pAlStates->Attribute("id");
								if(allowedStates.find(typeName+"_"+compName+"_"+aState)!=allowedStates.end()) continue;
								allowedStates[typeName+"_"+compName+"_"+aState] = allowedStateCount;
								//cout<<"\nadding allowed state: "<<typeName+"_"+compName+"_"+aState<<" to be "<<allowedStateCount;
								allowedStateCount++;
								possibleStateNames.push_back(aState);
								if(verbose) cout<<"~"<<aState;
							}
							isIntegerState = false;
						}

						// If we should parse as integer states
						else if (intComponentNameCount>0 && stringComponentNameCount==0) {

							//First make sure that all the state values were greater than zero
							if(throwStateValLessThanZeroError || maxStateValue<0) {
								cerr<<"Error when creating Component '"+compName+"' of MoleculeType: '"+
									typeName +"'.\n  State values cannot be negative Integers!" +
									" Check your moleculeTypes."<<endl;
								return false;
							}
							//Second, check that the max state value was less than the limit
							if(maxStateValue>10000) {
								cerr<<"Error when creating Component '"+compName+"' of MoleculeType: '"+
									typeName +"'.\n  State values that are Integers cannot exceed" +
									" a value of 10,000!  Recheck your moleculeTypes."<<endl;
								return false;
							}

							//Finally, we have the max value of our integers, and we know that they are correct!  So
							//Let's build the list of possible state Integer values, and get the heck out of this
							//annoying loop.
							string stateString=""; int allowedStateCount = 0;
							for(int sv = 0; sv<=maxStateValue; sv++) {
								stringstream out; out << sv;
								string stateString = out.str();
								allowedStates[typeName+"_"+compName+"_"+stateString] = sv;
								allowedStateCount++;
								possibleStateNames.push_back(stateString);
							}
							if(verbose) cout<<"~integer[0-"<<maxStateValue<<"]";

							//Now we remember that this is an integer state, and break because we have already
							//created all possible state values for this component
							isIntegerState=true;
						}

						else { /* no possible states... */ }


					}

					//register if the component has an integer state, if so, remember that
					if(isIntegerState) isIntegerComponent.push_back(true);
					else isIntegerComponent.push_back(false);

					//register the possible component states and names, with the default set
					//at the first allowed state given, or "" if no allowed states are provided
					possibleComponentStates.push_back(possibleStateNames);
					if(possibleStateNames.size()>=1) defaultCompState.push_back(possibleStateNames.at(0));
					else defaultCompState.push_back("");
				}
			}
			if(verbose) cout<<")"<<endl;
			if(verbose) {
				for(unsigned int is=0; is<identicalComponents.size(); is++)  {
					cout<<"\t\t\t-Identified Equivalent Components: ";
					for(unsigned int isvec=0; isvec<identicalComponents.at(is).size(); isvec++) {
						cout<<identicalComponents.at(is).at(isvec)<<" ";
					} cout<<endl;
				}
			}



			//Go back and set the first symmetric component label to be 'compName1' so we know
			//immediately that they are symmetric sites (have to add in the possible binding site
			//names as well!
			for(unsigned int k=0; k<firstSymSiteToAppend.size(); k++) {
				string originalCompLabel = compLabels.at(firstSymSiteToAppend.at(k));
				compLabels.at(firstSymSiteToAppend.at(k)) = compLabels.at(firstSymSiteToAppend.at(k))+"1";

				string oldKeyStart = typeName+"_"+originalCompLabel+"_";
				//cout<<":::"<<typeName+"_"+originalCompLabel+"1_"<<endl;
				map<string,int>::iterator it;
				for ( it=allowedStates.begin() ; it != allowedStates.end(); it++ ) {
					string mappedKey = (*it).first;
					if(mappedKey.size()<=oldKeyStart.size()) continue;
					//cout<<mappedKey.substr(0,50)<<"  "<<mappedKey.substr(0,oldKeyStart.size())<<endl;

					if(mappedKey.substr(0,oldKeyStart.size())==oldKeyStart) {
						//cout<<"Found!! :"<<mappedKey.substr(0,oldKeyStart.size())<<endl;
						allowedStates[typeName+"_"+originalCompLabel+"1_"+mappedKey.substr(oldKeyStart.size())] = (*it).second;
						//cout<<mappedKey.substr(oldKeyStart.size())<<endl;
					}
				}
			}



			//Create the moleculeType and register any symmetric sites we may have...
			MoleculeType *mt = new MoleculeType(typeName,compLabels,defaultCompState,possibleComponentStates,isIntegerComponent,isPopulation,s);
			mt->addEquivalentComponents(identicalComponents);

			//Finally, clear the states and binding site labels that we read in
			compLabels.clear();
			defaultCompState.clear();
			possibleComponentStates.clear();
			identicalComponents.clear();
			isIntegerComponent.clear();
			firstSymSiteToAppend.clear();
		}

		// prints out allowed state map
		//for ( std::map< string, int, std::less< int > >::const_iterator iter = allowedStates.begin();
		//      iter != allowedStates.end(); ++iter )
		//      cout << iter->first << '\t' << iter->second << '\n';


		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		cerr<<"Caught some unknown error when creating MoleculeTypes."<<endl;
		return false;
	}

	return false;
}





bool NFinput::initStartSpecies(
		TiXmlElement * pListOfSpecies,
		System * s,
		map <string,double> &parameter,
		map<string,int> &allowedStates,
		bool verbose)
{
	////map<string,int>::iterator iter;
	////  for( iter = allowedStates.begin(); iter != allowedStates.end(); iter++ ) {
	////    cout << "state: " << iter->first << ", value: " << iter->second << endl;
	////  }


	try {
		//A vector to hold molecules as we are creating the species
		vector < vector <Molecule *> > molecules;

		//A vector that maps binding site ids into a molecule location in the molecules vector
		//and the name of the binding site
		map <string, string> bSiteSiteMapping;
		map <string, int> bSiteMolMapping;

		vector <string> stateName;
		vector <double> stateValue;

		vector<string>::iterator snIter;


		//Loop through all the species
		TiXmlElement *pSpec;
		for ( pSpec = pListOfSpecies->FirstChildElement("Species"); pSpec != 0; pSpec = pSpec->NextSiblingElement("Species"))
		{
			//First get the species name and make sure it exists
			string speciesName;
			if(!pSpec->Attribute("id")) {
				cerr<<"Species tag without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			} else {
				speciesName = pSpec->Attribute("id");
			}



			//Get the number of molecules of this species to create
			string specCount;
			if(!pSpec->Attribute("concentration")) {
				cerr<<"Species "<<speciesName<<" does not have a 'concentration' attribute.  Quitting"<<endl;
				return false;
			} else {
				specCount = pSpec->Attribute("concentration");
			}

			//Try to parse out the number of this species, or look it up in the parameter map
			int specCountInteger=0;
			try {
				specCountInteger = NFutil::convertToInt(specCount);
			} catch (std::runtime_error &e1) {

				// if we cannot get it as an integer, try as a double (for instance, for notation
				// such as 2e4).  We cast it as an int, which will always round the number down
				// to the nearest whole integer.
				try {
					specCountInteger = (int) NFutil::convertToDouble(specCount);

				} catch (std::runtime_error &e1) {
					// if we cannot get it as an integer, try as a double (for instance, for notation
					// such as 2e4).  We cast it as an int, which will always round the number down
					// to the nearest whole integer.
					try {
						specCountInteger = (int) NFutil::convertToDouble(specCount);
					} catch (std::runtime_error &e1) {
						if(parameter.find(specCount)==parameter.end()) {
							cerr<<"Could not find parameter: "<<specCount<<" when creating species "<<speciesName<<". Quitting"<<endl;
							return false;
						}
						specCountInteger = (int)parameter.find(specCount)->second;
					}

				}
			}

			//Make sure we didn't try to create a negative number of molecules
			if(specCountInteger<0) {
				cerr<<"I cannot, in good conscience, make a negative number ("<<specCount<<") of species when creating species "<<speciesName<<". Quitting"<<endl;
				return false;
			}

			// Removed this next check!!
			// We need to instantiate a population species, even if it does not
			// currently have a positive species count.  --Justin.
			/*
			// If we're not going to make anything, well, then, don't make anything silly!  Stop here!
			if(specCountInteger==0) {
				if(verbose) cout<<"\t\tNot creating any instances of the Species: "<<speciesName<<" because you said I should make zero of them."<<endl;
				continue;
			}
			*/


			//Make sure we have some molecules in our list of species
			TiXmlElement *pListOfMol = pSpec->FirstChildElement("ListOfMolecules");
			if(!pListOfMol) {
				cerr<<"Species "<<speciesName<<" contains no molecules!  I think that was a mistake, on your part, so I'm done."<<endl;
				return false;
			}


			// Give our users a nice little message...
			if(verbose) cout<<"\t\tCreating "<<specCountInteger<<" instances of the Species: "<<speciesName<<endl;

			/////////////////////////////////////////////////////////////////////////////////////////////////
			// Now loop through the molecules
			bool found_population = false;
			bool found_particle = false;
			TiXmlElement *pMol;
			for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule"))
			{
				// only one molecule per population species!!  --justin
				if ( found_population )
				{
					cerr << "!!!Error.  More than one population molecule when creating species '"
					     << speciesName << "'. Quitting"<<endl;
					return false;
				}

				//First get the type of molecule and retrieve the moleculeType object from the system
				string molName, molUid;
				if(!pMol->Attribute("name") || ! pMol->Attribute("id"))  {
					cerr<<"!!!Error.  Invalid 'Molecule' tag found when creating species '"<<speciesName<<"'. Quitting"<<endl;
					return false;
				} else {
					molName = pMol->Attribute("name");
					molUid = pMol->Attribute("id");
				}

				//Skip anything that is null or trash molecule
				if(molName=="Null" || molName=="NULL" || molName=="null") {
					if(verbose) cout<<"\t\t\tSkipping a null molecule in species declaration"<<endl;
					continue;
				}
				if(molName=="Trash" || molName=="TRASH" || molName=="trash" ) {
					if(verbose) cout<<"\t\t\tSkipping a trash molecule in species declaration"<<endl;
					continue;
				}


				// Identify the moleculeType if we can (note that this call could potentially kill our code if we can't find the type);
				MoleculeType *mt = s->getMoleculeTypeByName(molName);
				if(verbose) cout<<"\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;

				// check if this is a population molecule type  --justin
				if ( mt->isPopulationType() )
					found_population = true;
				else
					found_particle = true;

				// we can't mix populations and agents!  --justin
				if ( found_population  &&  found_particle )
				{
					cerr << "!!!Error.  Found mixed population and agent molecule types when creating species '"
					     << speciesName << "'. Quitting"<<endl;
					return false;
				}

				vector <string> usedComponentNames;

				//Loop through the components of the molecule in order to set state values
				TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
				if(pListOfComp)
				{
					TiXmlElement *pComp;
					for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component"))
					{
					
						//Get the basic properties of the component
						string compId,compName,compBondCount;
						if(!pComp->Attribute("id") || !pComp->Attribute("name") || !pComp->Attribute("numberOfBonds")) {
							cerr<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of species '"<<speciesName<<"'. Quitting"<<endl;
							return false;
						} else {
							compId=pComp->Attribute("id");
							compName = pComp->Attribute("name");
							compBondCount = pComp->Attribute("numberOfBonds");
						}

					/*	// Get component index from name (first symmetric component ok),
						// then check if this is an integer component.  --justin
						bool isIntegerComp;
						int compIndex;
						compIndex = mt->getCompIndexFromName(compName);
						isIntegerComp = mt->isIntegerComponent( compIndex );
					*/

						//cout<<" parsing comp: "<<compName<<endl;

						//First, if the site is symmetric, we have to relabel it correctly...
						if(mt->isEquivalentComponent(compName)) {
							//cout<<"is eq"<<endl;
							int *eqCompClass; int n_eqComp;
							mt->getEquivalencyClass(eqCompClass,n_eqComp,compName);


							bool couldPlaceSymComp=false;
							for(int eq=0; eq<n_eqComp;eq++) {
								string eqCompNameToCompare=mt->getComponentName(eqCompClass[eq]);
								//cout<<"comparing to: "<<eqCompNameToCompare<<endl;
								bool foundMatch=false;
								for(unsigned int ucn=0;ucn<usedComponentNames.size(); ucn++) {
									if(usedComponentNames.at(ucn).compare(eqCompNameToCompare)==0) {
										foundMatch=true; break;
									}
								}
								if(!foundMatch) {
									//cout<<" not used, using."<<endl;
									usedComponentNames.push_back(eqCompNameToCompare);
									compName=eqCompNameToCompare;
									couldPlaceSymComp=true;
									break;
								} else {
									//cout<<" used, moving on."<<endl;
								}
							}
							if(!couldPlaceSymComp) {
								cout<<"Too many symmetric sites specified, when creating species: "<<speciesName<<endl;
								return false;
							}
						} else {
							for(unsigned int ucn=0;ucn<usedComponentNames.size(); ucn++) {
								if(usedComponentNames.at(ucn).compare(compName)==0) {
									cout<<"Specified the same component multiple times, when creating species: "<<speciesName<<endl;
									return false;
								}
							}
							usedComponentNames.push_back(compName);
						}


						//If it is a state, treat it as such
						string compStateValue;
						if(pComp->Attribute("state"))
						{
							//First grab the states value as a string
							compStateValue = pComp->Attribute("state");
							if(allowedStates.find(molName+"_"+compName+"_"+compStateValue)==allowedStates.end()) {
								cerr<<"You are trying to create a molecule of type '"<<molName<<"', but you gave an "<<endl;
								cerr<<"invalid state! The state you gave was: '"<<compStateValue<<"'.  Quitting now."<<endl;
								return false;
							} else {

								//State is a valid allowed state, so push it onto our list
								int stateValueInt = allowedStates.find(molName+"_"+compName+"_"+compStateValue)->second;
								stateName.push_back(compName);
								stateValue.push_back(stateValueInt);
							}
						}


						//finally, we have to add the b site mapping that will let us later
						//easily connect binding sites with the molecules involved
						bSiteSiteMapping[compId] = compName;
						bSiteMolMapping[compId] = molecules.size();
					}
				}
				else
				{
					if ( mt->getNumOfComponents() > 0 )	{
						cout<<"!!! warning: no list of components specified for molecule: '"<<molUid<<"' of species '"<<speciesName<<"'"<<endl;
					}
				}
				usedComponentNames.clear();


				//We dont' have to do this anymore, because we handled it earlier!
				//int eqClassCount = mt->getNumOfEquivalencyClasses();
				//int *currentCount = new int[eqClassCount];
				//for(int i=0; i<eqClassCount; i++) { currentCount[i]=1; }
				//string *eqCompNames = mt->getEquivalencyClassCompNames();

				//loop to create the actual molecules of this type
				vector <Molecule *> currentM;
				molecules.push_back(currentM);

				if ( !found_population )
				{
					for(int m=0; m<specCountInteger; m++)
					{
						Molecule *mol = mt->genDefaultMolecule();

						//for(int i=0; i<eqClassCount; i++) { currentCount[i]=1; }

						//Loop through the states and set the ones we need to set
						int k=0;
						for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
						{
							//if(mt->isEquivalentComponent((*snIter))) {
							//	int eqNum = mt->getEquivalencyClassNumber((*snIter));
							//	std::stringstream numStream; numStream << currentCount[eqNum];
							//	string postFix = numStream.str();
							//	m->setComponentState((*snIter)+postFix, (int)stateValue.at(k));
							//	currentCount[eqNum]++;
							//} else {
								mol->setComponentState((*snIter), (int)stateValue.at(k));
							//}
						}

						molecules.at(molecules.size()-1).push_back(mol);
					}
				}
				// handle population case (only create one instance of this molecule type) --Justin
				else
				{
					Molecule *mol = mt->genDefaultMolecule();
					// set population
					mol->setPopulation( specCountInteger );

					//Loop through the states and set the ones we need to set
					int k=0;
					for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
						mol->setComponentState((*snIter), (int)stateValue.at(k));

					molecules.at(molecules.size()-1).push_back(mol);

				}

				//delete [] currentCount;

				//Reset the states for the next wave...
				stateName.clear();
				stateValue.clear();
			}


			///////////////////////////////////////////////////////////////
			//Here is where we add the bonds to the molecules in this species
			TiXmlElement *pListOfBonds = pListOfMol->NextSiblingElement("ListOfBonds");
			if(pListOfBonds)
			{
				// Complain and quit if this is a population species (no bonds allowed!)
				//  --Justin, 4Mar2011
				if ( found_population )
				{
					cerr << "!! Attempt to create illegal bond in population species: "
					     << speciesName << ".  Quitting." << endl;
					return false;
				}

				//First get the information on the bonds in the complex
				TiXmlElement *pBond;
				for ( pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond"))
				{
					string bondId, bSite1, bSite2;
					if(!pBond->Attribute("id") || !pBond->Attribute("site1") || !pBond->Attribute("site2")) {
						cerr<<"!! Invalid Bond tag for species: "<<speciesName<<".  Quitting."<<endl;
						return false;
					} else {
						bondId = pBond->Attribute("id");
						bSite1 = pBond->Attribute("site1");
						bSite2 = pBond->Attribute("site2");
					}
					//cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;


					//Get the information on this bond that tells us which molecules to connect
					try {
						string bSiteName1 = bSiteSiteMapping.find(bSite1)->second;
						int bSiteMolIndex1 = bSiteMolMapping.find(bSite1)->second;
						string bSiteName2 = bSiteSiteMapping.find(bSite2)->second;
						int bSiteMolIndex2 = bSiteMolMapping.find(bSite2)->second;

						for(int j=0;j<specCountInteger;j++) {
							Molecule::bind( molecules.at(bSiteMolIndex1).at(j),bSiteName1.c_str(),
											molecules.at(bSiteMolIndex2).at(j),bSiteName2.c_str());
						}

					} catch (exception& e) {
						cout<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating species '"<<speciesName<<"'. Quitting"<<endl;
						return false;
					}
				}
			}

			//Tidy up and clear the lists for the next species
			vector< vector <Molecule *> >::iterator mIter;
			for(mIter = molecules.begin(); mIter != molecules.end(); mIter++ ) {
				(*mIter).clear();
			}

			molecules.clear();
			bSiteMolMapping.clear();
			bSiteSiteMapping.clear();
		}


		//s->printAllMoleculeTypes();


		//If we got here, then we are indeed successful
		return true;
	} catch (...) {
		cerr<<"Caught some unknown error when creating Species."<<endl;
		return false;
	}

	return false;
}







bool NFinput::initReactionRules(
		TiXmlElement * pListOfReactionRules,
		System * s,
		map <string,double> &parameter,
		map<string,int> &allowedStates,
		bool blockSameComplexBinding,
		bool verbose,
		int &suggestedTraversalLimit)
{


	try {

		//First, loop through all the rules
		TiXmlElement *pRxnRule;
		for ( pRxnRule = pListOfReactionRules->FirstChildElement("ReactionRule"); pRxnRule != 0; pRxnRule = pRxnRule->NextSiblingElement("ReactionRule"))
		{

			//First, scan the reaction rule for possible symmetries!!!
			map <string, component> symComps;
			map <string, component> symRxnCenter;

			if(!FindReactionRuleSymmetry(pRxnRule, s,
									parameter,
									allowedStates,
									symComps,
									symRxnCenter,
									verbose)) return false;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Begin with some basic parsing of the rules and reactant patterns
			//cout<<symComps.size()<<"  ----  "<<symRxnCenter.size()<<endl;

			//For each possible permuation of the reaction rule, let us create a separate reaction
			//to keep track of the result...
			vector < map <string,component> > permutations;
			generateRxnPermutations(permutations, symComps, symRxnCenter,verbose);


			for( unsigned int p=0; p<permutations.size(); p++)
			{
				map <string,component> symMap = permutations.at(p);


				//Grab the name of the rule
				string rxnName;
				if(!pRxnRule->Attribute("id")) {
					cerr<<"ReactionRule tag without a valid 'id' attribute.  Quiting"<<endl;
					return false;
				} else {
					rxnName = pRxnRule->Attribute("id");
					if(pRxnRule->Attribute("name")) {
						rxnName = pRxnRule->Attribute("name");
					}

					if(permutations.size()>1) {
						stringstream out; out << (p+1);
						rxnName = rxnName + "_sym" + out.str();
					}
				}
				if(verbose) cout<<"\t\tCreating Reaction Rule: "<<rxnName<<endl;

				// grab symmetry factor, if any..
				bool useSymmetryFactor = false;
				double symmetryFactor = 1.0;
				if (pRxnRule->Attribute("symmetry_factor"))
				{
					try {
						symmetryFactor = NFutil::convertToDouble(pRxnRule->Attribute("symmetry_factor"));
						if ( symmetryFactor <= 0.0 ) throw std::runtime_error("Symmetry Factor must be a positive value!");
						useSymmetryFactor = true;
						if (verbose) cout << "\t\t\tUsing symmetry factor = " << symmetryFactor << endl;
					} catch (std::runtime_error &e1) {
						cerr<<"Error!! Symmetry Factor for ReactionRule "<<rxnName<<" was not set properly.  quitting."<<endl;
						exit(1);
					}
				}

				//First, read in the template molecules using these data structures
				// maps reactant pattern ids to TemplateMolecule pointers
				map <string,TemplateMolecule *> reactants;
				// maps component ids to component objects
				map <string, component> comps;
				// points to TemplateMolecules for Reactants
				vector <TemplateMolecule *> templates;
				// points to TemplateMolecules for AddMoleculeTransforms
				vector <TemplateMolecule *> addmol_templates;



				//  Read in the Reactant Patterns for this rule
				TiXmlElement *pListOfReactantPatterns = pRxnRule->FirstChildElement("ListOfReactantPatterns");
				if(!pListOfReactantPatterns) {
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no reactant patterns!"<<endl;
					continue;
				}

				TiXmlElement *pReactant;
				for ( pReactant = pListOfReactantPatterns->FirstChildElement("ReactantPattern"); pReactant != 0; pReactant = pReactant->NextSiblingElement("ReactantPattern"))
				{
					const char *reactantName = pReactant->Attribute("id");
					if(!reactantName) {
						cerr<<"Reactant tag in reaction "<<rxnName<<" without a valid 'id' attribute.  Quiting"<<endl;
						return false;
					}
					if(verbose) cout<<"\t\t\tReading Reactant Pattern: "<<reactantName<<endl;

					TiXmlElement *pListOfMols = pReactant->FirstChildElement("ListOfMolecules");
					if(pListOfMols) {
						TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, allowedStates, reactantName, reactants, comps, symMap, verbose, suggestedTraversalLimit);
						if(tm==NULL) return false;
						templates.push_back(tm);
					}
					else {
						cerr<<"Reactant pattern "<<reactantName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
						return false;
					}
				}

				//Outputting all the templates for debugging purposes
				//map<const char*, TemplateMolecule *, strCmp>::iterator it;
				//	for ( it=reactants.begin() ; it != reactants.end(); it++ )
				//		cout << (*it).first << " => " << (*it).second->getMoleculeType()->getName() << endl;


				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Read in the list of operations we need to perform in this rule
				TiXmlElement *pListOfOperations = pRxnRule->FirstChildElement("ListOfOperations");
				if ( !pListOfOperations )
				{
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no operations!  This rule will do nothing!"<<endl;
					continue;
				}


				//Find new molecule creations before we creating the transformation set!
				// This has to be done first so we know how many extra mappingSets we'll
				// need for new product molecules in the TransformationSet
				//   --Justin, 2May2011

				//As new product molecules are identified, we'll create components and
				// add them to the 'comps' map. THen we'll be able to apply further
				// transformations to the new molecules later on.

				// I'm pretty sure we don't need to worry about symmetry here.
				//   --Justin

				// remember added product ids, to avoid adding something more than once
				vector <string>            addedProductPatterns;
				// remember MoleculeCreator objects
				vector <MoleculeCreator *> moleculeCreatorsList;

				TiXmlElement *pAdd;
				for ( pAdd = pListOfOperations->FirstChildElement("Add");
						pAdd != 0;  pAdd = pAdd->NextSiblingElement("Add") )
				{

					//Make sure all the information about the state change is here
					string id;       // the full molecule id
					string patt_id;  // the pattern id

					if ( !pAdd->Attribute("id") )
					{
						cerr << "A specified Add operation in ReactionClass: '" << rxnName << "' does not\n"
						     << "have a valid id attribute.  Quitting." << endl;
						return false;
					}

					id = pAdd->Attribute("id");


					//Now make sure we're not adding something twice...
					bool foundAddedPattern = false;
					for ( unsigned int appIndex=0; appIndex<addedProductPatterns.size(); appIndex++ )
					{
						if( addedProductPatterns.at(appIndex).compare(id)==0 )
						{
							foundAddedPattern = true;
							break;
						}
					}
					if ( !foundAddedPattern )
					{
						addedProductPatterns.push_back(id);
					}
					else continue;


					if(verbose) cout<<"\t\t\t***Identified addition of product: "+id+"."<<endl;

					// Some preprocessing of the ID to make sure we are getting the right thing
					//First, check if we are pointing to a molecule - if we are, we have to stop
					//this nonsense! (current BNGL returns add transforms for each molecule, not
					//per the entire pattern...)  [Michael]

					// Pointing to molecules is more general, because we might
					//  want to create a molecule and then bind it to an existing molecule.
					//  And if we need to create an entire species, we create each molecule and
					//  then perform bond operations. Going forward, we should create new molecules
					//  one at a time. --Justin, 2Mar2011

					// Extract pattern substring
					int M_position = id.find("M");
					if ( M_position >= 0 )
					{
						patt_id = id.substr(0,M_position-1);
					}
					else
					{
						cerr << "Error:: ReactionClass " << rxnName << " contains an Add transform\n"
						     << "with id that does not refer to a molecule." << endl;
						return false;
					}


					//Go get the product pattern we need which will specify how to make this new species
					TiXmlElement *pListOfProductPatterns = pRxnRule->FirstChildElement("ListOfProductPatterns");
					if ( !pListOfProductPatterns )
					{
						cerr << "Error:: ReactionRule " << rxnName << " contains no product patterns,\n"
						     << "but needs at least one to add a molecule creation rule!" << endl;
						return false;
					}

					bool found_mol = false;
					TiXmlElement * pProduct;
					for ( pProduct = pListOfProductPatterns->FirstChildElement("ProductPattern");
							pProduct != 0; pProduct = pProduct->NextSiblingElement("ProductPattern") )
					{
						//First extract out the product Id
						string productName;
						if( !pProduct->Attribute("id") )
						{
							cerr<<"Product pattern in ReactionRule "+rxnName+" does not have a valid id attribute!"<<endl;
							return false;
						}
						productName = pProduct->Attribute("id");

						// If this isn't the pattern we're looking for, go onto the next pattern
						if ( productName != patt_id ) continue;

						TiXmlElement *pListOfMols = pProduct->FirstChildElement("ListOfMolecules");

						if( pListOfMols==NULL )
						{
							cerr << "Reactant product pattern " << productName << " in reaction\n"
							     << rxnName << " without a valid 'ListOfMolecules'!  Quiting." << endl;
							return false;
						}

						TiXmlElement * pMolecule;
						for ( pMolecule = pListOfMols->FirstChildElement("Molecule");
								pMolecule != 0;  pMolecule = pMolecule->NextSiblingElement("Molecule") )
						{
							//First extract out the molecule Id
							string moleculeName;
							if ( !pMolecule->Attribute("id") )
							{
								cerr << "Product molecule in ReactionRule " << rxnName << " does not\n"
								     << "have a valid id attribute!" << endl;
								return false;
							}

							moleculeName = pMolecule->Attribute("id");

							// Is this the molecule we're looking for?
							if ( moleculeName != id ) continue;


							// This will create templates and components for the product molecule and its sites
							//  (also handles increment population transforms!)
							bool ok = NFinput::readProductMolecule(
									      pMolecule, s, parameter, allowedStates,
										  productName, moleculeCreatorsList,
										  comps, verbose );

							if( !ok )
							{
								cerr << "Some problem reading product molecule for Add transform in reaction "
								     << rxnName << "." << endl;
								return false;
							}

							// if we haven't found a problem yet, then we found the right molecule
							//  (even if the molecule is being ignored, like null
							found_mol = true;
							break;
						}
						if (found_mol) break;
					}

					// make sure we found the molecule!
					if ( !found_mol )
					{
						cerr << "Could not find product molecule '" << id << "' referenced in\n"
						     << "Add transform in reaction " << rxnName << "." << endl;
						return false;
					}
				}



				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Create the TransformationSet so that we can collect all the operations that are specified for this rule
				TransformationSet *ts;
				if ( moleculeCreatorsList.empty() )
				{
					// create transformation set without add molecule transformas
					ts = new TransformationSet(templates);
				}
				else
				{
					// get the templates for added molecules
					vector <MoleculeCreator *>::iterator  mc_iter;
					for ( mc_iter = moleculeCreatorsList.begin();
							mc_iter != moleculeCreatorsList.end();  ++mc_iter )
					{
						addmol_templates.push_back( (*mc_iter)->getTemplateMolecule() );
					}

					// create transformation set with add molecule transforms!
					ts = new TransformationSet( templates, addmol_templates );

					// add molecule creators!
					for ( mc_iter = moleculeCreatorsList.begin();
							mc_iter != moleculeCreatorsList.end();  ++mc_iter )
					{
						if( !ts->addAddMolecule( *mc_iter ) ) return false;
					}
				}

				// Should we use complex bookkeeping?
				ts->setComplexBookkeeping( blockSameComplexBinding );

				// Are we using a symmetry factor?
				if (useSymmetryFactor)
				{
					ts->setSymmetryFactor(symmetryFactor);
				}


				//Next extract out the state changes
				TiXmlElement *pStateChange;
				for ( pStateChange = pListOfOperations->FirstChildElement("StateChange"); pStateChange != 0; pStateChange = pStateChange->NextSiblingElement("StateChange"))
				{
					//Make sure all the information about the state change is here
					string site, finalState;
					if(!pStateChange->Attribute("site") || !pStateChange->Attribute("finalState")) {
						cerr<<"A specified state change operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid site or finalState attribute.  Quitting."<<endl;
						return false;
					} else {
						site = pStateChange->Attribute("site");
						finalState = pStateChange->Attribute("finalState");
					}

					//First grab the component that is going to change...
					component *c;
					int finalStateInt = 0;
					if(!lookup(c, site, comps, symMap)) return false;

					// Check if this is modifying a population (illegal!)
					if ( c->t->getMoleculeType()->isPopulationType() )
					{
						cerr << "Attempt to change state of a population type molecule in ReactionClass: '"
						     << rxnName << "'. Quitting." << endl;
						return false;
					}

					//handle both increment and decrement states first...
					if(finalState=="PLUS") {
						if(!ts->addIncrementStateTransform(c->t,c->symPermutationName)) return false;

						if(verbose) {
							cout<<"\t\t\t***Identified increment state of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
							cout<<") to new state value: oldStateValue+1"<<endl;
						}
					}
					else if(finalState=="MINUS") {
						if(!ts->addDecrementStateTransform(c->t,c->symPermutationName)) return false;

						if(verbose) {
							cout<<"\t\t\t***Identified decrement state of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
							cout<<") to new state value: oldStateValue-1"<<endl;
						}

					}
					else {

						//Here, we handle your typical state change operation
						try {

							string lookupname = c->symPermutationName;
							if(allowedStates.find(c->t->getMoleculeTypeName()+"_"+lookupname+"_"+finalState)==allowedStates.end()) {

								//if(c->t->getMoleculeType()->isEquivalentComponent(c->name)) {
								//	lookupname = c->name;
								//}
								//if(allowedStates.find(c->t->getMoleculeTypeName()+"_"+lookupname+"_"+finalState)==allowedStates.end()) {
									cout<<"Error! in NFinput, when looking up state: "<<c->t->getMoleculeTypeName()+"_"+c->symPermutationName+"_"+finalState<<endl;
									cout<<"Could not find this in the list of allowed states!  exiting!"<<endl;
									exit(1);
								//}
							}
							finalStateInt = allowedStates.find(c->t->getMoleculeTypeName()+"_"+lookupname+"_"+finalState)->second;
							//cout<<"found:"<<finalStateInt<<endl;
						} catch (exception& e) {
							cerr<<"Error in adding a state change operation in ReactionClass: '"+rxnName+"'."<<endl;
							cerr<<"It seems that the final state is not valid."<<endl;
							return false;
						}
						if(!ts->addStateChangeTransform(c->t,c->symPermutationName,finalStateInt)) return false;
						if(verbose) {
							cout<<"\t\t\t***Identified state change of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
							cout<<") to new state value: " + finalState<<endl;
						}
					}
				}


				//Extract out removal of bonds (Must do this before creation of bonds!  Otherwise
				//it is possible to add a bond before another was deleted in the same rule!  This
				//would give you an error!)
				TiXmlElement *pDeleteBond;
				for ( pDeleteBond = pListOfOperations->FirstChildElement("DeleteBond");
						pDeleteBond != 0;  pDeleteBond = pDeleteBond->NextSiblingElement("DeleteBond") )
				{
					//Make sure all the information about the unbinding operation change is here
					string site1,site2;
					if(!pDeleteBond->Attribute("site1") || !pDeleteBond->Attribute("site2")) {
						cerr<<"A specified binding operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid site1 or site2 attribute.  Quitting."<<endl;
						return false;
					} else {
						site1 = pDeleteBond->Attribute("site1");
						site2 = pDeleteBond->Attribute("site2");
						//Skip this if we are messing with a bond in the product pattern....
						// @TODO:  FIX THIS!  should reject adds in molecule species that are newly added!
						//if(site1.find("RP")>=0 || site2.find("RP")>=0) continue;
					}



					component *c1;
					component *c2;
					if(!lookup(c1, site1, comps, symMap)) return false;
					if(!lookup(c2, site2, comps, symMap)) return false;

					//Even though we had to make sure both ends exist, and we really only need one transformation
					//we give both templates so we can check for symmetric unbinding
					if(!ts->addUnbindingTransform(c1->t, c1->symPermutationName, c2->t, c2->symPermutationName)) return false;
					if(verbose) {
						cout<<"\t\t\t***Identified unbinding of site: "+c1->t->getMoleculeTypeName()+"("+c1->symPermutationName + ")";
						cout<<" to site " + c2->t->getMoleculeTypeName()+"("+c2->symPermutationName<<")"<<endl;
					}

				}

				//Next extract out the new bonds that are formed
				TiXmlElement *pAddBond;
				for ( pAddBond = pListOfOperations->FirstChildElement("AddBond");
						pAddBond != 0;  pAddBond = pAddBond->NextSiblingElement("AddBond") )
				{
					//cout<<"adding binding transform!"<<endl;
					//Make sure all the information about the binding operation is here
					string site1, site2;
					if( !pAddBond->Attribute("site1") || !pAddBond->Attribute("site2") )
					{
						cerr << "A specified binding operation in ReactionClass: '" << rxnName << "' does not\n"
						     << "have a valid site1 or site2 attribute.  Quitting." << endl;
						return false;
					}
					else
					{
						site1 = pAddBond->Attribute("site1");
						site2 = pAddBond->Attribute("site2");
					}

					// NFsim now is able to add bonds to newly created molecules. While processing
					//  AddMoleculeTransforms, we created components corresponding to new molecules
					//  and inserted maps from site id to component object into 'comps'
					// --Justin, 3Mar2011

					// Adding newly created molecules that are bound turns out to not be compatible with
					// generalized checks for null conditions because in TransformationSet::transform, the
					// null conditions are checked before new molecules are created.  You cannot create
					// new molecules before you check for null conditions, because otherwise if a null
					// condition is met, then the newly created molecules would have to be removed from
					// the system.  But if the new molecules are not created yet, then the check for null
					// conditions creates a segfault in some cases because it is looking into a mappingset
					// that has a null mapping.  Checking for null mappings is insuffecient, because mappings
					// are not required to be cleared after comparisons or usage, thus there can be dangling
					// mappings which will not match null, but still may be incorrect.  The work around implemented
					// here addresses this by creating a new type of transformation that is unique for binding
					// reactions taking place between newly created molecules.  For these binding reactions, a null
					// condition can never be met because molecularity constraints cannot be broken by a new molecule,
					// because that molecule did not exist before.  To check for bonds to new molecules, we use the
					// component lookup name, which if pointing to a reactant (designated in the XML by ..._RP_...),
					// then we know it existed, but if it is pointing to a product (designated by ..._PP_...), we
					// know that it was created by this rule.
					// --michael Aug2011

					component *c1;
					component *c2;

					if ( !lookup( c1, site1, comps, symMap) ) return false;
					if ( !lookup( c2, site2, comps, symMap) ) return false;


					// check if these are bonds to new product molecules
					string::size_type PPpos1 = site1.find("_PP"); string::size_type PPpos2 = site2.find("_PP");
					// if either is a new molecule, set the flag so that we use the special
					// newMoleculeBinding transform instead of the standard transform.
					bool isNewMoleculeBond = false;
					if(PPpos1!=string::npos || PPpos2!=string::npos) {
						isNewMoleculeBond = true;
					}


					// Check if this is binding a population (illegal!)
					if (    c1->t->getMoleculeType()->isPopulationType()
						 || c2->t->getMoleculeType()->isPopulationType() )
					{
						cerr << "Attempt to AddBond at site on a population type molecule in ReactionClass: '"
						     << rxnName << "'. Quitting." << endl;
						return false;
					}

					//Make sure we only block binding on the same complex if they were on separate reactants
					//if this is an internal binding, then we have to allow it, even if we have the flag
					//that blocks same complex binding...
					//cout<<"\n"<<site1<<endl;
					//cout<<site2<<endl;

					// Comment about binding to new product molecules
					// 1) Binding two new molecules:
					// 2) Binding a new molecule and an existing molecule:

					int underScore1 = site1.find_first_of("_");
					int underScore2 = site2.find_first_of("_");

					string reactantNum1 = site1.substr( 0, site1.find_first_of("_",underScore1+1) );
					string reactantNum2 = site2.substr( 0, site2.find_first_of("_",underScore2+1) );
					//cout << "reactant1: " << reactantNum1 << endl;
					//cout << "reactant2: " << reactantNum2 << endl;

					// Handling inter- and intra-complex binding is now part of a general procedure
					//  for handling reaction molecularity  --Justin, 4Mar2011
					/*if ( reactantNum1.compare(reactantNum2)==0 )
					{
						//this means that they were on the same reactant, so we should always add
						//this as a normal binding reaction...
						if ( !ts->addBindingTransform( c1->t, c1->symPermutationName, c2->t, c2->symPermutationName) )
							return false;
					}
					else
					{
						//Otherwise, we should check how we should add this reaction, depending on the input flags
						if ( !blockSameComplexBinding )
						{
							if ( !ts->addBindingTransform( c1->t, c1->symPermutationName, c2->t, c2->symPermutationName) )
								return false;
						}
						else
						{
							if ( !ts->addBindingSeparateComplexTransform( c1->t, c1->symPermutationName, c2->t, c2->symPermutationName) )
								return false;
						}
						if (verbose)
						{
							cout << "\t\t\t***Identified binding of site: " << c1->t->getMoleculeTypeName()
							     << "(" << c1->symPermutationName << ")" << " to site " << c2->t->getMoleculeTypeName()
							     << "(" << c2->symPermutationName << ")" << endl;
						}
					}*/

					// Add binding transform
					if(isNewMoleculeBond) {
						if ( !ts->addNewMoleculeBindingTransform( c1->t, c1->symPermutationName, c2->t, c2->symPermutationName) )
							return false;
					} else {
						if ( !ts->addBindingTransform( c1->t, c1->symPermutationName, c2->t, c2->symPermutationName) )
							return false;
					}

					if (verbose)
					{
						cout << "\t\t\t***Identified binding of site: " << c1->t->getMoleculeTypeName()
						     << "(" << c1->symPermutationName << ")" << " to site " << c2->t->getMoleculeTypeName()
						     << "(" << c2->symPermutationName << ")" << endl;
						if(isNewMoleculeBond)
							cout << "\t\t\t   (this bond connects a new molecule that was synthesized by this rule)" << endl;
					}
				}


				//Next extract out anything that is destroyed
				TiXmlElement *pDelete;
				for ( pDelete = pListOfOperations->FirstChildElement("Delete");
						pDelete != 0;  pDelete = pDelete->NextSiblingElement("Delete") )
				{
					//Make sure all the information about the state change is here
					string id,delMolKeyword;
					if(!pDelete->Attribute("id") || !pDelete->Attribute("DeleteMolecules")) {
						cerr<<"A specified delete operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid id attribute or a DeleteMolecules attribute.  Quitting."<<endl;
						return false;
					} else {
						try {
							id = pDelete->Attribute("id");
							delMolKeyword = pDelete->Attribute("DeleteMolecules");


							//Here we can check if we are referencing a Molecule for deletion (As opposed to
							//the entire complex)
							int M_position = id.find_first_of("M");
							if(M_position>=0) {
								if(verbose) {
									cout<<"\t\t\t***Identified deletion of a molecule: "+id<<". ";
									if(delMolKeyword.compare("0")==0) {
										cout<<"\n\t\t\t\tDeleteMolecules keyword is turned off."<<endl;
									} else {
										cout<<"\n\t\t\t\tDeleteMolecules keyword is turned on, so only this specific molecule will be removed."<<endl;
									}
								}


								//Pointing to just a single molecule, so retrieve that molecule
								//cout<<"pointing to just a single molecule, so retrieve it."<<endl;
								component c = comps.find(id)->second;

								if(delMolKeyword.compare("0")==0) {
									//Pointing to a single molecule, no DeleteMolecules keyword, so delete it conditional
									cout<<"No delete molecule keyword, so delete this molecule, conditional on"
											"the fact that it cannot create multiple additional species."<<endl;

									//Give a warning, because the behavior of this function is very strange indeed!!!
									cerr<<"\nERROR!  You created a reaction ("+rxnName+") that deletes a molecule, but you did not use"<<endl;
									cerr<<"the 'DeleteMolecules' keyword.  Thus, conforming with BNGL specification, this rule will not"<<endl;
									cerr<<"delete the molecule IF the deletion creates two disjoint species that are no longer connected."<<endl;
									cerr<<"NFsim doesn't believe in this weird behavior, so NFsim enforces the use of the DeleteMolecules "<<endl;
									cerr<<"keyword in such situations, as in this rule.  So go update your rule now.\n"<<endl;
									return false;

									if ( c.mt->isPopulationType() )
									{
										if(!ts->addDeleteMolecule(c.t,TransformationFactory::DECREMENT_POPULATION)) return false;
									}
									else
									{
										if(!ts->addDeleteMolecule(c.t,TransformationFactory::DELETE_MOLECULES_NO_KEYWORD)) return false;
									}
								} else {
									//Pointing to a single molecule, DeleteMolecules keyword is on, so delete it regardless
									//cout<<"Using the DeleteMolecules keyword, so delete only the molecules"
									//	" that are being pointed to..."<<endl;
									if ( c.t->getMoleculeType()->isPopulationType() )
									{
										if(!ts->addDecrementPopulation(c.t)) return false;
									}
									else
									{
										if(!ts->addDeleteMolecule(c.t,TransformationFactory::DELETE_MOLECULES)) return false;
									}
								}
							}

							else if(comps.find(id)!=comps.end()) {
								if(verbose) {
									cout<<"\t\t\t***Identified deletion of the complete pattern: "+id<<". ";
									if(delMolKeyword.compare("0")==0) {
										cout<<"\n\t\t\t\tDeleteMolecules keyword is turned off, so the full species will be removed."<<endl;
									} else {
										cout<<"\n\t\t\t\tDeleteMolecules keyword is turned on, so only the molecules in the pattern will be removed."<<endl;
									}
								}

								// get component
								component c = comps.find(id)->second;
								if(delMolKeyword.compare("0")==0) {
									//Pointing to the whole complex, no DeleteMolecules keyword, so delete it all!
									//cout<<"Pointing to the whole complex, without DeleteMolecules keyword, so delete it and all connected..."<<endl;
									component c = comps.find(id)->second;
									if ( c.t->getMoleculeType()->isPopulationType() )
									{   // we're dealing with a population here, create a decrement population transform
										if(!ts->addDecrementPopulation(c.t)) return false;
									}
									else
									{	// we're dealing with a particle here, create regular delete transform
										if(!ts->addDeleteMolecule(c.t,TransformationFactory::COMPLETE_SPECIES_REMOVAL)) return false;
									}
								}
								else
								{
									//Pointing to the whole complex, DeleteMolecules keyword is turned on.  According to
									//BNGL logic, this will only delete the molecule(s) specified in the pattern, and not
									//the entire species.  Fortunately, BNG should detect that the keyword is on, and pass
									//deletion of individual molecules, not of the entire species.  so if we get here, then
									//there was an error in BNG!!

									cerr<<"\nERROR!  You created a reaction ("+rxnName+") that deletes a molecule, and you did use"<<endl;
									cerr<<"the 'DeleteMolecules' keyword.  Unfortunately, BNG provided NFsim with information about"<<endl;
									cerr<<"deleting an entire species.  This error should never happen, but if it does, contact us"<<endl;
									cerr<<"through the NFsim or BioNetGen websites to report this error.  "<<endl;

									return false;
								}




							} else {
								cerr<<"Error in adding an delete molecule operation in ReactionClass: '"+rxnName+"'."<<endl;

								cerr<<"Invalid read!  Looking for "<<id<<" but could not find it!"<<endl;
								exit(1);
							}
							//cout<<"Templates.size() "<<templates.size()<<endl;
							//c.t->printDetails();

						} catch (exception& e) {
							cerr<<"Error in adding an delete molecule operation in ReactionClass: '"+rxnName+"'."<<endl;
							cerr<<"It seems that I couldn't find the molecule to delete that you are refering to. (I was looking for ID: "<<pDelete->Attribute("id")<<endl;
							return false;
						}
					}

				}


				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// With the transforations now set, Let's actually create the reaction (remember to finalize the TransformationSet!
				//We can't finalize the transformation set here anymore! we have to do it just before we create
				//the reaction because we still have to add function pointers!!
				//ts->finalize();
				ReactionClass *r = 0;
				bool totalRateFlag=false;

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//  Read in the rate law for this reaction
				TiXmlElement *pRateLaw = pRxnRule->FirstChildElement("RateLaw");
				if(!pRateLaw){
					cerr<<"!!Error:: ReactionRule "<<rxnName<<" contains no rate law specification!"<<endl;
					return false;
				}

				if( !pRateLaw->Attribute("totalrate") ) {
					cerr<<"\n!!Error! This XML file was generated using an older version of BioNetGen that does not support the 'TotalRate' convention!"<<endl;
					cerr<<"You should upgrade your BioNetGen distribution now, or download the latest NFsim package, and regenerate this XML file."<<endl;
				} else {
					try {
						int rf = NFutil::convertToInt(pRateLaw->Attribute("totalrate"));
						if(rf>0) totalRateFlag=true;
						if(verbose) cout<<"\t\t\tTotal rate flag = "<<totalRateFlag<<endl;
					} catch (std::runtime_error &e1) {
						//cerr<<e1.what()<<endl;
						cerr<<"Error!! totalrate flag for ReactionRule "<<rxnName<<" was not set properly.  quitting."<<endl;
						exit(1);
					}
				}

				if(!pRateLaw->Attribute("id") || !pRateLaw->Attribute("type")) {
					cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification: cannot read 'id' or 'type' attribute!"<<endl;
						return false;
				}
				else {
					string rateLawName = pRateLaw->Attribute("id");
					string rateLawType = pRateLaw->Attribute("type");

					if(verbose) cout<<"\t\t\tRate Law for Reaction is: "<<rateLawType<<endl;
					if(rateLawType=="Ele")
					{
						//Create the Elementary Reaction...
						ts->finalize();

						// decide if we need a basic or population reaction
						bool found_pop = false;
						for ( unsigned int ir = 0;  ir < ts->getNreactants();  ++ir )
						{
							if ( ts->getTemplateMolecule(ir)->getMoleculeType()->isPopulationType() )
							{
								found_pop = true;
								break;
							}
						}

						r = new BasicRxnClass(rxnName,0,"",ts,s);



						//Make sure that the rate constant exists
						TiXmlElement *pListOfRateConstants = pRateLaw->FirstChildElement("ListOfRateConstants");
						if(!pListOfRateConstants) {
							cerr<<"Elementry Rate Law definition for "<<rxnName<<" does not have listOfRateConstants specified!  Quiting"<<endl;
							return false;
						}
						TiXmlElement *pRateConstant = pListOfRateConstants->FirstChildElement("RateConstant");
						if(!pRateConstant) {
							cerr<<"Elementry Rate Law definition for "<<rxnName<<" does not have RateConstants specified!  Quiting"<<endl;
							return false;
						} else {

							//Get the rate constant value
							string rateValue;
							if(!pRateConstant->Attribute("value")) {
								cerr<<"Elementry Rate Law definition for "<<rxnName<<" does not have a valid RateConstant value!  Quiting"<<endl;
								return false;
							} else {
								rateValue = pRateConstant->Attribute("value");
							}

							//Try to parse it into a double value or look it up in the parameter map
							double rate=0; bool usedParam = false; string rateValueParameterName = "";
							try {
								rate = NFutil::convertToDouble(rateValue);
							} catch (std::runtime_error &e1) {
								if(parameter.find(rateValue)==parameter.end()) {
									cerr<<"Could not find parameter: "<<rateValue<<" when reading rate for rxn "<<rxnName <<". Quitting"<<endl;
									return false;
								}
								rate = parameter.find(rateValue)->second;
								usedParam = true;
								rateValueParameterName = rateValue;
							}
							if(verbose) {
								cout<<"\t\t\t\t...setting elementary rate to be: "<<rateValue;
								if(usedParam) cout<<" (has value: "<<rate<<")";
								cout<<endl;
							}

							r->setBaseRate(rate,rateValueParameterName);
						}

						pRateConstant = pRateConstant->NextSiblingElement("RateConstant");
						if(pRateConstant!=NULL) {
							cout<<"\n\n!!Warning:: Multiple RateConstant tags present for RateLaw definition of "<<rxnName<<"."<<endl;
							cout<<"Only the first RateConstant given will be used..."<<endl;
						}
					}
					else if(rateLawType=="Function") {

						//make sure the function has a name
						if(!pRateLaw->Attribute("name")) {
							cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: cannot read 'name' attribute!"<<endl;
							return false;
						}

						//check whether or not we have any arguments to determine if we are
						//a local or global function...
						bool isGlobal = true;
						vector <string> funcArgs;
						TiXmlElement *pListOfArgs = pRateLaw->FirstChildElement("ListOfArguments");
						if(pListOfArgs) {
							TiXmlElement *pArg;
							for ( pArg = pListOfArgs->FirstChildElement("Argument"); pArg != 0; pArg = pArg->NextSiblingElement("Argument"))
							{
								if(!pArg->Attribute("id") || !pArg->Attribute("type") || !pArg->Attribute("value") ) {
									cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: Argument tag \n";
									cerr<<"must have id, type and value attributes defined!"<<endl;
									return false;
								}
								string argId = pArg->Attribute("id");
								string argType = pArg->Attribute("type");
								string argValue = pArg->Attribute("value");
								funcArgs.push_back(argId);

								//cout<<"found argument:"<<argId<<" of type "<<argType<<" which points to "<<argValue<<endl;
								isGlobal=false;

								// Let's find out what object this argument points to..
								if ( argValue.find("_M")!=string::npos ){
									// this appears to be a molecule reference
									// ..make sure the molecule reference is good
									if( reactants.find(argValue)!=reactants.end() ) {
										// found good molecule reference!
										ts->addLocalFunctionReference(reactants.find(argValue)->second,argId,LocalFunction::MOLECULE);
										if(verbose) {cout<<"\t\t\t\tScope of tag "<< argId <<" is MOLECULE (argValue="<<argValue<<")"<<endl; }
									}
									else {
										cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function:\n"
												<<" cannot find referent "<<argId<<" for local function tag "<<argValue<<"."<<endl;
										return false;
									}
								}
								else if ( argValue.find("_RP")!=string::npos ) {
									// this appears to be a pattern ("species") reference
									// ..make sure the pattern reference is good
									if ( comps.find(argValue)!=comps.end() ) {
										// found good pattern reference!
										component c = comps.find(argValue)->second;
										ts->addLocalFunctionReference(c.t,argId,LocalFunction::SPECIES);
										if(verbose) {cout<<"\t\t\t\tScope of tag "<<argId<<" is SPECIES  (argValue="<<argValue<<")"<<endl; }
									}
									else {
										cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function:\n"
												<<" cannot find referent "<<argId<<" for local function tag "<<argValue<<"."<<endl;
										return false;
									}
								}

							}
						}

						//cout<<"found args (in vector):"<<endl;
						//for(int i=0; i<funcArgs.size(); i++) {
						//	cout<<funcArgs.at(i)<<endl;
						//}

						if(isGlobal)
						{
							string functionName = pRateLaw->Attribute("name");
							GlobalFunction *gf = s->getGlobalFunctionByName(functionName);
							if(gf!=NULL) {
								ts->finalize();
								r=new FunctionalRxnClass(rxnName,gf,ts,s);
							} else {
								CompositeFunction *cf = s->getCompositeFunctionByName(functionName);
								if(cf==NULL) {
									cerr<<"When parsing reaction: '"<<rxnName<<"', could not identify function: '"<<functionName<<"' in\n";
									cerr<<"the system.  Therefore, I must abort."<<endl;
									return false;
								}
								ts->finalize();
								r=new FunctionalRxnClass(rxnName,cf,ts,s);
							}
						} else {

							string functionName = pRateLaw->Attribute("name");
							LocalFunction *lf = s->getLocalFunctionByName(functionName);
							if(lf!=NULL) {
								cout<<"Error!! call a local function through a composite function always!"<<endl;
								cout<<"DOR rxn should never directly call a local function."<<endl;
								exit(1);
							} else {

								ts->finalize();

								CompositeFunction *cf = s->getCompositeFunctionByName(functionName);

								r=new DORRxnClass(rxnName,1,"",ts,cf,funcArgs,s);
							}
						}
					}
					else if(rateLawType=="FunctionProduct") {

						//make sure the function1 has a name
						if(!pRateLaw->Attribute("name1")) {
							cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: cannot read 'name' attribute!"<<endl;
							return false;
						}

						//make sure the function2 has a name
						if(!pRateLaw->Attribute("name2")) {
							cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: cannot read 'name' attribute!"<<endl;
							return false;
						}

						//find argument1
						vector <string> funcArgs1;
						TiXmlElement *pListOfArgs1 = pRateLaw->FirstChildElement("ListOfArguments1");
						if(pListOfArgs1)
						{
							TiXmlElement *pArg;
							for ( pArg = pListOfArgs1->FirstChildElement("Argument"); pArg != 0; pArg = pArg->NextSiblingElement("Argument"))
							{
								if(!pArg->Attribute("id") || !pArg->Attribute("type") || !pArg->Attribute("value") ) {
									cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: Argument tag \n";
									cerr<<"must have id, type and value attributes defined!"<<endl;
									return false;
								}
								string argId = pArg->Attribute("id");
								string argType = pArg->Attribute("type");
								string argValue = pArg->Attribute("value");
								funcArgs1.push_back(argId);

								//Get the template molecule we are referring to
								//add a reference to it with the given name
								if(comps.find(argValue)!=comps.end()){
									component c = comps.find(argValue)->second;
									ts->addLocalFunctionReference(c.t,argId,LocalFunction::SPECIES);
									if(verbose) {cout<<"\t\t\t\tScope is SPECIES"<<endl; }
									//exit(1); // set all local functions to allow species scope here!!
								}

								if(reactants.find(argValue)!=reactants.end())
								{
									ts->addLocalFunctionReference(reactants.find(argValue)->second,argId,LocalFunction::MOLECULE);
									if(verbose) {cout<<"\t\t\t\tScope is MOLECULE"<<endl; }
								}
							}
						}

						//find argument2
						vector <string> funcArgs2;
						TiXmlElement *pListOfArgs2 = pRateLaw->FirstChildElement("ListOfArguments2");
						if(pListOfArgs2)
						{
							TiXmlElement *pArg;
							for ( pArg = pListOfArgs2->FirstChildElement("Argument"); pArg != 0; pArg = pArg->NextSiblingElement("Argument"))
							{
								if(!pArg->Attribute("id") || !pArg->Attribute("type") || !pArg->Attribute("value") ) {
									cerr<<"!!Error:: ReactionRule "<<rxnName<<" rate law specification Function: Argument tag \n";
									cerr<<"must have id, type and value attributes defined!"<<endl;
									return false;
								}
								string argId = pArg->Attribute("id");
								string argType = pArg->Attribute("type");
								string argValue = pArg->Attribute("value");
								funcArgs2.push_back(argId);

								//Get the template molecule we are referring to
								//add a reference to it with the given name
								if(comps.find(argValue)!=comps.end()){
									component c = comps.find(argValue)->second;
									ts->addLocalFunctionReference(c.t,argId,LocalFunction::SPECIES);
									if(verbose) {cout<<"\t\t\t\tScope is SPECIES"<<endl; }
									//exit(1); // set all local functions to allow species scope here!!
								}

								if(reactants.find(argValue)!=reactants.end())
								{
									ts->addLocalFunctionReference(reactants.find(argValue)->second,argId,LocalFunction::MOLECULE);
									if(verbose) {cout<<"\t\t\t\tScope is MOLECULE"<<endl; }
								}
							}
						}

						if ( funcArgs1.size() < 1 )
						{
							cout<<"Error!! FunctionProduct ratelaw is missing argument for function1."<<endl;
							exit(1);
						}

						if ( funcArgs2.size() < 1 )
						{
							cout<<"Error!! FunctionProduct ratelaw is missing argument for function2."<<endl;
							exit(1);
						}

						//cout << "n_funcArgs1: " << funcArgs1.size() << endl;
						//cout << "n_funcArgs2: " << funcArgs2.size() << endl;

						string functionName1 = pRateLaw->Attribute("name1");
						LocalFunction *lf1 = s->getLocalFunctionByName(functionName1);
						if(lf1 != NULL) {
							cout<<"Error!! call a local function through a composite function always!"<<endl;
							cout<<"DOR rxn should never directly call a local function."<<endl;
							exit(1);
						}

						string functionName2 = pRateLaw->Attribute("name2");
						LocalFunction *lf2 = s->getLocalFunctionByName(functionName2);
						if(lf2 != NULL) {
							cout<<"Error!! call a local function through a composite function always!"<<endl;
							cout<<"DOR rxn should never directly call a local function."<<endl;
							exit(1);
						}

						ts->finalize();
						CompositeFunction *cf1 = s->getCompositeFunctionByName(functionName1);
						CompositeFunction *cf2 = s->getCompositeFunctionByName(functionName2);

						r=new DOR2RxnClass(rxnName,1,"",ts,cf1,cf2,funcArgs1,funcArgs2,s);

					}
					else if(rateLawType=="MM") {
						//Make sure that the rate constant exists
						TiXmlElement *pListOfRateConstants = pRateLaw->FirstChildElement("ListOfRateConstants");
						if(!pListOfRateConstants) {
							cerr<<"Michaelis-Menten Rate Law definition for "<<rxnName<<" does not have a ListOfRateConstants tag!  Quiting"<<endl;
							return false;
						} else {
							double kcat=0;
							double Km=0;

							//We should always get the catalytic rate (kcat) first!
							TiXmlElement *pRateConstant = pListOfRateConstants->FirstChildElement("RateConstant");
							if(!pRateConstant) {
								cerr<<"Michaelis-Menten Law definition for "<<rxnName<<" does not have kcat value defined!  Quiting"<<endl;
								return false;
							} else {
								//Get the rate constant value
								string kcatName;
								if(!pRateConstant->Attribute("value")) {
									cerr<<"Michaelis-Menten Rate Law definition for "<<rxnName<<" does not have a valid kcat RateConstant 'value'!  Quiting"<<endl;
									return false;
								} else {
									kcatName = pRateConstant->Attribute("value");
									bool usedParam = false;
									try {
										kcat = NFutil::convertToDouble(kcatName);
									} catch (std::runtime_error &e1) {
										if(parameter.find(kcatName)==parameter.end()) {
											cerr<<"Could not find parameter: "<<kcatName<<" when reading kcat rate for rxn "<<rxnName <<". Quitting"<<endl;
											return false;
										}
										kcat = parameter.find(kcatName)->second;
										usedParam = true;
									}
									if(verbose) {
										cout<<"\t\t\t\t...setting kcat rate to be: "<<kcatName;
										if(usedParam) cout<<" (has value: "<<kcat<<")";
										cout<<endl;
									}
								}
							}

							//Now, do the same thing for the Michaelis constant (Km)
							pRateConstant = pRateConstant->NextSiblingElement("RateConstant");
							if(!pRateConstant) {
								cerr<<"Michaelis-Menten Law definition for "<<rxnName<<" does not have Km value defined!  Quiting"<<endl;
								return false;
							} else {
								//Get the rate constant value
								string KmName;
								if(!pRateConstant->Attribute("value")) {
									cerr<<"Michaelis-Menten Rate Law definition for "<<rxnName<<" does not have a valid Km RateConstant 'value'!  Quiting"<<endl;
									return false;
								} else {
									KmName = pRateConstant->Attribute("value");
									bool usedParam = false;
									try {
										Km = NFutil::convertToDouble(KmName);
									} catch (std::runtime_error &e1) {
										if(parameter.find(KmName)==parameter.end()) {
											cerr<<"Could not find parameter: "<<KmName<<" when reading Km rate for rxn "<<rxnName <<". Quitting"<<endl;
											return false;
										}
										Km = parameter.find(KmName)->second;
										usedParam = true;
									}
									if(verbose) {
										cout<<"\t\t\t\t...setting Km rate to be: "<<KmName;
										if(usedParam) cout<<" (has value: "<<Km<<")";
										cout<<endl;
									}
								}
							}

							//Finally, we can make the actual reaction
							ts->finalize();
							r = new MMRxnClass(rxnName,kcat,Km,ts,s);
							
						}
					}
					else if(rateLawType=="Sat") {
						cerr<<"!! Nfsim cannot, and will not, interpret a Rate Law of 'type': 'Sat' as BioNetGen.\n";
						cerr<<"  once could.  You should instead use a Michaelis-Menten rate law (use type 'MM')\n";
						cerr<<"  that also takes the parameters (kcat, Km) but is more accurate.\n"<<endl;
						cerr<<"  But for now, I'm aborting..."<<endl;
						return false;
					}

					////  To extend NFsim to parse more rate law types, add an extra else if clause here to catch the rate law
					////  by its name.

					else {
						cerr<<"!! I cannot yet interpret a Rate Law of 'type': "<<rateLawType<<".\n";
						cerr<<"  Remember, an elementry reaction needs to be specified as 'Ele' (case sensitive),"<<endl;
						cerr<<"  a Michaelis-Menten reaction should be specified as 'MM', and"<<endl;
						cerr<<"  a functional reaction should be specified as 'Functional'."<<endl;
						cerr<<"Quitting."<<endl;
					//	return false;
					}
				} // end to the else statement that parses the rate law.

				if(r==0) {
					cout<<"\n!! Warning!! Unable to create a reaction for some reason!!\n\n"<<endl;
				} else {
					//Finally, add the completed rxn rule to the system
					s->addReaction(r);
					r->setTotalRateFlag(totalRateFlag);
					comps.clear();
				}

			} //end loop through all permutations

		} //end loop through all reaction rules

		//If we got here, then by golly, I think we have a new reaction rule
		return true;

	} catch (...) {
		cout<<"caught something when creating reaction rules.."<<endl;
		return false;
	}

	return false;
}


bool NFinput::readObservableForTemplateMolecules(TiXmlElement *pObs,
		string observableName,
		vector <TemplateMolecule *> &tmList,
		vector <string> &stochRelation,
		vector <int> &stochQuantity,
		System *s,
		map <string,double> &parameter,
		map<string,int> &allowedStates,
		int obsType,
		bool verbose,
		int &suggestedTraversalLimit) {


	TiXmlElement *pListOfPatterns = pObs->FirstChildElement("ListOfPatterns");
	if(!pListOfPatterns) {
		cout<<"\n\n!! Warning:: Observable "<<observableName<<" contains no patterns!"<<endl;
		return true;
	}


	map <string, TemplateMolecule *> templates;

	//Read the pattern for symmetry
	TiXmlElement *pPattern;
	for ( pPattern = pListOfPatterns->FirstChildElement("Pattern"); pPattern != 0; pPattern = pPattern->NextSiblingElement("Pattern"))
	{
		const char *patternName = pPattern->Attribute("id");
		if(!patternName) {
			cerr<<"Pattern tag in observable "<<observableName<<" without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		}

		string relation = "";
		int quantity = -1;

		if(pPattern->Attribute("relation")) {
			if(pPattern->Attribute("quantity")) {
				relation = pPattern->Attribute("relation");
				try {
				quantity = NFutil::convertToInt(pPattern->Attribute("quantity"));
				} catch (std::runtime_error e) {
					cerr<<"I couldn't parse the 'quantity' value when creating pattern: "<<patternName<<endl;
					cerr<<e.what()<<endl;
					return false;
				}
			} else {
				cerr<<"Error when creating observable: "<<observableName<<": relation attribute for stoichiometric\n";
				cerr<<"observable was found without a matching quantity attribute."<<endl;
				return false;
			}
		}

		//cout<<"   --- reading pattern "<<patternName<<" for symmetry"<<endl;
		TiXmlElement *pListOfMols = pPattern->FirstChildElement("ListOfMolecules");
		if(pListOfMols) {

			//We only need to create a new template molecule for each case of symmetry
			//if the observable is a Molecules type
			if(obsType==Observable::MOLECULES)
			{
				map <string, component> comps;
				map <string, component> symComps;
				map <string, component> symMap;

				//First, read the pattern for symmetry - if symmetries exist
				if(!readPatternForSymmetry(pListOfMols, s, patternName, comps, symComps, verbose)) return false;

				//Get the permutations that were created from the symmetries
				vector <map<string,component> > permutations;
				if(!generateRxnPermutations(permutations,symComps,symComps,verbose)) return false;

				//For each valid permutation, create a template molecule that can match it
				for( unsigned int p=0; p<permutations.size(); p++)
				{
					map <string,component> symMap = permutations.at(p);
					map <string,TemplateMolecule *> allTemplatesMap;
					TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, allowedStates, patternName, allTemplatesMap, comps, symMap, verbose, suggestedTraversalLimit);
					if(tm==NULL) return false;
					tmList.push_back(tm);

					if(!relation.empty()) {
						cerr<<"Error when creating observable: "<<observableName<<": a stoichiometric observable found for\n";
						cerr<<"observable of type Molecules.  Currently, NFsim only handles stoichiometric Species\n";
						cerr<<"observables for effeciency.  Rewrite this observable as type 'Species'."<<endl;
						return false;
					}
				}
			}

			//Otherwise, we only have to match once, so we do this for species
			else if(obsType==Observable::SPECIES) {

				map <string,TemplateMolecule *> allTemplatesMap;
				map <string, component> comps;
				map <string, component> symMap;

				TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, allowedStates, patternName, allTemplatesMap, comps, symMap, verbose, suggestedTraversalLimit);
				if(tm==NULL) return false;
				tmList.push_back(tm);
				stochRelation.push_back(relation);
				stochQuantity.push_back(quantity);
				if(verbose && !relation.empty()) {
					cout<<"\t\t\t\t\tHas stoichiometric constraint: '"<<relation<<" "<<quantity<<"'"<<endl;
				}
			}

			//if the observable was not a valid type, we don't read anything
			else {}

		}
		else {
			cerr<<"Observable pattern in observable "<<observableName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
			return false;
		}
	}


	return true;
}


bool NFinput::initObservables(
		TiXmlElement * pListOfObservables,
		System * s,
		map <string,double> &parameter,
		map<string,int> &allowedStates,
		bool verbose,
		int &suggestedTraversalLimit)
{
	try {
		//We will parse this in a similar manner to parsing species, except to say that we don't create
		//actual molecules, just template molecules.

		//First, we loop over all instances of the <Observable> tag
		TiXmlElement *pObs;
		for ( pObs = pListOfObservables->FirstChildElement("Observable"); pObs != 0; pObs = pObs->NextSiblingElement("Observable"))
		{

			string observableId="", observableName="", observableType="";

			//First get the observable name, id, and type
			if(!pObs->Attribute("id") || !pObs->Attribute("name") || !pObs->Attribute("type")) {
				cerr<<"Observable tag without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			} else {
				observableId = pObs->Attribute("id");
				observableName = pObs->Attribute("name");
				observableType = pObs->Attribute("type");
			}
			if(verbose) cout<<"\t\tCreating Observable: '"<<observableName<<"' of type: '"<<observableType<<"'"<<endl;

			//Depending on the type of observable, we have to create it in a particular way
			NFutil::trim(observableType);
			if(observableType.compare("Molecules")==0)
			{
				//First, handle Molecules observables, which are counted for every instance
				//of the given pattern throughout the entire system (so when using identical component
				//names, a single molecule can match the pattern in potentially many ways.  This is
				//handled by cloning each of the ways that a molecule can match as a different template
				//molecule

				//First, read in the pattern as a list of template molecules, which creates the needed
				//symmetric templateMolecules as mentioned above
				vector <TemplateMolecule *> tmList;
				vector <string> stochRelation;
				vector <int> stochQuantity;
				if(!readObservableForTemplateMolecules(pObs,observableName,tmList,stochRelation,stochQuantity,s,parameter,allowedStates,Observable::MOLECULES,verbose,suggestedTraversalLimit)) {return false;}

				//Create the observable, which in this case, is a MoleculesObservable
				MoleculesObservable *mo = new MoleculesObservable(observableName,tmList);

				//Add the observable to each molecule type that will have to check in with this observable
				//Generally, there is just one - but if there are multiple patterns, then we have to match
				//each one separately...
				vector <MoleculeType *> addedMolTypes;
				for(unsigned int k=0; k<tmList.size(); k++) {
					bool alreadyAdded = false;
					for(unsigned int mtc=0; mtc<addedMolTypes.size(); mtc++) {
						if(addedMolTypes.at(mtc)==tmList.at(k)->getMoleculeType()) {
							alreadyAdded = true;
							break;
						}
					}
					if(alreadyAdded) continue;
					tmList.at(k)->getMoleculeType()->addMolObs(mo);
					addedMolTypes.push_back(tmList.at(k)->getMoleculeType());
				}

				//Finally, add the observable to the system so that we can keep track of it for output
				s->addObservableForOutput(mo);

			}
			else if(observableType.compare("Species")==0) {

				cout<<"\nWarning!! Creating a Species Observable!  NFsim does not keep track of Species"<<endl;
				cout<<"by default, so this type of observable requires complex bookkeeping to be turned on"<<endl;
				cout<<"which makes NFsim run slower.  Be sure that you need this observable!\n"<<endl;

				if(!s->isUsingComplex()) {
					cerr<<"You have not turned on complex book keeping, so I can not run you simulation"<<endl;
					cerr<<"because your simulation has a Species observable."<<endl;
					cerr<<"Rerun NFsim with the -cb flag."<<endl;
					return false;
				}

				//First, read in the pattern as a list of template molecules, which creates the needed
				//symmetric templateMolecules as mentioned above
				vector <TemplateMolecule *> tmList;
				vector <string> stochRelation;
				vector <int> stochQuantity;
				if(!readObservableForTemplateMolecules(pObs,observableName,tmList,stochRelation,stochQuantity,s,parameter,allowedStates,Observable::SPECIES,verbose, suggestedTraversalLimit)) {return false;}

				SpeciesObservable *so = new SpeciesObservable(observableName,tmList,stochRelation,stochQuantity);
				s->addObservableForOutput(so);

			}
			else
			{
				cerr<<"Cannot create an Observable of type '"<<observableType<<"'."<<endl;
				cerr<<"The only valid types in NFsim are: 'Molecules' and 'Species', case sensitive."<<endl;
				return false;
			}

		}

		//Getting here means success!
		return true;

	} catch (...) {
		//oh no, what happened this time?
		return false;
	}

	return false;
}









TemplateMolecule *NFinput::readPattern(
		TiXmlElement * pListOfMol,
		System * s,
		map <string,double> &parameter,
		map <string,int> &allowedStates,
		string patternName,
		map <string , TemplateMolecule *> &templates,
		map <string, component> &comps,
		map <string, component> &symMap,
		bool verbose,
		int &suggestedTraversalLimit)
{
	try {

		//A vector to hold each template molecule as we compose the entire template molecule
		vector < TemplateMolecule * > tMolecules;

		//Maps that map binding site ids into a molecule location in the molecules vector
		//and the name of the binding site
		map <string, string> bSiteSiteMapping;
		map <string, int> bSiteMolMapping;

		//vectors that keep track of the states and thier specified values as we create the templates
		vector <string> stateName;
		vector <double> stateValue;

		vector <string> emptyBondSite;
		vector <string> occupiedBondSite;

		//An iterator
		vector <string>::iterator strVecIter;


		// Now loop through the molecules in the list
		TiXmlElement *pMol;
		for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule"))
		{
			//First get the type of molecule and retrieve the moleculeType object from the system
			string molName, molUid;
			if(!pMol->Attribute("name") || ! pMol->Attribute("id")) {
				cerr<<"!!!Error.  Invalid 'Molecule' tag found when creating pattern '"<<patternName<<"'. Quitting"<<endl;
				return NULL;
			} else {
				molName = pMol->Attribute("name");
				molUid = pMol->Attribute("id");
			}


			//Generate an error for any null or trash molecule
			if(molName=="Null" || molName=="NULL" || molName=="null") {
				cerr<<"\n\nError!  You cannot include a 'null' molecule in any reactant or observable pattern!"<<endl;
				cerr<<"Null is a keyword in NFsim used for degradation, so cannot be used.  Check your"<<endl;
				cerr<<"rules and patterns and remove any usage of 'Null' or 'null' or 'NULL'."<<endl;
				exit(1);
			}
			if(molName=="Trash" || molName=="TRASH" || molName=="trash") {
				cerr<<"\n\nError!  You cannot include a 'trash' molecule in any reactant or observable pattern!"<<endl;
			    cerr<<"Trash is a keyword in NFsim used for degradation, so cannot be used.  Check your"<<endl;
			    cerr<<"rules and patterns and remove any usage of 'Trash' or 'trash' or 'TRASH'."<<endl;
				exit(1);
			}


			//Get the moleculeType and create the actual template
			MoleculeType *moltype = s->getMoleculeTypeByName(molName);
			TemplateMolecule *tempmol = new TemplateMolecule(moltype);
			if(verbose) cout<<"\t\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;

			//Create a comp element that matches onto this molecule so we can retrieve
			//any pointers to this molecule (for instance, to allow deletion of this molecule)
			component c(tempmol,"");
			comps.insert(pair <string, component> (molUid,c));


			//Loop through the components of the molecule in order to set state values
			TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
			if(pListOfComp)
			{
				TiXmlElement *pComp;
				for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component"))
				{
					//Get the basic components of this molecule
					string compId, compName, compBondCount;
					if(!pComp->Attribute("id") || !pComp->Attribute("name") || !pComp->Attribute("numberOfBonds")) {
						cerr<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
						return NULL;
					} else {
						compId = pComp->Attribute("id");
						compName = pComp->Attribute("name");
						compBondCount = pComp->Attribute("numberOfBonds");
					}

					//Set up basic component info so we can go back to it if need be...
					component c(tempmol,compName);
					comps.insert(pair <string, component> (compId,c));


					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////

					// For debugging: does this component actually exist?

//					cout<<"Here is the sym Map: "<<endl;
//					map <string,component>::iterator mapIter;
//					for(mapIter=symMap.begin();mapIter!=symMap.end(); mapIter++) {
//						cout<<mapIter->first<<"   "<<mapIter->second.name<<"  "<<mapIter->second.uniqueId<<endl;
//					}
//					component *symC;
//					cout<<compId<<endl;


					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////
					// Handle equivalent components off reaction center differently
					// it is off reaction center if 1) it is an eq component and 2) it is not in the symMap
					if(symMap.find(compId)==symMap.end() && moltype->isEquivalentComponent(compName)) {
						//cout<<"we should treat this as a symmetric constraint"<<endl;
						//cout<<"equivalent type! :"<<compName<<endl;
						int stateConstraint = -1;
						if(pComp->Attribute("state")) {
							string compStateValue = pComp->Attribute("state");
							if(compStateValue!="*" && compStateValue!="?") {
								if(allowedStates.find(molName+"_"+compName+"_"+compStateValue)==allowedStates.end()) {
									cerr<<"You are trying to give a pattern of type '"<<molName<<"', but you gave an "<<endl;
									cerr<<"invalid state! The state you gave was: '"<<compStateValue<<"'.  Quitting now."<<endl;
									return false;
								} else {
									//State is a valid allowed state, so push it onto our list
									stateConstraint = allowedStates.find(molName+"_"+compName+"_"+compStateValue)->second;
								}
							}
						}

						//Check it as a binding site...
						int bondConstraint = TemplateMolecule::NO_CONSTRAINT;
						if(pComp->Attribute("numberOfBonds")) {
							string numOfBonds = pComp->Attribute("numberOfBonds");
							int numOfBondsInt = -1;

							//const int MUST_BE_OCCUPIED = -2;
							//const int EITHER_WAY_WORKS = -3;
							if(numOfBonds.compare("+")==0) {
								bondConstraint = TemplateMolecule::OCCUPIED;
							} else if(numOfBonds.compare("*")==0) {
								bondConstraint = TemplateMolecule::NO_CONSTRAINT;
							} else if(numOfBonds.compare("?")==0) {
								bondConstraint = TemplateMolecule::NO_CONSTRAINT;
							} else {
								try {
									numOfBondsInt = NFutil::convertToInt(numOfBonds);
								} catch (std::runtime_error e) {
									cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
									cerr<<e.what()<<endl;
									return false;
								}

								if(numOfBondsInt==0) {
									bondConstraint = TemplateMolecule::EMPTY;
								} else if (numOfBondsInt==1) {
									bondConstraint = TemplateMolecule::OCCUPIED;
									bSiteSiteMapping[compId] = compName;
									bSiteMolMapping[compId] = tMolecules.size();
								} else {
									cerr<<"I can only handle a site that has 0 or 1 bonds in pattern: "<<patternName<<endl;
									cerr<<"You gave me "<<numOfBondsInt<<" instead for component "<<compName<<endl;
									return false;
								}
							}
						}
						tempmol->addSymCompConstraint(compName,compId,bondConstraint,stateConstraint);


					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////
					} else {
						//cout<<"can be mapped as a symmetric component, or is not a symmetric component"<<endl;
						//Read in a state, if it is in fact has an associated state
						if(pComp->Attribute("state")) {
							string compStateValue = pComp->Attribute("state");

							//Make sure the given state is allowed (we allow for wildcards...)
							if(compStateValue!="*" && compStateValue!="?") {
								if(allowedStates.find(molName+"_"+compName+"_"+compStateValue)==allowedStates.end()) {
									cerr<<"You are trying to give a pattern of type '"<<molName<<"', but you gave an "<<endl;
									cerr<<"invalid state! The state you gave was: '"<<compStateValue<<"'.  Quitting now."<<endl;
									return false;
								} else {
									//State is a valid allowed state, so push it onto our list
									int stateValueInt = allowedStates.find(molName+"_"+compName+"_"+compStateValue)->second;

									//Make sure we catch symmetric components in the case of a state change...
									component *symC;  //cout<<"state value: "<< compId;
									if(!lookup(symC, compId, comps, symMap)) {
										cerr<<"Could not find the symmetric component when creating a component state, but there\n";
										cerr<<"should have been one!!  I don't know what to do, so I'll quit."<<endl;
										return false;
									}
									stateName.push_back(symC->symPermutationName);
									stateValue.push_back(stateValueInt);
								}
							}
						}

						//Check it as a binding site...
						if(pComp->Attribute("numberOfBonds")) {
							string numOfBonds = pComp->Attribute("numberOfBonds");
							int numOfBondsInt = -1;

							//Only try to parse this bond as a number if the number
							//of bonds is not the '+' character.  The '+' character implies
							//that the site is occupied without explicitly specifying who it is
							//bound to.  Now also handles the wild card character (implying that
							//it can be bound or unbound - it doesn't matter.
							const int MUST_BE_OCCUPIED = -2;
							const int EITHER_WAY_WORKS = -3;
							if(numOfBonds.compare("+")==0) {
								numOfBondsInt = MUST_BE_OCCUPIED;
							} else if(numOfBonds.compare("*")==0) {
								numOfBondsInt = EITHER_WAY_WORKS;
							} else if(numOfBonds.compare("?")==0) {
								numOfBondsInt = EITHER_WAY_WORKS;
							} else {
								try {
									numOfBondsInt = NFutil::convertToInt(numOfBonds);
								} catch (std::runtime_error e) {
									cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
									cerr<<e.what()<<endl;
									return false;
								}
							}

							//cout<<"bond value: "<< compId <<"    -  " <<numOfBondsInt<<"\n";

							//Look up this site in case we have some symmetry going on...
							component *symC;
							if(!lookup(symC, compId, comps, symMap)) {
								cerr<<"Could not find the symmetric component when creating a binding site, but there\n";
								cerr<<"should have been one!!  I don't know what to do, so I'll quit."<<endl;
								return false;
							}

							if(numOfBondsInt==MUST_BE_OCCUPIED) {
								occupiedBondSite.push_back(symC->symPermutationName);
							} else if(numOfBondsInt==EITHER_WAY_WORKS) {
								//add nothing if either way works of course!  There
								//is no constraint (the two ways are either bonded or not bonded)
							} else if(numOfBondsInt==0) {
								emptyBondSite.push_back(symC->symPermutationName);
							} else if (numOfBondsInt==1) {
								bSiteSiteMapping[compId] = symC->symPermutationName;
								bSiteMolMapping[compId] = tMolecules.size();
							} else {
								cerr<<"I can only handle a site that has 0 or 1 bonds in pattern: "<<patternName<<endl;
								cerr<<"You gave me "<<numOfBondsInt<<" instead for component "<<compName<<endl;
								return false;
							}
						}

					} // end if statement over symmetric components
					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////
					//////////////////////////////////////////////////////

				} //end loop over components
			} //end if statement for compenents to exist


			//Loop through the states and set the constraints we need to set
			int k=0;
			for(strVecIter = stateName.begin(); strVecIter != stateName.end(); k++, strVecIter++ )
			{
				tempmol->addComponentConstraint(*strVecIter,(int)stateValue.at(k));
			}
			for(strVecIter = emptyBondSite.begin(); strVecIter != emptyBondSite.end(); strVecIter++ )
			{
				tempmol->addEmptyComponent(*strVecIter);
			}
			for(strVecIter = occupiedBondSite.begin(); strVecIter != occupiedBondSite.end(); strVecIter++ )
			{
				tempmol->addBoundComponent(*strVecIter);
			}



			//tempmol->printDetails();


			//Update our data storage with the new template and empty out the things we don't need
			templates.insert(pair <string, TemplateMolecule *> (molUid,tempmol));
			tMolecules.push_back(tempmol);
			stateName.clear();
			stateValue.clear();
			emptyBondSite.clear();
			occupiedBondSite.clear();
		}

		//Here is where we add the bonds to the template molecules in the pattern
		TiXmlElement *pListOfBonds = pListOfMol->NextSiblingElement("ListOfBonds");
		if(pListOfBonds)
		{
			//First get the information on the bonds in the complex
			TiXmlElement *pBond;
			for ( pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond"))
			{
				//First, grab the bond information that we need to set things up
				string bondId, bSite1, bSite2;
				if(!pBond->Attribute("id") || !pBond->Attribute("site1") || !pBond->Attribute("site2")) {
					cerr<<"!! Invalid Bond tag for pattern: "<<patternName<<".  Quitting."<<endl;
					return false;
				} else {
					bondId = pBond->Attribute("id");
					bSite1 = pBond->Attribute("site1");
					bSite2 = pBond->Attribute("site2");
				}

				//if(verbose)cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;

				//Get the information on this bond that tells us which molecules to connect
				try {

//					cout<<"here"<<endl;
//					cout<<"bSite1: "<<bSite1<<endl;
//					cout<<"bSite2: "<<bSite2<<endl;
//
//					cout<<"bSiteSiteMapping"<<endl;
//					for ( std::map< string, string>::const_iterator iter = bSiteSiteMapping.begin();
//					iter != bSiteSiteMapping.end(); ++iter )
//						cout << iter->first << '\t' << iter->second << '\n';
//					cout<<"bSiteMolMapping"<<endl;
//					for ( std::map< string, string>::const_iterator iter = bSiteSiteMapping.begin();
//					iter != bSiteSiteMapping.end(); ++iter )
//						cout << iter->first << '\t' << iter->second << '\n';



					//First look up the info from the component maps
					if(		bSiteSiteMapping.find(bSite1)!=bSiteSiteMapping.end() &&
							bSiteMolMapping.find(bSite1)!=bSiteMolMapping.end() &&
							bSiteSiteMapping.find(bSite2)!=bSiteSiteMapping.end() &&
							bSiteMolMapping.find(bSite2)!=bSiteMolMapping.end()
							) {

						string bSiteName1 = bSiteSiteMapping.find(bSite1)->second;
						int bSiteMolIndex1 = bSiteMolMapping.find(bSite1)->second;
						string bSiteName2 = bSiteSiteMapping.find(bSite2)->second;
						int bSiteMolIndex2 = bSiteMolMapping.find(bSite2)->second;
						TemplateMolecule::bind(tMolecules.at(bSiteMolIndex1),bSiteName1.c_str(),bSite1,
								tMolecules.at(bSiteMolIndex2),bSiteName2.c_str(),bSite2);

						//Erase the bonds to make sure we don't add them again
						bSiteSiteMapping.erase(bSite1);
						bSiteMolMapping.erase(bSite1);
						bSiteSiteMapping.erase(bSite2);
						bSiteMolMapping.erase(bSite2);
					} else {

						cerr<<"here"<<endl;
						cerr<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating pattern '"<<patternName<<"'. "<<endl;
						cerr<<"This may be caused because you are adding two bonds to one binding site or because you listed"<<endl;
						cerr<<"a binding site at the end of the pattern that does not exist.  Quitting"<<endl;
						return false;
					}
				} catch (exception& e) {

					cerr<<e.what()<<endl;

					cerr<<"now here"<<endl;

					cerr<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating pattern '"<<patternName<<"'. "<<endl;
					cerr<<"This may be caused because you are adding two bonds to one binding site or because you listed"<<endl;
					cerr<<"a binding site at the end of the pattern that does not exist.  Quitting"<<endl;
					return false;
				}
			}
		}

		//Now, we have to loop through all the bonds that we did not explicitly use, but that we set to bonded
		//this must be a constraint on our pattern.
		try {
			map<string,string>::iterator strMapit;
			for(strMapit = bSiteSiteMapping.begin(); strMapit != bSiteSiteMapping.end(); strMapit++ ) {
				string bSiteId = (*strMapit).first;
				string bSiteName = (*strMapit).second;
				int bSiteMolIndex = bSiteMolMapping.find(bSiteId)->second;

				//skip symmetric sites here, because we already considered them earlier...
				if(tMolecules.at(bSiteMolIndex)->getMoleculeType()->isEquivalentComponent(bSiteName.c_str())) {
					continue;
				}

				tMolecules.at(bSiteMolIndex)->addBoundComponent(bSiteName);
			}
		} catch (exception& e) {
			cerr<<"Something went wacky when I was parsing pattern '"<<patternName<<"'."<<endl;
			cerr<<"It happened as I was trying to add an occupied binding site to a template molecule."<<endl;
			return false;
		}


		//Print out all the templates we made...
		//for(int k=0; k<tMolecules.size(); k++) {
		//	tMolecules.at(k)->printDetails(cout);
		//}


		//Now we have to find disjointed sets - that is whenever we have a Template
		//Molecule that is connected, but not explicitly through bonds, we have to
		//connect them via the connectedTo specification


		//cout<<"checking for disjoint sets..."<<endl;
		vector <vector <TemplateMolecule *> > sets;
		vector <int> uniqueSetId;
		int setCount = TemplateMolecule::getNumDisjointSets(tMolecules,sets,uniqueSetId);

//		cout<<"Unique Set Ids for the templates: "<<endl;
//		for(unsigned int i=0; i<uniqueSetId.size(); i++) {
//			cout<<uniqueSetId.at(i)<<endl;
//		}


		if(setCount>1) {
			// Possibly, we might want to enforce complex bookkeeping for such reactions....
			//if(!s->isUsingComplex()) {
			//	cout.flush();
			//	cerr<<"Disjoint pattern found, but complex bookkeeping is turned off!"<<endl;
			//	cerr<<"Rerun with the -cb flag"<<endl;
			//	exit(1);
			//}

			cout<<"\nFound disjoint sets in a pattern. (As in A().B(), with no explicit connection through components)\n";
			cout<<"Warning!  These type of patterns can be dangerous!!  They also make NFsim run slower!\n";
			cout<<"If you can express this pattern without this syntax, it is highly advised!\n"<<endl;

			//Add the connected-to connections
			int tm1=0; int tm2=0;

			//connect them in order, 0 to 1, then 1 to 2, then 2 to 3...
			for(int cSet=0; cSet<(setCount-1); cSet++) {

				//cout<<"Matching up set: "<<cSet<<" to "<<cSet+1<<endl;
				for(unsigned int i=0; i<tMolecules.size(); i++) {
					if(uniqueSetId.at(i)==cSet) { tm1=i; break; }
				}
				for(unsigned int i=0; i<tMolecules.size(); i++) {
					if(uniqueSetId.at(i)==(cSet+1)) { tm2=i; break; }
				}
				int ctIndex1=tMolecules.at(tm1)->getN_connectedTo();
				int ctIndex2=tMolecules.at(tm2)->getN_connectedTo();
				tMolecules.at(tm1)->addConnectedTo(tMolecules.at(tm2),ctIndex2);
				tMolecules.at(tm2)->addConnectedTo(tMolecules.at(tm1),ctIndex1);
			}


			//for(unsigned int i=0; i<tMolecules.size(); i++) {
			//	tMolecules.at(i).addConnectedTo(tMolecules.at())
			//	tMolecules.at(i)->printDetails();
			//}

			//cout<<"traversing...  let's see if we got everyone:"<<endl;
			//vector <TemplateMolecule *> tmList;
			//TemplateMolecule::traverse(tMolecules.at(1),tmList);
			//for(unsigned int i=0; i<tmList.size(); i++) {
			//	tmList.at(i)->printDetails();
			//}


		}


		//The number of templateMolecules in this pattern will give us the depth of the traversal that we could
		//ever really encounter in the system, so we should suggest that this be the traversal limit if it
		//is higher than what has already been suggested.
		if((int)tMolecules.size()>suggestedTraversalLimit) {
			suggestedTraversalLimit = (int)tMolecules.size();
		}

		//Grab the first template molecule from the list, and arbitrarily set this as the root
		if(tMolecules.empty()){
			cerr<<"You have a pattern named "<<patternName<<" that doesn't include any actual patterns!  (Or I just couldn't find any)"<<endl;
			cerr<<"Therefore, I see no other choice than to quit until you fix the problem."<<endl;
			return false;
		}
		TemplateMolecule *finalTemplate = tMolecules.at(0);

		tMolecules.clear();
		bSiteMolMapping.clear();
		bSiteSiteMapping.clear();


		//Add a pointer to this reactant as a component so that we can get to it on deletes
		component c(finalTemplate, "");
		comps.insert(pair <string, component> (patternName,c));


		if(verbose) { cout<<"\t\t\t\t => Final processed pattern: "<<finalTemplate->getPatternString()<<endl; }
		return finalTemplate;
	} catch (...) {
		//Here's our final catch all!
		return NULL;
	}

	return false;
}


bool NFinput::readProductPattern(
		TiXmlElement * pListOfMol,
		System * s, map <string,double> &parameter,
		map<string,int> &allowedStates,
		string patternName,
		vector <MoleculeType *> &productMoleculeTypes,
		vector < vector <int> > &stateInformation,
		vector < vector <int> > &bindingSiteInformation,
		bool verbose)
{
	//cout<<"reading the product pattern!"<<endl;
	try {

		//Variables to remember the index of things
		const int proMolTypeIndex = 0;
		const int stateIndex = 1, stateValue=2;
		const int bSiteIndex = 1, partnerMolTypeIndex = 2, partnerBsiteIndex = 3;

		//First, set up the stateInformation vector
		{
			//v1 maps index into productMoleculetypes, v2 maps the stateIndex
			//v3 maps the state value
			vector <int> v1,v2,v3;
			stateInformation.push_back(v1);
			stateInformation.push_back(v2);
			stateInformation.push_back(v3);
		}

		//And here we go with the bindingSiteInformation vector
		{
			//v1 maps index into productMoleculetypes, v2 maps the bindingSiteIndex
			//v3 maps the partner's index into productMoleculeTypes, and v4 maps the
			vector <int> v1,v2,v3, v4;
			bindingSiteInformation.push_back(v1);
			bindingSiteInformation.push_back(v2);
			bindingSiteInformation.push_back(v3);
			bindingSiteInformation.push_back(v4);
		}

		//Some mappings to keep track of binding sites
		map <string, string> bSiteSiteMapping;
		map <string, int> bSiteMolMapping;

		//Map the ids of the molecules with thier location in the productMoleculeTypes vector
		map <string,int> uniqueIdMapping;

		// Now loop through the molecules in the list
		TiXmlElement *pMol;
		bool foundTrash = false;
		for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule"))
		{
			//First get the type of molecule and retrieve the moleculeType object from the system
			string molName, molUid;
			if(!pMol->Attribute("name") || ! pMol->Attribute("id")) {
				cerr<<"!!!Error.  Invalid 'Molecule' tag found when creating product pattern '"<<patternName<<"'. Quitting"<<endl;
				return false;
			} else {
				molName = pMol->Attribute("name");
				molUid = pMol->Attribute("id");
			}

			//Skip anything that is a null molecule
			if(molName=="Null" || molName=="NULL" || molName=="null") {
				if(verbose) cout<<"\t\t\t\tskipping a null molecule in product pattern..."<<endl;
				foundTrash=true;
				continue;
			}
			//Skip anything that is trash molecule
			if(molName=="Trash" || molName=="TRASH" || molName=="trash") {
				if(verbose) cout<<"\t\t\t\tskipping a trash molecule in product pattern..."<<endl;
				foundTrash=true;
				continue;
			}

			//Retrieve the moleculeType and remember it
			MoleculeType *moltype = s->getMoleculeTypeByName(molName);
			productMoleculeTypes.push_back(moltype);

			int currentMoleculeTypeIndex = productMoleculeTypes.size()-1;
			uniqueIdMapping.insert(pair<string,int>(molUid,currentMoleculeTypeIndex));

			if(verbose) cout<<"\t\t\t\tIncluding Product Molecule of type: "<<molName<<" with local id: " << molUid<<endl;

			//Loop through the components of the molecule in order to remember state values
			TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
			if(pListOfComp)
			{
				TiXmlElement *pComp;
				for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component"))
				{
					//Get the basic components of this molecule
					string compId, compName, compBondCount;
					if(!pComp->Attribute("id") || !pComp->Attribute("name") || !pComp->Attribute("numberOfBonds")) {
						cerr<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
						return false;
					} else {
						compId = pComp->Attribute("id");
						compName = pComp->Attribute("name");
						compBondCount = pComp->Attribute("numberOfBonds");
					}


					//Read in a state, if it is in fact a state
					if(pComp->Attribute("state")) {
						string compStateValue = pComp->Attribute("state");

						//Make sure the given state is allowed
						if(allowedStates.find(molName+"_"+compName+"_"+compStateValue)==allowedStates.end()) {
							cerr<<"You are trying to give a pattern of type '"<<molName<<"', but you gave an "<<endl;
							cerr<<"invalid state! The state you gave was: '"<<compStateValue<<"'.  Quitting now."<<endl;
							return false;
						} else {
							//State is a valid allowed state, so push it onto our lists to remember it
							int stateValueInt = allowedStates.find(molName+"_"+compName+"_"+compStateValue)->second;
							stateInformation.at(proMolTypeIndex).push_back(currentMoleculeTypeIndex);
							stateInformation.at(stateIndex).push_back(moltype->getCompIndexFromName(compName));
							stateInformation.at(stateValue).push_back(stateValueInt);
						}
					}

					//Otherwise, parse it as a binding site - here we just check to make sure the number of
					//bonds specified is 0 or 1, and remember the site if it was set to 1.  This way we ensure
					//that you can't specify the site has zero bonds, but later try to add a bond there.
					if(pComp->Attribute("numberOfBonds")) {
						string numOfBonds = pComp->Attribute("numberOfBonds");
						int numOfBondsInt = -1;
						//string bondConstraint;

						//const int MUST_BE_OCCUPIED = -2;
						//const int EITHER_WAY_WORKS = -3;
						if(numOfBonds.compare("+")==0) {
							//bondConstraint = TemplateMolecule::OCCUPIED;
							//numOfBondsInt=1;
						} else if(numOfBonds.compare("*")==0) {
							//bondConstraint = TemplateMolecule::NO_CONSTRAINT;
						} else if(numOfBonds.compare("?")==0) {
							//bondConstraint = TemplateMolecule::NO_CONSTRAINT;
						} else {
							try {
								numOfBondsInt = NFutil::convertToInt(numOfBonds);
							} catch (std::runtime_error e) {
								cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
								cerr<<e.what()<<endl;
								return false;
							}
							if (numOfBondsInt==1) {
								bSiteSiteMapping[compId] = compName;
								bSiteMolMapping[compId] = currentMoleculeTypeIndex;
							} else if (numOfBondsInt!=0) {
								cerr<<"I can only handle a site that has 0 or 1 bonds in pattern: "<<patternName<<endl;
								cerr<<"You gave me "<<numOfBondsInt<<" instead for component "<<compName<<endl;
								return false;
							}
						}


					}
				} //end loop over components
			} //end if statement for compenents to exist
		} // end loop over molecules



		//Here is where we find the bonds that need to be created
		TiXmlElement *pListOfBonds = pListOfMol->NextSiblingElement("ListOfBonds");
		if(pListOfBonds)
		{
			//First get the information on the bonds in the complex
			TiXmlElement *pBond;
			for ( pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond"))
			{
				//First, grab the bond information that we need to set things up
				string bondId, bSite1, bSite2;
				if(!pBond->Attribute("id") || !pBond->Attribute("site1") || !pBond->Attribute("site2")) {
					cerr<<"!! Invalid Bond tag for pattern: "<<patternName<<".  Quitting."<<endl;
					return false;
				} else {
					bondId = pBond->Attribute("id");
					bSite1 = pBond->Attribute("site1");
					bSite2 = pBond->Attribute("site2");
				}

				//if(verbose)cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;

				//Get the information on this bond that tells us which molecules to connect
				try {

					//First look up the info and add the bond
					string bSiteName1 = bSiteSiteMapping.find(bSite1)->second;
					int bSiteMolTypeIndex1 = bSiteMolMapping.find(bSite1)->second;
					string bSiteName2 = bSiteSiteMapping.find(bSite2)->second;
					int bSiteMolIndex2 = bSiteMolMapping.find(bSite2)->second;

					//Add the information to the list
					bindingSiteInformation.at(proMolTypeIndex).push_back(bSiteMolTypeIndex1);
					unsigned int bSiteIndex1 = productMoleculeTypes.at(bSiteMolTypeIndex1)->getCompIndexFromName(bSiteName1);
					bindingSiteInformation.at(bSiteIndex).push_back(bSiteIndex1);
					bindingSiteInformation.at(partnerMolTypeIndex).push_back(bSiteMolIndex2);
					unsigned int bSiteIndex2 = productMoleculeTypes.at(bSiteMolIndex2)->getCompIndexFromName(bSiteName2);
					bindingSiteInformation.at(partnerBsiteIndex).push_back(bSiteIndex2);

					//Erase the bonds to make sure we don't add them again
					bSiteSiteMapping.erase(bSite1);
					bSiteMolMapping.erase(bSite1);
					bSiteSiteMapping.erase(bSite2);
					bSiteMolMapping.erase(bSite2);
				} catch (exception& e) {
					cerr<<"here we are"<<endl;

					cerr<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating product pattern '"<<patternName<<"'. "<<endl;
					cerr<<"This may be caused because you are adding two bonds to one binding site or because you listed"<<endl;
					cerr<<"a binding site at the end of the pattern that does not exist.  Quitting"<<endl;
					return false;
				}
			}
		}


		//Make sure we actually did something!
		if(productMoleculeTypes.empty()){
			if(foundTrash) {
				if(verbose) cout<<"\t\t\t\tWarning: You have an add molecule rule, but only a Trash or Null pattern listed..."<<endl;
				return true;
			}
			cerr<<"You have a product pattern named "<<patternName<<" that doesn't include any actual patterns!  (Or I just couldn't find any)"<<endl;
			cerr<<"Therefore, I see no other choice than to quit until you fix the problem."<<endl;
			return false;
		}

		bSiteMolMapping.clear();
		bSiteSiteMapping.clear();
		return true;
	} catch (...) {
		//Here's our final catch all!
		return false;
	}

	return false;
}



bool NFinput::readProductMolecule(
		TiXmlElement * pMol,
		System * s,
		map <string,double> & parameter,
		map<string,int> & allowedStates,
		string patternName,
		vector <MoleculeCreator *> & moleculeCreatorsList,
		map <string, component> & comps,
		bool verbose )
{
	//cout<<"reading the product molecule!"<<endl;
	try
	{
		// this will hold the template molecule
		TemplateMolecule * tempmol;
		// this holds a pointer to the molecule type
		MoleculeType     * moltype;
		// this will hold the state information
		vector < pair <int,int> >  component_states;


		//First get the type of molecule and retrieve the moleculeType object from the system
		string molName, molUid;

		if( !pMol->Attribute("name") || !pMol->Attribute("id") )
		{
			cerr << "!!!Error.  Invalid 'Molecule' tag found when creating product molecule '"
			     << patternName << "'. Quitting" << endl;
			return false;
		}

		molName = pMol->Attribute("name");
		molUid  = pMol->Attribute("id");


		//Skip anything that is a null molecule
		if(molName=="Null" || molName=="NULL" || molName=="null")
		{
			if(verbose) cout<<"\t\t\t\treadProductMolecule is ignoring " << molName << " molecule..."<<endl;
			return true;
		}

		//Skip anything that is trash molecule
		if(molName=="Trash" || molName=="TRASH" || molName=="trash")
		{
			if(verbose) cout<<"\t\t\t\treadProductMolecule is ignoring " << molName << " molecule..."<<endl;
			return true;
		}


		// Retrieve the moleculeType and remember it
		moltype = s->getMoleculeTypeByName( molName );

		// Create the TemplateMolecule
		tempmol = new TemplateMolecule( moltype );

		// Create a component object for this molecule
		component c_mol(tempmol,"");
		comps.insert( pair <string, component> (molUid,c_mol) );


	    // Write some info if verbose
		if(verbose) cout << "\t\t\t\tIncluding Product Molecule of type: " << molName
		                 << " with local id: " << molUid << endl;


		// BEGIN building the symmetric component list that we need to consider
		// we have to remember the symmetric sites- the input file is going to provide us
		// with the generic site name, but to actually create the molecule, those site names
		// must be mapped onto a specific site name.  So here, we use the molecule type
		// to create us a vector of queues of names that, when given a generic site component, such
		// as 'r', will return us the next available specific site name, such as 'r3' (if we
		// have used up 'r1' and 'r2' already).  This also allows us to check at the end
		// if we used up all the sites, because if we didn't, then this new molecule has
		// not been completely specified.  --michael

		//step 1: generate the vector of queues.  The way this is set up is the vector is as
		//long as the number of equivalence classes (which is the number of unique site names)
		//and for each type of symmetric site, we list all the individual names of the specific
		//sites.
		int numberOfSymmetricSites = moltype->getNumOfEquivalencyClasses();
		vector < queue <string> > symmetricSitesNotUsed;
		for(int i=0; i<numberOfSymmetricSites; i++)
		{
			queue <string> q;
			symmetricSitesNotUsed.push_back(q);
		}




		//step 2: loop over the molecules components and populate the vector of queues
        for(int k=0; k<moltype->getNumOfComponents(); k++) {

        	string compName = moltype->getComponentName(k);

        	// if the component is symmetric...
        	if(moltype->isEquivalentComponent(k)) {

        		string genericCompName = moltype->getEquivalenceClassComponentNameFromComponentIndex(k);
        		//cout<<"generic compname: "<< genericCompName <<endl;

        		// then we determine the component equivalence class number and add it to the queue
        		int eqClassNumber = moltype->getEquivalenceClassNumber(k);
        		symmetricSitesNotUsed.at(eqClassNumber).push(compName);
        		//cout<<"eq class number: "<< moltype->getEquivalenceClassNumber(k) <<endl;
        	}
        }


		// END building the symmetric component list that we need to consider
        // /////////////////////////////////////////////////////////////////////////////


		//Loop through the components of the molecule in order to remember state values
		TiXmlElement * pListOfComp = pMol->FirstChildElement("ListOfComponents");
		if(pListOfComp)
		{
			TiXmlElement *pComp;
			for ( pComp = pListOfComp->FirstChildElement("Component");
					pComp != 0;  pComp = pComp->NextSiblingElement("Component") )
			{

				//Get the basic components of this molecule (ignore any bonds)
				string compId, compName, compStateValue;
				int componentIndex, stateValueIndex;

				if( !pComp->Attribute("id") || !pComp->Attribute("name") || !pComp->Attribute("numberOfBonds") )
				{
					cerr << "!!!Error.  Invalid 'Component' tag found when creating '" << molUid
					     << "' of pattern '"<< patternName << "'. Quitting." << endl;
					return false;
				}

				compId   = pComp->Attribute("id");
				compName = pComp->Attribute("name");

				// we need to make sure this isn't a symmetric component
				if(moltype->isEquivalentComponent(compName))
				{
					//cout<<"Is EQ component!"<<endl;
					int eqClass = moltype->getEquivalencyClassNumber(compName);
					//cout<<"class: "<<moltype->getEquivalencyClassNumber(compName)<<endl;
					if (symmetricSitesNotUsed.at(eqClass).size()>0) {
						compName = (symmetricSitesNotUsed.at(eqClass)).front();
						(symmetricSitesNotUsed.at(eqClass)).pop();
					} else {
						cerr <<"In NFinput::readProductMolecule(...): some symmetric components were specified more times than they were originally declared!"<<endl;
						cerr <<"This must be an error in BioNetGen, as all molecules to be created must have all component sites specified!\n\n"<<endl;
						return false;
					}
				}

				// Add empty component site to the templateMolecule
				tempmol->addEmptyComponent( compName );

				// Construct a component object for ths site
				component c_site(tempmol,compName);
				comps.insert(pair <string, component> (compId,c_site));

				// Read in a state, if it is in fact a state
				if(pComp->Attribute("state"))
				{
					compStateValue = pComp->Attribute("state");

					//Make sure the given state is allowed
					if( allowedStates.find(molName+"_"+compName+"_"+compStateValue)==allowedStates.end() )
					{
						cerr << "Attempt to assign invalid component state value '" << compStateValue
						     << "' in product molecule '" << molName << "'.  Quitting now." << endl;
						return false;
					}

					// State is a valid allowed state: get component index
					componentIndex  = moltype->getCompIndexFromName(compName);
					// Get state value index
					stateValueIndex = allowedStates.find(molName+"_"+compName+"_"+compStateValue)->second;
					// Push onto component states vector
					component_states.push_back( pair <int,int> (componentIndex, stateValueIndex) );

					// add component state value constraint to the template molecule
					tempmol->addComponentConstraint( compName, compStateValue );
				}

			} //end loop over components

		} //end if statement for compenents to exist


		// //// quick check to see if we used up all symmetric component sites
		for(unsigned int i=0; i<symmetricSitesNotUsed.size(); i++) {
			if(symmetricSitesNotUsed.at(i).size()>0) {
				cerr <<"In NFinput::readProductMolecule(...): some symmetric components exist, but were not defined!"<<endl;
				cerr <<"This must be an error in BioNetGen, as all molecules to be created must have all component sites specified!\n\n"<<endl;
				return false;
			}
		}






		// don't forget:  add molecule creator to list
		moleculeCreatorsList.push_back( new MoleculeCreator( tempmol, moltype, component_states ) );

		// return the templateMolecule
		return true;
	}
	catch (...)
	{
		//Here's our final catch all!
		return false;
	}


	return false;
}










