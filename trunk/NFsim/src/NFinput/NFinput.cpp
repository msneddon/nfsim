#include "NFinput.hh"





using namespace NFinput;
using namespace std;



component::component(TemplateMolecule *t, string name)
{
	this->t=t;
	this->name = name;
}

component::component(MoleculeType *mt, string name)
{
	this->mt=mt;
	this->name = name;
}

component::~component()
{
	t=NULL;
}



System * NFinput::initializeFromXML(
		string filename,
		bool verbose)
{
	if(!verbose) cout<<"reading xml file ("+filename+")  [";
	if(verbose) cout<<"\tTrying to read xml model specification file: '"<<filename<<"'"<<endl;
	
	
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
			s=new System("nameless");
			if(verbose) cout<<"\tNo System name given, so I'm calling your system: "<<s->getName()<<endl;
		}
		else  {
			modelName=pModel->Attribute("id");
			s=new System(modelName);
			if(verbose) cout<<"\tCreating system: "<<s->getName()<<endl;
		}
		
		
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
		if(!initParameters(pListOfParameters, parameter, verbose))
		{
			cout<<"\n\nI failed at parsing your Parameters.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}
		
		if(!verbose) cout<<"-";
		else if(pListOfFunctions) cout<<"\n\tReading list of Global Functions..."<<endl;
		if(pListOfFunctions)
		{
			if(!initGlobalFunctions(pListOfFunctions, s, parameter, verbose)) {
				cout<<"\n\nI failed at parsing your Global Functions.  Check standard error for a report."<<endl;
				delete s;
				return NULL;
			}
		}
		
		
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of MoleculeTypes..."<<endl;
		map<string,int> allowedStates;
		if(!initMoleculeTypes(pListOfMoleculeTypes, s, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your MoleculeTypes.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}
		
		
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Species..."<<endl;
		if(!initStartSpecies(pListOfSpecies, s, parameter, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your species.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}

		
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Reaction Rules..."<<endl;
		
		if(!initReactionRules(pListOfReactionRules, s, parameter, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your reaction rules.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}
		
		
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Observables..."<<endl;
		if(!initObservables(pListOfObservables, s, parameter, allowedStates, verbose))
		{
			cout<<"\n\nI failed at parsing your observables.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}
		
		
		/////////////////////////////////////////
		// Parse is finally over!  Now we just have to take care of some final details.
	
		//Finish up the output message
		if(!verbose) cout<<"-]\n";
		
		//Prepare the simulation
		s->prepareForSimulation();
		
		
		//Output the long form
		if(verbose) {
			cout<<"\n\nparse appears to be succussful.  Here, check your system:\n";
			s->printAllMoleculeTypes();
			s->printAllReactions();
			cout<<"-------------------------\n";
		}
		
		return s;
	}
	else
	{
		cout<<"\nError reading the file.  I could not find / open it, or it is not valid xml."<<endl;
	}
	
	
	return 0;
}




bool NFinput::initParameters(TiXmlElement *pListOfParameters, map <string,double> &parameter, bool verbose)
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
			if(verbose) cout<<"\t\t Identified parameter:\t"<<paramName<<"\tValue:"<<d<<endl;
		}
		return true;
	} catch (...) {
		cerr<<"Undefined exception thrown while parsing the parameters."<<endl;
		return false;
	}
}



bool NFinput::initGlobalFunctions(
	TiXmlElement * pListOfFunctions, 
	System * system,
	map <string,double> &parameter,
	bool verbose)
{
	try {
		vector <string> argNames;
		vector <string> argTypes;
		vector <string> paramNames;
		vector <double> paramValues;
		
		//Loop through the Function tags...
		TiXmlElement *pFunction;
		for ( pFunction = pListOfFunctions->FirstChildElement("Function"); pFunction != 0; pFunction = pFunction->NextSiblingElement("Function")) 
		{
			//Check if MoleculeType tag has a name...
			if(!pFunction->Attribute("id")) {
				cerr<<"!!!Error:  Function tag must contain the id attribute.  Quitting."<<endl;
				return false;	
			}
				
			//Read in the Function Name
			string funcName = pFunction->Attribute("id");
			if(verbose) cout<<"\t\tReading and Creating Function: "+funcName+"(";
			
			
			//Get the list of arguments for this function
			TiXmlElement *pListOfArgs = pFunction->FirstChildElement("ListOfArgs");
			if(pListOfArgs) {
				//Loop through each arguement
				TiXmlElement *pArg; bool firstArg = true;
				for ( pArg = pListOfArgs->FirstChildElement("Arg"); pArg != 0; pArg = pArg->NextSiblingElement("Arg")) 
				{
					if(verbose && !firstArg) cout<<", ";
					firstArg=false;
					//Check again for errors by making sure the argument has a name
					string argName, argType;
					if(!pArg->Attribute("id") || !pArg->Attribute("name") || !pArg->Attribute("type")) {
						if(verbose) cout<<" ?? ...\n";
						cerr<<"!!!Error:  Arg tag in Function: '" + funcName + "' must contain the all attributes including id, name, and type.  Quitting."<<endl;
						return false;
					}
					
					argName = pArg->Attribute("name");
					argType = pArg->Attribute("type");
					argNames.push_back(argName);
					argTypes.push_back(argType);
					
					if(verbose) cout<<argName;
				}
			}
			if(verbose) cout<<")"<<endl;
			
			
			//Get the list of Parameter Constants for this function
			TiXmlElement *pListOfParamConst = pFunction->FirstChildElement("ListOfParameterConstants");
			if(pListOfParamConst) {
				
				//Loop through each arguement
				TiXmlElement *pParamConst;
				for ( pParamConst = pListOfParamConst->FirstChildElement("ParameterConstant"); pParamConst != 0; pParamConst = pParamConst->NextSiblingElement("ParameterConstant")) 
				{
					//Check again for errors by making sure the parameter has a name
					if(!pParamConst->Attribute("id")) {
						cerr<<"!!!Error:  ParameterConstant tag in Function: '" + funcName + "' must contain a proper id.  Quitting."<<endl;
						return false;
					}
								
					string paramName = pParamConst->Attribute("id");
					if(parameter.find(paramName)==parameter.end()) {
						cerr<<"Could not find parameter constant: "<<paramName<<" when creating function "<<funcName<<". Quitting"<<endl;
						return false;
					}
					paramNames.push_back(paramName);
					paramValues.push_back(parameter.find(paramName)->second);
				}
			}
			
			
			//Read in the actual function definition
			TiXmlElement *pDefinition = pFunction->FirstChildElement("Definition");
			if(pDefinition) {
				if(!pDefinition->GetText()) {
					cerr<<"!!!Error:  Definition tag must actually contain a string for the function.  Quitting."<<endl;
										return false;
				}
				string functionDefintion = pDefinition->GetText();
				GlobalFunction *gf = new GlobalFunction(funcName, functionDefintion, argNames, argTypes, paramNames, paramValues);
				system->addGlobalFunction(gf);
				if(verbose) cout<<"\t\t\t= "<<functionDefintion<<endl;
			} else {
				cerr<<"!!!Error:  Definition tag for a function must exist!  Quitting."<<endl;
				return false;
			}
		
			argNames.clear();
			argTypes.clear();
			paramNames.clear();
			paramValues.clear();
		}
		
		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		return false;
	}
	
	
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
		//Organized as each vector in this vector contains the set of names of components that are identical

		
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
			if(typeName.compare("Null")==0 || typeName.compare("NULL")==0){
				if(verbose) cout<<"\t\tSkipping Moleculetype of name: '" + typeName + "'"<<endl;
				continue;
			}
			if(typeName=="Trash" || typeName=="trash" || typeName=="TRASH") {
					if(verbose) cout<<"\t\tSkipping Moleculetype of name: '" + typeName + "'"<<endl;
				continue;
			}
			if(verbose) cout<<"\t\tReading and Creating MoleculeType: "+typeName+"(";
	
			
			//Get the list of components in the moleculeType
			TiXmlElement *pListOfComp = pMoTypeEl->FirstChildElement("ListOfComponentTypes");
			if(pListOfComp)
			{
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
						if(compLabels.size()!=0) cout<<","+compName;
						else cout<<compName;
					}
					
					
					//First check if the component Name already exists, if so, we gotta do more!
					//This means that one of the sites are symmetric, so we must handle it correctly
					for(vector<string>::iterator it = compLabels.begin(); it != compLabels.end(); it++ ) {
						if((*it)==compName) {
							
							string newCompName = compName;
							string num="0"; bool matchedSiteName=false;
							for(unsigned int is=0; is<identicalComponents.size(); is++)
							{
								if(identicalComponents.at(is).at(0)==compName)
								{
									unsigned int lastIndex = identicalComponents.at(is).size();
									std::stringstream lastIndexStream; lastIndexStream << lastIndex+1;
									num = lastIndexStream.str();
									newCompName = compName+num;
									identicalComponents.at(is).push_back(newCompName);
									matchedSiteName = true;
								}
							}
							if(!matchedSiteName) {
								num="2";
								newCompName = compName+num;
								vector <string> v;
								v.push_back(compName);
								v.push_back(newCompName);
								identicalComponents.push_back(v);
							}
							compName = newCompName;
							if(verbose) cout<<num;
						}
					}
					compLabels.push_back(compName);
					
					
					//Take a look at the allowed states
					vector <string> possibleStateNames;
					TiXmlElement *pListOfAllowedStates = pComp->FirstChildElement("ListOfAllowedStates");
					if(pListOfAllowedStates)  {
						
						//Look at the list of allowed states
						TiXmlElement *pAlStates;
						int allowedStateCount = 0;
						for ( pAlStates = pListOfAllowedStates->FirstChildElement("AllowedState"); pAlStates != 0; pAlStates = pAlStates->NextSiblingElement("AllowedState")) 
						{
							if(!pAlStates->Attribute("id")) {
								cerr<<"!!!Error:  AllowedState tag in ComponentType '"+compName+"' of MoleculeType: '" + 
									typeName + "' must contain the id attribute.  Quitting."<<endl;
								return false;
							}
							string aState = pAlStates->Attribute("id");
							if(allowedStates.find(typeName+"_"+compName+"_"+aState)!=allowedStates.end()) continue;
							allowedStates[typeName+"_"+compName+"_"+aState] = allowedStateCount;
							allowedStateCount++;
							possibleStateNames.push_back(aState);
							if(verbose) cout<<"~"<<aState;
						}
					}
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
			
			
			//Create the moleculeType and register any symmetric sites we may have...
			MoleculeType *mt = new MoleculeType(typeName,compLabels,defaultCompState,possibleComponentStates,s);
			mt->addEquivalentComponents(identicalComponents);
			
			//Finally, clear the states and binding site labels that we read in
			compLabels.clear();
			defaultCompState.clear();
			possibleComponentStates.clear();
			identicalComponents.clear();
		}
				
		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		return false;
	}
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
				if(parameter.find(specCount)==parameter.end()) {
					cerr<<"Could not find parameter: "<<specCount<<" when creating species "<<speciesName<<". Quitting"<<endl;
					return false;
				}
				specCountInteger = (int)parameter.find(specCount)->second;
			}
			
			//Make sure we didn't try to create a negative number of molecules
			if(specCountInteger<0) {
				cerr<<"I cannot, in good conscience, make a negative number ("<<specCount<<") of species when creating species "<<speciesName<<". Quitting"<<endl;
				return false;
			}
			
			//If we're not going to make anything, well, then, don't make anything silly!  Stop here!
			if(specCountInteger==0) {
				if(verbose) cout<<"\t\tNot creating any instances of the Species: "<<speciesName<<" because you said I should make zero of them."<<endl;
				continue;
			}

			
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
			TiXmlElement *pMol;
			for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule")) 
			{
				//First get the type of molecule and retrieve the moleculeType object from the system
				string molName, molUid;
				if(!pMol->Attribute("name") || ! pMol->Attribute("id"))  {
					cerr<<"!!!Error.  Invalid 'Molecule' tag found when creating species '"<<speciesName<<"'. Quitting"<<endl;
					return false;
				} else {
					molName = pMol->Attribute("name");
					molUid = pMol->Attribute("id");
				}

				//Skip over null molecules
				if(molName=="Null" || molName=="NULL") {
					if(verbose) cout<<"\t\t\tSkipping Molecule of type: "<<molName<<" with local id: " << molUid<<endl;
					continue;
				}
				
				// Identify the moleculeType if we can (note that this call could potentially kill our code if we can't find the type);
				MoleculeType *mt = s->getMoleculeTypeByName(molName);
				if(verbose) cout<<"\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;

				
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
					cout<<"!!! warning: no list of components specified for molecule: '"<<molUid<<"' of species '"<<speciesName<<"'"<<endl;
				}

				//loop to create the actual molecules of this type
				vector <Molecule *> currentM;
				molecules.push_back(currentM);
				for(int m=0; m<specCountInteger; m++)
				{
					Molecule *m = mt->genDefaultMolecule();

					//Loop through the states and set the ones we need to set
					int k=0;
					for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
						m->setComponentState((*snIter).c_str(), (int)stateValue.at(k));
					molecules.at(molecules.size()-1).push_back(m);
				}
				
				//Reset the states for the next wave...
				stateName.clear();
				stateValue.clear();
			}

			
			///////////////////////////////////////////////////////////////
			//Here is where we add the bonds to the molecules in this species
			TiXmlElement *pListOfBonds = pListOfMol->NextSiblingElement("ListOfBonds");
			if(pListOfBonds)
			{
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
		return false;
	}
}




bool NFinput::FindReactionRuleSymmetry(
		TiXmlElement * pRxnRule, 
		System * s, 
		map <string,double> &parameter, 
		map<string,int> &allowedStates,
		map <string, component> &symComps,
		map <string, component> &symRxnCenter,
		bool verbose) 
{
	try {
		map <string, component> comps;
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////					
		//Grab the name of the rule
		string rxnName;
		if(!pRxnRule->Attribute("id")) {
			cerr<<"ReactionRule tag without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		} else {
			rxnName = pRxnRule->Attribute("id");
		}
		if(verbose) cout<<"\n\t\tReading Reaction Rule: "<<rxnName<<" to find symmetries...  ";
		
		
		TiXmlElement *pListOfReactantPatterns = pRxnRule->FirstChildElement("ListOfReactantPatterns");
		if(!pListOfReactantPatterns) {
			cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no reactant patterns!"<<endl;
			return true;
		}
		
		TiXmlElement *pReactant;
		for ( pReactant = pListOfReactantPatterns->FirstChildElement("ReactantPattern"); pReactant != 0; pReactant = pReactant->NextSiblingElement("ReactantPattern")) 
		{
			const char *reactantName = pReactant->Attribute("id");
			if(!reactantName) {
				cerr<<"\nReactant tag in reaction "<<rxnName<<" without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			}
			//if(verbose) cout<<"\t\t\tReading Reactant Pattern: "<<reactantName<<endl;
		
			//int NumOfSymComps = symComps.size();
			TiXmlElement *pListOfMols = pReactant->FirstChildElement("ListOfMolecules");
			if(pListOfMols) {
				if(!readPatternForSymmetry(pListOfMols, s, reactantName, comps, symComps, verbose)) return false;
			}
			else {
				cerr<<"\nReactant pattern "<<reactantName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
				return false;
			}
			
			//cout<<"("<<symComps.size() - NumOfSymComps<<")"<<endl;
		}
					
		
		//Read in the list of operations we need to perform in this rule
		TiXmlElement *pListOfOperations = pRxnRule->FirstChildElement("ListOfOperations");
		if(!pListOfOperations) {
			cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no operations!  This rule will do nothing!"<<endl;
			return true;
		}
		
		//First extract out the state changes
		TiXmlElement *pStateChange;
		for ( pStateChange = pListOfOperations->FirstChildElement("StateChange"); pStateChange != 0; pStateChange = pStateChange->NextSiblingElement("StateChange")) 
		{
			//Make sure all the information about the state change is here
			string site, finalState;
			if(!pStateChange->Attribute("site") || !pStateChange->Attribute("finalState")) {
				cerr<<"\nA specified state change operation in ReactionClass: '"+rxnName+"' does not "<<endl;
				cerr<<"have a valid site or finalState attribute.  Quitting."<<endl;
				return false;
			} else {
				site = pStateChange->Attribute("site");
				finalState = pStateChange->Attribute("finalState");
			}
			
			if(comps.find(site)!=comps.end()) {
				component c = comps.find(site)->second;
				MoleculeType *mt = c.mt;
				
				if(mt->isEquivalentComponent(c.name)) {
					symRxnCenter.insert(pair <string, component> (site,c));
					symComps.erase(site);
				}
			} else {
					cerr<<"\nError in ReactionClass: '"+rxnName+"'."<<endl;
					cerr<<"It seems that I couldn't find the states you are refering to."<<endl;
					return false;
			}
		
//			try {
//				component c = comps.find(site)->second;
//				int finalStateInt = allowedStates.find(c.t->getMoleculeTypeName()+"_"+c.name+"_"+finalState)->second;
//				ts->addStateChangeTransform(c.t,c.name,finalStateInt);
//			} catch (exception& e) {
//				cerr<<"Error in adding a state change operation in ReactionClass: '"+rxnName+"'."<<endl;
//				cerr<<"It seems that either I couldn't find the state, or the final state is not valid."<<endl;
//				return false;
//			}
		}
		
		
		//Search for symmetric sites in the bonds that are formed...
		TiXmlElement *pAddBond;
		for ( pAddBond = pListOfOperations->FirstChildElement("AddBond"); pAddBond != 0; pAddBond = pAddBond->NextSiblingElement("AddBond")) 
		{
			//Make sure all the information about the binding operation is here
			string site1, site2;
			if(!pAddBond->Attribute("site1") || !pAddBond->Attribute("site2")) {
				cerr<<"\nA specified binding operation in ReactionClass: '"+rxnName+"' does not "<<endl;
				cerr<<"have a valid site1 or site2 attribute.  Quitting."<<endl;
				return false;
			} else {
				site1 = pAddBond->Attribute("site1");
				site2 = pAddBond->Attribute("site2");
			}
			
			if(comps.find(site1)!=comps.end() && comps.find(site2)!=comps.end()) {
				component c1 = comps.find(site1)->second;
				component c2 = comps.find(site2)->second;
										
				MoleculeType *mt1 = c1.mt;
				MoleculeType *mt2 = c2.mt;
										
				if(mt1->isEquivalentComponent(c1.name)) {
					symRxnCenter.insert(pair <string, component> (site1,c1));
					symComps.erase(site1);
				}
				if(mt2->isEquivalentComponent(c2.name)) {
					symRxnCenter.insert(pair <string, component> (site2,c2));
					symComps.erase(site2);
				}
			} else {
				cerr<<"\nError in adding a binding operation in ReactionClass: '"+rxnName+"'."<<endl;
				cerr<<"It seems that either I couldn't find the binding sites you are refering to."<<endl;
				return false;
			}
		}
		
		//Next extract out removal of bonds
		TiXmlElement *pDeleteBond;
		for ( pDeleteBond = pListOfOperations->FirstChildElement("DeleteBond"); pDeleteBond != 0; pDeleteBond = pDeleteBond->NextSiblingElement("DeleteBond")) 
		{
			//Make sure all the information about the unbinding operation change is here
			string site1,site2;
			if(!pDeleteBond->Attribute("site1") || !pDeleteBond->Attribute("site2")) {
				cerr<<"\nA specified binding operation in ReactionClass: '"+rxnName+"' does not "<<endl;
				cerr<<"have a valid site1 or site2 attribute.  Quitting."<<endl;
				return false;
			} else {
				site1 = pDeleteBond->Attribute("site1");
				site2 = pDeleteBond->Attribute("site2");
			}
			
			if(comps.find(site1)!=comps.end() && comps.find(site2)!=comps.end()) {
				component c1 = comps.find(site1)->second;
				component c2 = comps.find(site2)->second;
										
				MoleculeType *mt1 = c1.mt;
				MoleculeType *mt2 = c2.mt;
										
				if(mt1->isEquivalentComponent(c1.name)) {
					symRxnCenter.insert(pair <string, component> (site1,c1));
					symComps.erase(site1);
				}
				if(mt2->isEquivalentComponent(c2.name)) {
					symRxnCenter.insert(pair <string, component> (site2,c2));
					symComps.erase(site2);
				}
				
			} else {
				cerr<<"\nError in adding an unbinding operation in ReactionClass: '"+rxnName+"'."<<endl;
				cerr<<"It seems that I couldn't find the binding sites you are refering to."<<endl;
				return false;
			}
		}
		
		
		if(verbose)
		if(symComps.size()>0 || symRxnCenter.size()>0) {
			cout<<"Found "<< symRxnCenter.size() <<" equivalent components in the rxn center and ";
			cout<<symComps.size()<<" outside rxn center."<<endl;
		} else {
			cout<<"No symmetries found.\n";
		}
		
		return true;
			
	} catch (...) {
		cout<<"caught something.."<<endl;
		return false;
	}
}





bool isValid(vector <vector <component> > &symRxnCenterComp, vector <int> &currentPos) {
	for(unsigned int i=0; i<symRxnCenterComp.size(); i++) {
		for(unsigned int j=i+1; j<symRxnCenterComp.size(); j++) {
			if(symRxnCenterComp.at(i).at(currentPos.at(i)).symPermutationName
					== symRxnCenterComp.at(j).at(currentPos.at(j)).symPermutationName) return false;
		}
	}
	return true;
}


void dumpState(vector <vector <component> > &symRxnCenterComp, vector <int> &currentPos) {
	cout<<"( ";
	for(unsigned int s=0; s<symRxnCenterComp.size(); s++)
		cout<<symRxnCenterComp.at(s).at(currentPos.at(s)).symPermutationName<<" ";
	cout<<")"<<endl;
}

void createSymMap(map<string,component> & symMap,
		vector <string> &uniqueId,
		vector <vector <component> > &symRxnCenterComp, 
		vector <int> &currentPos) 
{
	for(unsigned int s=0; s<symRxnCenterComp.size(); s++)
	{
		component c = symRxnCenterComp.at(s).at(currentPos.at(s));
		component newComp(c.mt, c.name);
		newComp.symPermutationName = c.symPermutationName;
		symMap.insert(pair <string, component> (uniqueId.at(s),newComp));
	}
}


bool NFinput::generateRxnPermutations(vector<map<string,component> > &permutations, 
		map<string,component> &symComps, 
		map<string,component> &symRxnCenter)
{
	//First, make sure we have some symmetric sites.  If not, just return and
	//carry on as normal...
	//if(symComps.size()==0 && symRxnCenter.size()==0) {
	if(symRxnCenter.size()==0) {
		map <string,component> m;
		permutations.push_back(m);
		return true;
	}
	cout<<"generating permutations..."<<endl;
	
	
	//Vectors to hold the info we need as we are generating things...
	vector <vector <component> > symRxnCenterComp;
	vector <string> uniqueId;
	vector <int> currentPos;
	
	
	//First, translate the map into vectors that have the 
	map<string, component>::iterator it;
	for ( it=symRxnCenter.begin() ; it != symRxnCenter.end(); it++)
	{
		//Create a new set to hold all the names of this equivalent site.
		vector <component> v;
		
		//Get the generalized component for this site or state
		component c = (*it).second;
		int *eq; int n_eq; //here we get the number of equivalent sites
		//if(c.type==component::BSITE) {
		c.mt->getEquivalencyClass(eq,n_eq, c.name);
		//} else if(c.type==component::STATE) {
		//	c.mt->getEquivalencyStateClass(eq,n_eq, c.name);
		//}
			
		
		//Loop through the equivalent sites or states and add it to the vector
		for(int e=0; e<n_eq; e++) {
			component newComp(c.mt, c.name);
			//if(c.type==component::BSITE) {
			string name(c.mt->getComponentName(eq[e]));
			newComp.symPermutationName=name;
			//}
			//else if(c.type==component::STATE) {
			//	string name(c.mt->getStateName(eq[e]));
			//	newComp.symPermutationName=name;
			//}
			v.push_back(newComp);
		}
		
		//finally put that list of equivalent sites on the main vector
		symRxnCenterComp.push_back(v);
		currentPos.push_back(0);
		uniqueId.push_back((*it).first);
	}
	
	// Output for debugging...
	for(unsigned int i=0; i<symRxnCenterComp.size(); i++) {
		cout<<"sym class "<<i<<": ";
		for(unsigned int j=0; j<symRxnCenterComp.at(i).size(); j++)
			cout<<symRxnCenterComp.at(i).at(j).symPermutationName<<" ";
		cout<<endl;
	}
	 
	
	int dumpCounter = 0;
	int activeIndex = 0; bool finished=false;
	if(isValid(symRxnCenterComp, currentPos)) {
		dumpCounter++;
		cout<<dumpCounter<<": ";
		dumpState(symRxnCenterComp, currentPos);
		
		map<string,component> symMap;
		createSymMap(symMap,uniqueId,symRxnCenterComp, currentPos);
		permutations.push_back(symMap);
	} else {
		cout<<"Invalid Configuration! : ";
		dumpState(symRxnCenterComp, currentPos);
	}
	while(true) {
		bool looped = false;
		while(((unsigned int)currentPos.at(activeIndex)+1)>=symRxnCenterComp.at(activeIndex).size()) {
			looped=true;
			currentPos.at(activeIndex)=0;
			activeIndex++;
			if((unsigned int)activeIndex>=symRxnCenterComp.size()) { finished=true; break; }
		}
		if(finished) break;
		currentPos.at(activeIndex) = currentPos.at(activeIndex)+1;
		if(looped) { activeIndex=0; } //Reset to the beginning of the list
		if(isValid(symRxnCenterComp, currentPos)) {
			dumpCounter++;
			cout<<dumpCounter<<": ";
			dumpState(symRxnCenterComp, currentPos);
			
			map<string,component> symMap;
			createSymMap(symMap,uniqueId,symRxnCenterComp, currentPos);
			permutations.push_back(symMap);
		} else {
			cout<<"Invalid Configuration! : ";
			dumpState(symRxnCenterComp, currentPos);
		}
	}
	
	
	return true;
}


bool lookup(component *&c, string id, map<string,component> &comps, map<string,component> &symMap) {
	try {
		if(symMap.find(id)!=symMap.end()) {
			component symC = symMap.find(id)->second;
			c = (&(comps.find(id)->second));
			c->symPermutationName=symC.symPermutationName;
		} else {
			if(comps.find(id)!=comps.end()) {
				c = (&(comps.find(id)->second));
				c->symPermutationName = c->name;
			} else {
				cerr<<"It seems that I couldn't find the binding sites or states you are refering to."<<endl;
				cerr<<"Could not find the component that matches the id: "<<id<<endl;
				return false;
			}
		}
	} catch (exception &e) {
		cerr<<"There was some problem when looking up the location of a particular component."<<endl;
		cerr<<"Could not find the component that matches the id: "<<id<<endl;
		return false;
	}
	return true;
}




bool NFinput::initReactionRules(
		TiXmlElement * pListOfReactionRules, 
		System * s, 
		map <string,double> &parameter, 
		map<string,int> &allowedStates, 
		bool verbose)
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
			generateRxnPermutations(permutations, symComps, symRxnCenter);
			
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
					
					if(permutations.size()>1) {
						stringstream out; out << (p+1);
						rxnName = rxnName + "_sym" + out.str();
					}
				}
				if(verbose) cout<<"\t\tCreating Reaction Rule: "<<rxnName<<endl;
				
				
				//First, read in the template molecules using these data structures
				map <string,TemplateMolecule *> reactants;
				map <string, component> comps;
				vector <TemplateMolecule *> templates;
			
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
						TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, allowedStates, reactantName, reactants, comps, symMap, verbose);
						if(tm==NULL) return false;
						templates.push_back(tm);
//////////////////////////////////////////tm->printDetails();
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
				// Create the TransformationSet so that we can collect all the operations that are specified for this rule
				TransformationSet *ts = new TransformationSet(templates);
				
			
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Read in the list of operations we need to perform in this rule
				TiXmlElement *pListOfOperations = pRxnRule->FirstChildElement("ListOfOperations");
				if(!pListOfOperations) {
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no operations!  This rule will do nothing!"<<endl;
					continue;
				}
			
				//First extract out the state changes
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
				
					if(verbose) cout<<"\t\t\t***Identified state change to site: "+site+" going to new state value: " + finalState<<endl;

					component *c;
					int finalStateInt = 0;
					if(!lookup(c, site, comps, symMap)) return false;
					try {
						finalStateInt = allowedStates.find(c->t->getMoleculeTypeName()+"_"+c->symPermutationName+"_"+finalState)->second;
					} catch (exception& e) {
						cerr<<"Error in adding a state change operation in ReactionClass: '"+rxnName+"'."<<endl;
						cerr<<"It seems that the final state is not valid."<<endl;
						return false;
					}
					ts->addStateChangeTransform(c->t,c->symPermutationName,finalStateInt);
				}
			
			
				//Next extract out the new bonds that are formed
				TiXmlElement *pAddBond;
				for ( pAddBond = pListOfOperations->FirstChildElement("AddBond"); pAddBond != 0; pAddBond = pAddBond->NextSiblingElement("AddBond")) 
				{
					//Make sure all the information about the binding operation is here
					string site1, site2;
					if(!pAddBond->Attribute("site1") || !pAddBond->Attribute("site2")) {
						cerr<<"A specified binding operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid site1 or site2 attribute.  Quitting."<<endl;
						return false;
					} else {
						site1 = pAddBond->Attribute("site1");
						site2 = pAddBond->Attribute("site2");
					}
					
					if(verbose) cout<<"\t\t\t***Identified binding of site: "+site1+" to binding of site " + site2<<endl;
					
					component *c1;
					component *c2;
					
					if(!lookup(c1, site1, comps, symMap)) return false;
					if(!lookup(c2, site2, comps, symMap)) return false;
					ts->addBindingTransform(c1->t, c1->symPermutationName, c2->t, c2->symPermutationName);
				}
				
				
				
				//Next extract out removal of bonds
				TiXmlElement *pDeleteBond;
				for ( pDeleteBond = pListOfOperations->FirstChildElement("DeleteBond"); pDeleteBond != 0; pDeleteBond = pDeleteBond->NextSiblingElement("DeleteBond")) 
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
					}
					
					if(verbose) cout<<"\t\t\t***Identified unbinding of site: "+site1+" to site " + site2<<endl;
					
					component *c1;
					component *c2;
					if(!lookup(c1, site1, comps, symMap)) return false;
					if(!lookup(c2, site2, comps, symMap)) return false;
						
					//Even though we had to make sure both ends exist, we really only need one transformation
					ts->addUnbindingTransform(c1->t, c1->symPermutationName);
				}
				
				
				//Next extract out anything that is destroyed
				TiXmlElement *pDelete;
				for ( pDelete = pListOfOperations->FirstChildElement("Delete"); pDelete != 0; pDelete = pDelete->NextSiblingElement("Delete")) 
				{
					//Make sure all the information about the state change is here
					string id;
					if(!pDelete->Attribute("id")) {
						cerr<<"A specified delete operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid id attribute.  Quitting."<<endl;
						return false;
					} else {
						try {
							id = pDelete->Attribute("id");
							if(verbose) cout<<"\t\t\t***Identified deletion of pattern: "+id<<"."<<endl;
							//cout<<"id: "<<id<<endl;
							component c = comps.find(id)->second;
							//cout<<"Templates.size() "<<templates.size()<<endl;
							//c.t->printDetails();
							ts->addDeleteMolecule(c.t);
						} catch (exception& e) {
							cerr<<"Error in adding an delete molecule operation in ReactionClass: '"+rxnName+"'."<<endl;
							cerr<<"It seems that I couldn't find the molecule to delete that you are refering to. (I was looking for ID: "<<pDelete->Attribute("id")<<endl;
							return false;
						}
					}
					
					
				}
				
				//Finally, figure out any new creations
				TiXmlElement *pAdd;
				for ( pAdd = pListOfOperations->FirstChildElement("Add"); pAdd != 0; pAdd = pAdd->NextSiblingElement("Add")) 
				{
					
					//Make sure all the information about the state change is here
					string id;
					if(!pAdd->Attribute("id")) {
						cerr<<"A specified add operation in ReactionClass: '"+rxnName+"' does not "<<endl;
						cerr<<"have a valid id attribute.  Quitting."<<endl;
						return false;
					} else {
						id = pAdd->Attribute("id");
						if(verbose) cout<<"\t\t\t***Identified addition of product pattern: "+id+"."<<endl;	
										
						SpeciesCreator *sc = NULL;
						
						//Go get the product pattern we need which will specify how to make this new species
						TiXmlElement *pListOfProductPatterns = pRxnRule->FirstChildElement("ListOfProductPatterns");
						if(!pListOfProductPatterns) {
							cerr<<"Error:: ReactionRule "<<rxnName<<" contains no product patterns, but needs at least one to add a species creation rule!"<<endl;
							return false;
						}
						TiXmlElement * pProduct;
						for ( pProduct = pListOfProductPatterns->FirstChildElement("ProductPattern"); pProduct != 0; pProduct = pProduct->NextSiblingElement("ProductPattern")) 
						{
							//First extract out the product Id
							string productName;
							if(!pProduct->Attribute("id")) {
								cerr<<"Product pattern in ReactionRule "+rxnName+" does not have a valid id attribute!"<<endl;
								return false;
							} else {
								productName = pProduct->Attribute("id");
							}
							
							//When we find the product pattern with the correct Id, then lets create it
							if(productName==id) {
								TiXmlElement *pListOfMols = pProduct->FirstChildElement("ListOfMolecules");
									if(pListOfMols) {
										
										vector <MoleculeType *> productMoleculeTypes;
										vector < vector <int> > stateInformation;
										vector < vector <int> > bindingSiteInformation;
										
										bool ok = NFinput::readProductPattern(pListOfMols,s,parameter,allowedStates, productName, 
												productMoleculeTypes, stateInformation, bindingSiteInformation, verbose);
										
										if(!ok) {
											cout<<"Could not read the list of product patterns for addition of a molecule in reaction "<<rxnName<<endl;
					
										} else if(productMoleculeTypes.size()>0){
											sc = new SpeciesCreator(productMoleculeTypes,stateInformation,bindingSiteInformation);
											ts->addAddMolecule(sc);
										}
									}
									else {
										cerr<<"Reactant product pattern "<<productName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
										return false;
									}
							}
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// With the transforations now set, Let's actually create the reaction (remember to finalize the TransformationSet!
				ts->finalize();
				ReactionClass *r;
				
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//  Read in the rate law for this reaction
				TiXmlElement *pRateLaw = pRxnRule->FirstChildElement("RateLaw");
				if(!pRateLaw){
					cerr<<"!!Error:: ReactionRule "<<rxnName<<" contains no rate law specification!"<<endl;
					return false;
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
						r = new BasicRxnClass(rxnName,0,ts);
						
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
							double rate=0; bool usedParam = false;
							try {
								rate = NFutil::convertToDouble(rateValue);
							} catch (std::runtime_error &e1) {
								if(parameter.find(rateValue)==parameter.end()) {
									cerr<<"Could not find parameter: "<<rateValue<<" when reading rate for rxn "<<rxnName <<". Quitting"<<endl;
									return false;
								}
								rate = parameter.find(rateValue)->second;
								usedParam = true;
							}
							if(verbose) {
								cout<<"\t\t\t\t...setting elementary rate to be: "<<rateValue;
								if(usedParam) cout<<" (has value: "<<rate<<")";
								cout<<endl;
							}
							
							r->setBaseRate(rate);
						}
						
						pRateConstant = pRateConstant->NextSiblingElement("RateConstant");
						if(pRateConstant!=NULL) {
							cout<<"\n\n!!Warning:: Multiple RateConstant tags present for RateLaw definition of "<<rxnName<<"."<<endl;
							cout<<"Only the first RateConstant given will be used..."<<endl;
						}	
					}
					else if(rateLawType=="Functional") {
						//Make sure that the rate constant exists
						TiXmlElement *pFunction = pRateLaw->FirstChildElement("Function");
						if(!pFunction) {
							cerr<<"Functional Rate Law definition for "<<rxnName<<" does not have Function specified!  Quiting"<<endl;
							return false;
						} else {
							//Get the rate constant value
							string functionName; 
							if(!pFunction->Attribute("name")) {
								cerr<<"Elementry Rate Law definition for "<<rxnName<<" does not have a valid function 'name'!  Quiting"<<endl;
								return false;
							} else {
								functionName = pFunction->Attribute("name");
								GlobalFunction *gf = s->getGlobalFunctionByName(functionName);
								if(gf==NULL) {
									cerr<<"When parsing reaction: '"<<rxnName<<"', could not identify function: '"<<functionName<<"' in\n";
									cerr<<"the system.  Therefore, I must abort."<<endl;
									exit(1);
								}
								
								//Create the Functional Reaction from the found function...
								r = new FunctionalRxnClass(rxnName,gf,ts);
							}
						}
					}
					
					
					
					////  To extend NFsim to parse more rate law types, add an extra else if clause here to catch the rate law
							
					else {
						cerr<<"!! I cannot yet interpret a Rate Law of 'type': "<<rateLawType<<".\n";
						cerr<<"  Remember, an elementry reaction needs to be specified as 'Ele' (case sensitive)."<<endl;
						cerr<<"Quitting."<<endl;
						return false;
					}
				} // end to the else statement that parses the rate law.
				
				//Finally, add the completed rxn rule to the system
				s->addReaction(r);
				comps.clear();
				
			} //end loop through all permutations
			
		} //end loop through all reaction rules
	
		//If we got here, then by golly, I think we have a new reaction rule
		return true;
		
	} catch (...) {
		cout<<"caught something when creating reaction rules.."<<endl;
		return false;
	}
}





bool NFinput::initObservables(
		TiXmlElement * pListOfObservables, 
		System * s, 
		map <string,double> &parameter, 
		map<string,int> &allowedStates, 
		bool verbose)
{
	try {
		//We will parse this in a similar manner to parsing species, except to say that we don't create
		//actual molecules, just template molecules.
		TiXmlElement *pObs;
		for ( pObs = pListOfObservables->FirstChildElement("Observable"); pObs != 0; pObs = pObs->NextSiblingElement("Observable")) 
		{
			//First get the observable name 			
			string observableId, observableName, observableType;
			if(!pObs->Attribute("id") || !pObs->Attribute("name") || !pObs->Attribute("type")) {
				cerr<<"Observable tag without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			} else {
				observableId = pObs->Attribute("id");
				observableName = pObs->Attribute("name");
				observableType = pObs->Attribute("type");
			}
			
			if(verbose) cout<<"\t\tCreating Observable: "<<observableName<<endl;
			
			//Now enter into the list of patterns
			TiXmlElement *pListOfPatterns = pObs->FirstChildElement("ListOfPatterns");
			if(!pListOfPatterns) {
				cout<<"\n\n!! Warning:: Observable "<<observableName<<" contains no patterns!"<<endl;
				continue;
			}
			
			//Now enter into the specific pattern making sure that it exists
			TiXmlElement *pPattern = pListOfPatterns->FirstChildElement("Pattern");
			if(!pPattern) {
				cout<<"\n\n!! Warning:: Observable "<<observableName<<" contains no patterns!"<<endl;
				continue;
			}
			
			//Make sure (for the time being) there is only one pattern
			if(pPattern->NextSiblingElement("Pattern")!=0) {
				cout<<"\n\n!!Warning:: Observable "<<observableName<<" contains multiple patterns!"<<endl;
				cout<<"                This is not yet supported, so only the first pattern will be read."<<endl;
			}
			
			//Go into the list of molecules that make up this pattern
			TiXmlElement *pListOfMol = pPattern->FirstChildElement("ListOfMolecules");
			if(!pListOfMol) {
				cout<<"\n\n!!Warning:: Observable "<<observableName<<" contains no molecules in its pattern!"<<endl;
				continue;
			}
			
			//Let the other function (readPattern) gather and create our template 
			map <string, TemplateMolecule *> templates;
			map <string, component> comps;
			map <string, component> symMap;
			TemplateMolecule *tempmol = readPattern(pListOfMol, s, parameter, allowedStates, observableName, templates, comps, symMap, verbose);
			if(tempmol==NULL) {
				cerr<<"Somehow, I couldn't parse out your pattern for observable "<<observableName<<" so I'll stop here."<<endl;
				return false;
			}
			
			
			//Finally, let's create the observable
			MoleculeType *moltype = tempmol->getMoleculeType();
			Observable *o  = new Observable(observableName.c_str(),tempmol);
			moltype->addObservable(o);
//////////////////////			tempmol->printDetails();
		}
		
		//Getting here means success!
		return true;
		
	} catch (...) {
		//oh no, what happened this time?
		return false;
	}
}
	



bool NFinput::readPatternForSymmetry(
		TiXmlElement * pListOfMol, 
		System * s,
		string patternName,
		map <string, component> &comps,
		map <string, component> &symComps,
		bool verbose)
{
	TiXmlElement *pMol;
	for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule")) 
	{
		//First get the type of molecule and retrieve the moleculeType object from the system
		string molName, molUid;
		if(!pMol->Attribute("name") || ! pMol->Attribute("id")) {
			cerr<<"!!!Error.  Invalid 'Molecule' tag found when creating pattern '"<<patternName<<"'. Quitting"<<endl;
			return false;
		} else {
			molName = pMol->Attribute("name");
			molUid = pMol->Attribute("id");
		}
			
		//Skip anything that is a null molecule
		if(molName=="Null" || molName=="NULL" || molName=="null") continue;
		if(molName=="Trash" || molName=="trash" || molName=="TRASH") continue;
			
		//Get the moleculeType and create the actual template
		MoleculeType *moltype = s->getMoleculeTypeByName(molName);
		
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
					return false;
				} else {
					compId = pComp->Attribute("id");
					compName = pComp->Attribute("name");
					compBondCount = pComp->Attribute("numberOfBonds");
				}
						
				
				//Declare and remember this component...
				component c(moltype, compName);
				comps.insert(pair <string, component> (compId,c));
				if(moltype->isEquivalentComponent(compName)) {
						symComps.insert(pair <string, component> (compId,c));
				} else {/*cout<<"no"<<endl;*/ }
				
				//Make sure the number of binding sites makes sense
				if(pComp->Attribute("numberOfBonds")) {
					string numOfBonds = pComp->Attribute("numberOfBonds");
					int numOfBondsInt = -1;
					try {
						numOfBondsInt = NFutil::convertToInt(numOfBonds);
					} catch (std::runtime_error e) {
						cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
						cerr<<e.what()<<endl;
						return false;
					}
				}
								
						
			} //end loop over components
		} //end if statement for compenents to exist
	
	}
	return true;
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
		bool verbose)
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
		
			//Skip anything that is a null molecule
			if(molName=="Null" || molName=="NULL" || molName=="null") continue;
			if(molName=="Trash" || molName=="trash" || molName=="TRASH") continue;
		
			//Get the moleculeType and create the actual template
			MoleculeType *moltype = s->getMoleculeTypeByName(molName);
			TemplateMolecule *tempmol = new TemplateMolecule(moltype);
			if(verbose) cout<<"\t\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;
		
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
												
					
					
					
					
					//Read in a state, if it is in fact has an associated state
					if(pComp->Attribute("state")) {
						string compStateValue = pComp->Attribute("state");
						
						//Make sure the given state is allowed
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
					
					//Check it as a binding site...
					if(pComp->Attribute("numberOfBonds")) {
						string numOfBonds = pComp->Attribute("numberOfBonds");
						int numOfBondsInt = -1;
						try {
							numOfBondsInt = NFutil::convertToInt(numOfBonds);
						} catch (std::runtime_error e) {
							cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
							cerr<<e.what()<<endl;
							return false;
						}
							
						//cout<<"bond value: "<< compId <<" " <<numOfBondsInt;
							
						//Look up this site in case we have some symmetry going on...
						component *symC;
						if(!lookup(symC, compId, comps, symMap)) {
							cerr<<"Could not find the symmetric component when creating a binding site, but there\n";
							cerr<<"should have been one!!  I don't know what to do, so I'll quit."<<endl;
							return false;
						}
							
						if(numOfBondsInt==0) {
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
				} //end loop over components
			} //end if statement for compenents to exist
			
			
			//Loop through the states and set the constraints we need to set
			int k=0;
			for(strVecIter = stateName.begin(); strVecIter != stateName.end(); k++, strVecIter++ )
			{
				tempmol->addStateValue((*strVecIter).c_str(),(int)stateValue.at(k));
			}
			for(strVecIter = emptyBondSite.begin(); strVecIter != emptyBondSite.end(); strVecIter++ )
			{
				tempmol->addEmptyBindingSite((*strVecIter).c_str());
			}
			
			//Update our data storage with the new template and empty out the things we don't need
			templates.insert(pair <string, TemplateMolecule *> (molUid,tempmol));
			tMolecules.push_back(tempmol);	
			stateName.clear();
			stateValue.clear();
			emptyBondSite.clear();
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
						TemplateMolecule::bind(tMolecules.at(bSiteMolIndex1),bSiteName1.c_str(),tMolecules.at(bSiteMolIndex2),bSiteName2.c_str());
						
						//Erase the bonds to make sure we don't add them again
						bSiteSiteMapping.erase(bSite1);
						bSiteMolMapping.erase(bSite1);
						bSiteSiteMapping.erase(bSite2);
						bSiteMolMapping.erase(bSite2);
					} else {
						cerr<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating pattern '"<<patternName<<"'. "<<endl;
						cerr<<"This may be caused because you are adding two bonds to one binding site or because you listed"<<endl;
						cerr<<"a binding site at the end of the pattern that does not exist.  Quitting"<<endl;
						return false;
					}
				} catch (exception& e) {
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
				tMolecules.at(bSiteMolIndex)->addOccupiedBindingSite(bSiteName.c_str());
			}
		} catch (exception& e) {
			cerr<<"Something went wacky when I was parsing pattern '"<<patternName<<"'."<<endl;
			cerr<<"It happened as I was trying to add an occupied binding site to a template molecule."<<endl;
			return false;
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
		
		return finalTemplate;
	} catch (...) {
		//Here's our final catch all!
		return NULL;
	}
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
}





bool NFinput::parseArguments(int argc, const char *argv[], map<string,string> &argMap)
{
	for(int a=1; a<argc; a++)
	{
		string s(argv[a]);
		
		//Find the strings that start with a minus, these are the flags
		if(s.compare(0,1,"-")==0)
		{
			string sFlag = s.substr(1,s.size()-1);
			if(sFlag.empty()) {
				cout<<"Warning: possible error in arguments.  You gave a '-' with a space"<<endl;
				cout<<"directly following. This is not a valid argument."<<endl<<endl;
				continue;
			}
			
			
			//See if the flag has some other input value that follows
			string sVal;
			if((a+1)<argc) {
				sVal=argv[a+1];
			}
			if(sVal.compare(0,1,"-")==0) {
				sVal = "";
			}
			
			//cout<<"found:  '"<<sFlag<<"' with arg: '"<<sVal<<"' "<<endl;
			if(argMap.find(sFlag)!=argMap.end()) {
				cout<<"Found two values for the same flag: '"<<sFlag<<"' so I am stopping"<<endl;
				return false;
			}
			
			argMap[sFlag] = sVal;
		}
	}
	return true;
}



int NFinput::parseAsInt(map<string,string> &argMap,string argName,int defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			int intVal = NFutil::convertToInt(strVal);
			return intVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as an integer."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}

double NFinput::parseAsDouble(map<string,string> &argMap,string argName,double defaultValue)
{
	if(argMap.find(argName)!=argMap.end()) {
		string strVal = argMap.find(argName)->second;
		try {
			double doubleVal = NFutil::convertToDouble(strVal);
			return doubleVal;
		} catch (std::runtime_error e) {
			cout<<endl<<"!!  Warning: I couldn't parse your flag '-"+argName+" "+strVal+"' as a double."<<endl;
			cout<<"!!  Using the default value of "<<defaultValue<<endl<<endl;;
		}
	}
	return defaultValue;
}




















