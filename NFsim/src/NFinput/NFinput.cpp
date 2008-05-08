#include "NFinput.hh"









using namespace NFinput;
using namespace std;




System * NFinput::initializeFromXML(string filename)
{
	bool verbose = true;
	
	if(!verbose) cout<<"reading xml file ("+filename+")  [";
	if(verbose) cout<<"\tTrying to read xml model specification file: "<<filename<<endl;
	
	
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
		
		
		delete s;
		return NULL;
		
		if(!verbose) cout<<"-";
		else cout<<"\n\tReading list of Observables..."<<endl;
		if(!initObservables(pListOfObservables, s, parameter))
		{
			cout<<"\n\nI failed at parsing your observables.  Check standard error for a report."<<endl;
			delete s;
			return NULL;
		}
		
		/////////////////////////////////////////
	
		
		
		
		if(verbose) cout<<"\n\n-------------------------\n";
		for(int m=0; m<s->getNumOfMoleculeTypes(); m++)
			s->getMoleculeType(m)->printDetails();
		
		s->prepareForSimulation();
		
		s->registerOutputFileLocation("/home/msneddon/Desktop/new_xml/exampleOut.txt");
		s->outputAllObservableNames();
		s->outputAllObservableCounts();
		
		
		
		s->printAllReactions();	
		s->sim(20,200);
		s->printAllReactions();
		
		
		
		if(!verbose) cout<<"-]";
		else 
		
		return s;
	}
	else
	{
		cout<<"could not open.."<<endl;
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



/**
 * 
 * The strategy is to look at one MoleculeType at a time, make sure that moleculeType contains
 * valid information, then, and only then, create it.
 * 
 * 
 */
bool NFinput::initMoleculeTypes(TiXmlElement * pListOfMoleculeTypes, System * s,  map<string,int> &allowedStates, bool verbose) 
{
	try {
		vector <string> stateLabels;
		vector <string> bSiteLabels;
		
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
						if(stateLabels.size()!=0 || bSiteLabels.size()!=0) cout<<","+compName;
						else cout<<compName;
					}
					
					//Decide whether this component is a binding site or a state and act accordingly
					TiXmlElement *pListOfAllowedStates = pComp->FirstChildElement("ListOfAllowedStates");
					if(pListOfAllowedStates)  {
						//First check if the component Name already exists, if so, we can't handle that yet!
						for(vector<string>::iterator it = bSiteLabels.begin(); it != bSiteLabels.end(); it++ ) {
							if((*it)==compName) {
								cerr<<"!!!Error:  Binding Site Name: '"+compName+"' of MoleculeType: '"+typeName+"' was used more than once!\n";
								cerr<<"I don't know how to handle identical sites!! So I'm quitting. "<<endl;
								return false;
							}
						}
						stateLabels.push_back(compName);
						
						//Look at the list of allowed states, although I never check them in the future
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
							
							if(verbose) cout<<"~"<<aState;
						}
					}
					else {
						//First check if the component Name already exists, if so, we can't handle that yet!
						for(vector<string>::iterator it = bSiteLabels.begin(); it != bSiteLabels.end(); it++ ) {
							if((*it)==compName) {
								cerr<<"!!!Error:  Binding Site Name: '"+compName+"' of MoleculeType: '"+typeName+"' was used more than once!\n";
								cerr<<"I don't know how to handle identical sites!! So I'm quitting. "<<endl;
								return false;
							}
						}
						bSiteLabels.push_back(compName);
					}
				}
			}
			if(verbose) cout<<")"<<endl;
			
			
			//Now, actually create the molecule, starting with the binding sites
			unsigned int numOfBsites = bSiteLabels.size();
			string * bSiteNames = new string [numOfBsites];
			for(unsigned int b=0; b<numOfBsites; b++) {
				bSiteNames[b] = bSiteLabels.at(b);	
				if(false) cout<<"bSiteNames["<<b<<"] = "<<bSiteNames[b]<<endl;
			}
			
			//Here is where we create the states get setup
			unsigned int numOfStates = stateLabels.size();
			string * stateNames = new string [numOfStates];
			int * stateValues = new int [numOfStates];
			for(unsigned int i=0;i<numOfStates; i++)
			{
				stateNames[i] = stateLabels.at(i);
				stateValues[i] = 0;
				if(false) cout<<"stateNames["<<i<<"] = "<<stateNames[i]<<endl;
			}
			
			//With everything good to go, let's create the moleculeType
			new MoleculeType(typeName,stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);
	
			//Finally, clear the states and binding site labels that we read in
			bSiteLabels.clear();
			stateLabels.clear();
		}
				
		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		return false;
	}
}





bool NFinput::initStartSpecies(TiXmlElement * pListOfSpecies, System * s, map <string,double> &parameter, map<string,int> &allowedStates, bool verbose) 
{
	//map<string,int>::iterator iter;   
	//  for( iter = allowedStates.begin(); iter != allowedStates.end(); iter++ ) {
	//    cout << "state: " << iter->first << ", value: " << iter->second << endl;
	//  }
	
	
	
	
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
				
				// Identify the moleculeType if we can (this call could potentially kill our code if we can't find the type);
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
						else
						{
							//Add the b site mapping that will let us later easily connect binding sites
							//with the molecules involved
							bSiteSiteMapping[compId] = compName;
							bSiteMolMapping[compId] = molecules.size();
						}
					}

					//loop to create the actual molecules of this type
					vector <Molecule *> currentM;
					molecules.push_back(currentM);
					for(int m=0; m<specCountInteger; m++)
					{
						Molecule *m = new Molecule(mt);

						//Loop through the states and set the ones we need to set
						int k=0;
						for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
							m->setState((*snIter).c_str(), (int)stateValue.at(k));
						molecules.at(molecules.size()-1).push_back(m);
					}
					
					//Reset the states for the next wave...
					stateName.clear();
					stateValue.clear();
					
				}
				else
				{
					cout<<"!!! warning: no list of components specified for molecule: '"<<molUid<<"' of species '"<<speciesName<<"'"<<endl;
				}
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


bool NFinput::initReactionRules(TiXmlElement * pListOfReactionRules, System * s, map <string,double> &parameter, map<string,int> &allowedStates, bool verbose)
{
	try {
		
	//First, loop through all the rules
	TiXmlElement *pRxnRule;
	for ( pRxnRule = pListOfReactionRules->FirstChildElement("ReactionRule"); pRxnRule != 0; pRxnRule = pRxnRule->NextSiblingElement("ReactionRule")) 
	{
		//Grab the name of the rule
		string rxnName;
		if(!pRxnRule->Attribute("id")) {
			cerr<<"ReactionRule tag without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		} else {
			rxnName = pRxnRule->Attribute("id");
		}
		if(verbose) cout<<"\t\tCreating Reaction Rule: "<<rxnName<<endl;
		
		
		//First, read in the template molecules
		map <const char*,TemplateMolecule *, strCmp> reactants;
		vector <TemplateMolecule *> templates;
		
		
		//////////////////////////////////////////////////////////////////////
		//  Read in the Reactant Patterns for this rule
		TiXmlElement *pListOfReactantPatterns = pRxnRule->FirstChildElement("ListOfReactantPatterns");
		if(!pListOfReactantPatterns)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no reactant patterns!"<<endl;
			continue;
		}
		
		TiXmlElement *pReactant;
		for ( pReactant = pListOfReactantPatterns->FirstChildElement("ReactantPattern"); pReactant != 0; pReactant = pReactant->NextSiblingElement("ReactantPattern")) 
		{
			const char *reactantName = pReactant->Attribute("id");
			if(!reactantName)
			{
				cout<<"Reactant tag in reaction "<<rxnName<<" without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			}
			cout<<"\t\t\tReading Reactant Pattern: "<<reactantName<<endl;
			
			TiXmlElement *pListOfMols = pReactant->FirstChildElement("ListOfMolecules");
			if(pListOfMols)
			{
				TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, reactantName, reactants);
				if(tm==NULL) return false;
				templates.push_back(tm);
			}
			else
			{
				cout<<"Reactant pattern "<<reactantName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
				return false;
			}
			
			
					
			
			
			
		}
		
			map<const char*, TemplateMolecule *, strCmp>::iterator it;
				for ( it=reactants.begin() ; it != reactants.end(); it++ )
						cout << (*it).first << " => " << (*it).second->getMoleculeType()->getName() << endl;
							
				
		
		/////////////////////
		// Create the Reaction
		ReactionClass *r = new ReactionClass(rxnName.c_str(),templates,0);
		
		
		return false;
		
		
		//////////////////////////////////////////////////////////////////////
		//  Read in the Product Patterns for this rule
		TiXmlElement *pListOfProductPatterns = pRxnRule->FirstChildElement("ListOfProductPatterns");
		if(!pListOfProductPatterns)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no product patterns!"<<endl;
			continue;
		}
				
		
		TiXmlElement *pProduct;
		for ( pProduct = pListOfProductPatterns->FirstChildElement("ProductPattern"); pProduct != 0; pProduct = pProduct->NextSiblingElement("ProductPattern")) 
		{
			const char *productName = pProduct->Attribute("id");
			if(!productName)
			{
				cout<<"Product tag in reaction "<<rxnName<<" without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			}
			cout<<"\t\t\tReading Product Pattern: "<<productName<<endl;
			
			TiXmlElement *pListOfMols = pProduct->FirstChildElement("ListOfMolecules");
			
			if(pListOfMols)
			{
				//TemplateMolecule *tm = readPattern(pListOfMols, s, parameter, reactantName, reactants);
				//if(tm==NULL) return false;
				//templates.push_back(tm);
				
				
				if(!addTransformations(pListOfMols, 
						s, 
						parameter, 
						productName,
						reactants,
						r)) return false;
				
			}
			else
			{
				cout<<"Product pattern "<<productName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
				return false;
			}
			
		}
		
		
		
		
		//////////////////////////////////////////////////////////////////////
		//  Read in the rate law for this reaction
		TiXmlElement *pRateLaw = pRxnRule->FirstChildElement("RateLaw");
		if(!pRateLaw)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Error:: ReactionRule "<<rxnName<<" contains no rate law specification!"<<endl;
			return false;
		}
		const char *rateLawName = pRateLaw->Attribute("id");
		const char *rateLawType = pRateLaw->Attribute("type");
		if(!rateLawType)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Error:: ReactionRule "<<rxnName<<" rate law specification: cannot read 'type' attribute!"<<endl;
			return false;
		}
		else
		{
			cout<<"\t\t\tRate Law for Reaction is: "<<rateLawType<<endl;
			if(strcmp(rateLawType,"Ele")==0)
			{
				TiXmlElement *pListOfRateConstants = pRateLaw->FirstChildElement("ListOfRateConstants");
				if(!pListOfRateConstants)
				{
					cout<<"Elementry Rate Law definition for "<<rxnName<<" does not have listOfRateConstants specified!  Quiting"<<endl;
					return false;
				}
				
				TiXmlElement *pRateConstant = pListOfRateConstants->FirstChildElement("RateConstant");
				if(!pRateConstant)
				{
					cout<<"Elementry Rate Law definition for "<<rxnName<<" does not have RateConstants specified!  Quiting"<<endl;
					return false;
				}
				else
				{
					const char *value = pRateConstant->Attribute("value");
					if(!value)
					{
						cout<<"Elementry Rate Law definition for "<<rxnName<<" does not have a valid RateConstant value!  Quiting"<<endl;
						return false;
					}
					
					
					char *p; double rate = strtod(value,&p);
					if(*p != '\n' && *p != '\0')
					{
						try {
							rate = parameter.find(value)->second;
						} catch(exception& e) {
							cout<<"!!!!Invalid rate value for elementary reaction: '"<<rxnName<<". Quitting"<<endl;
							return false;
						}
					}
					cout<<"\t\t\t...setting elementary rate to be: "<<rate<<endl;
					r->setRate(rate);
				}
				
				pRateConstant = pRateConstant->NextSiblingElement("RateConstant");
				if(pRateConstant!=NULL)
				{
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Multiple RateConstant tags present for RateLaw definition of "<<rxnName<<"."<<endl;
					cout<<"The first RateConstant given will be used..."<<endl;
				}
				
			}
			////  To extend NFsim to parse more rate law types, add an extra if statement here to catch the rate law
			
			else
			{
				cout<<"!!!!!!!!!!!!!!!!!!!!!!!!  I cannot yet interpret a Rate Law of 'type': "<<rateLawType<<".\n";
				cout<<"  Remember, an elementry reaction needs to be specified as 'Ele' (case sensitive)."<<endl;
				cout<<"Quitting."<<endl;
				return false;
			}
		}
		
		
		//Finally, add the completed rxn rule to the system
		s->addReaction(r);
	}
	
	return true;
	
	
	} catch (...) {
		return false;
	}
}

bool NFinput::initObservables(TiXmlElement * pListOfObservables, System * s, map <string,double> &parameter)
{
	//We will parse this in a similar manner to parsing species, except to say that we don't create
	//actual molecules, just template molecules.
	TiXmlElement *pObs;
	for ( pObs = pListOfObservables->FirstChildElement("Observable"); pObs != 0; pObs = pObs->NextSiblingElement("Observable")) 
	{
		//First get the observable name 
		const char *observableName = pObs->Attribute("id");
		if(!observableName)
		{
			cout<<"Observable tag without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		}
		
		cout<<"\t\tCreating Observable: "<<observableName<<endl;
		
		//Now enter into the list of patterns
		TiXmlElement *pListOfPatterns = pObs->FirstChildElement("ListOfPatterns");
		if(!pListOfPatterns)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Observable "<<observableName<<" contains no patterns!"<<endl;
			continue;
		}
		
		//Now enter into the specific pattern making sure that it exists
		TiXmlElement *pPattern = pListOfPatterns->FirstChildElement("Pattern");
		if(!pPattern)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Observable "<<observableName<<" contains no patterns!"<<endl;
			continue;
		}
		if(pPattern->NextSiblingElement("Pattern")!=0)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Observable "<<observableName<<" contains multiple patterns!"<<endl;
			cout<<"                                   This is not yet supported, so only the first pattern will be read."<<endl;
		}
		
		//Go into the list of molecules that make up this pattern
		TiXmlElement *pListOfMol = pPattern->FirstChildElement("ListOfMolecules");
		if(!pListOfMol)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Observable "<<observableName<<" contains no molecules in its pattern!"<<endl;
			continue;
		}
		
		//Let the other function gather and create our template 
		map <const char*, TemplateMolecule *, strCmp> templates;
		TemplateMolecule *tm = readPattern(pListOfMol, s, parameter, observableName, templates);
		if(tm==NULL) return false;
		
		
		
		
		
		
		MoleculeType *mt = tm->getMoleculeType();
		Observable *o  = new Observable(observableName,tm);
		mt->addObservable(o);
	}
	
	return true;
}
	






TemplateMolecule *NFinput::readPattern(
		TiXmlElement * pListOfMol, 
		System * s, map <string,double> &parameter, 
		const char *patternName,
		map <const char*, TemplateMolecule *, strCmp> &templates)
{
	//A vector to hold each template molecule as we compose the entire template molecule
	vector < TemplateMolecule * > tMolecules;
	
	//Maps that map binding site ids into a molecule location in the molecules vector
	//and the name of the binding site
	map <const char*, const char *, strCmp> bSiteSiteMapping;
	map <const char*, int, strCmp> bSiteMolMapping;
	
	//vectors that keep track of the states and thier specified values as we create the templates
	vector <const char *> stateName;
	vector <double> stateValue;
	
	vector <const char *> emptyBondSite;
	
	//An iterator
	vector<const char *>::iterator snIter;
	
	//a character to use in string comparisons
	char *p;
	
	
	// Now loop through the molecules in the list
	TiXmlElement *pMol;
	for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule")) 
	{
		//First get the type of molecule and retrieve the moleculeType object from the system
		const char *molName = pMol->Attribute("name");
		const char *molUid = pMol->Attribute("id");
		if(!molName || ! molUid)
		{
			cout<<"!!!Error.  Invalid 'Molecule' tag found when creating pattern '"<<patternName<<"'. Quitting"<<endl;
			return NULL;
		}
		
		//Skip anything that is a null molecule
		if(strcmp(molName,"Null")==0) continue;
		
		string s2(molName);
		MoleculeType *mt = s->getMoleculeTypeByName(s2);
		cout<<"\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;
		
		//Loop through the components of the molecule in order to set state values
		TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
		if(pListOfComp)
		{
			TiXmlElement *pComp;
			for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component")) 
			{
				const char *compId = pComp->Attribute("id");
				const char *compName = pComp->Attribute("name");
				const char *compBondCount = pComp->Attribute("numberOfBonds");
				const char *compStateValue = pComp->Attribute("state");
				if(!compId || !compName || !compBondCount)
				{
					cout<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
					return 0;
				}
				if(compStateValue)
				{
					double s = strtod(compStateValue,&p);
					if(*p != '\n' && *p != '\0')
					{
						try {
							s = parameter.find(compStateValue)->second;
						} catch(exception& e) {
							cout<<"!!!!Invalid state value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
							return NULL;
						}
					}
					stateName.push_back(compName);
					stateValue.push_back(s);
				}
				else
				{
					const char* numOfBonds = pComp->Attribute("numberOfBonds");
					if(numOfBonds!=NULL)
					{
						double b = strtod(numOfBonds,&p);
						if(*p != '\n' && *p != '\0')
						{
							try {
								b = parameter.find(numOfBonds)->second;
							} catch(exception& e) {
								cout<<"!!!!Invalid numberOfBonds value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
								return NULL;
							}
						}
						if(b==0) emptyBondSite.push_back(compName);
						else
						{
							//Add the b site mapping that will let us later easily connect binding sites
							//with the molecules involved
							bSiteSiteMapping[compId] = compName;
							bSiteMolMapping[compId] = tMolecules.size();
						}
					}
				}
			}
		}
		
		//create the actual templateMolecule
		TemplateMolecule *m = new TemplateMolecule(mt);
				
		//Loop through the states and set the constraints we need to set
		int k=0;
		for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
			m->addStateValue((const char *)(*snIter),(int)stateValue.at(k));
		for(snIter = emptyBondSite.begin(); snIter != emptyBondSite.end(); snIter++ )
			m->addEmptyBindingSite((const char *)(*snIter));
			
		templates.insert(pair <const char *, TemplateMolecule *> (molUid,m));
		tMolecules.push_back(m);	
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
			const char *bondId = pBond->Attribute("id");
				const char *bSite1 = pBond->Attribute("site1");
				const char *bSite2 = pBond->Attribute("site2");
				if(!bondId || !bSite1 || !bSite2)
				{
					cout<<"!! Invalid Bond tag for pattern: "<<patternName<<".  Quitting."<<endl;
					return false;	
				}
				cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;
				
				
				//Get the information on this bond that tells us which molecules to connect
				try {
					const char *bSiteName1 = bSiteSiteMapping.find(bSite1)->second;
					int bSiteMolIndex1 = bSiteMolMapping.find(bSite1)->second;
					const char *bSiteName2 = bSiteSiteMapping.find(bSite2)->second;
					int bSiteMolIndex2 = bSiteMolMapping.find(bSite2)->second;
					TemplateMolecule::bind(tMolecules.at(bSiteMolIndex1),bSiteName1,tMolecules.at(bSiteMolIndex2),bSiteName2);
				} catch (exception& e) {
					cout<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
					return false;
				}
			}
		}
	
	TemplateMolecule *finalTemplate = tMolecules.at(0);
	
	tMolecules.clear();
	bSiteMolMapping.clear();
	bSiteSiteMapping.clear();
	return finalTemplate;
}





bool NFinput::addTransformations(TiXmlElement * pListOfProducts, 
		System * s, 
		map <string,double> &parameter, 
		const char *patternName,
		map <const char*, TemplateMolecule *, strCmp> &reactants,
		ReactionClass *r)
{
	
	cout<<"Adding Transformations..."<<endl;
	
	
	//Maps that map binding site ids into a molecule location in the molecules vector
	//and the name of the binding site
	map <const char*, const char *, strCmp> bSiteSiteMapping;
	map <const char*, const char *, strCmp> bSiteMolMapping;
	
	//Map that keeps track of unbinding reactions that were already added, so we don't
	//create two unbinding transformations when only one is desired
	//vector <TemplateMolecule *> alreadyUnboundTemplate;
	
	
	
	//Loop through the Molecules that make up this product
	TiXmlElement *pMol;
	for ( pMol = pListOfProducts->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule")) 
	{
		const char *molName = pMol->Attribute("name");
		const char *molUid = pMol->Attribute("id");
		if(!molName || ! molUid)
		{
			cout<<"!!!Error.  Invalid 'Molecule' tag found when analyzing product pattern '"<<patternName<<"'. Quitting."<<endl;
			return false;
		}
		
		
		cout<<"Reading molecule in product : " <<molName <<" has id: "<<molUid<<endl;
		
		
		TemplateMolecule *tm = reactants.find(molUid)->second;
		////////////// ADD ERROR CHECK HERE ///////////
		
		
		//Now loop through the particular components of this pattern to see what we have to change
		TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
		if(pListOfComp)
		{
			TiXmlElement *pComp;
			for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component")) 
			{
				const char *compId = pComp->Attribute("id");
				const char *compName = pComp->Attribute("name");
				const char *compBondCount = pComp->Attribute("numberOfBonds");
				const char *compStateValue = pComp->Attribute("state");
				if(!compId || !compName || !compBondCount)
				{
					cout<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
					return false;
				}
				
				//We have a state value we must check
				if(compStateValue)
				{
					char *p;
					double state_value = strtod(compStateValue,&p);
					if(*p != '\n' && *p != '\0')
					{
						try {
							state_value = parameter.find(compStateValue)->second;
						} catch(exception& e) {
							cout<<"!!!!Invalid state value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
							return false;
						}
					}
					
					//Make sure it is a transformation we have to make
					if(!tm->isStateValue(compName, (int)state_value))
					{
						cout<<"State Value is different, adding..."<<endl;
						Transformation::genStateChangeTransform(tm,compName,(int)state_value,r);
					} else { cout<<"State value is equivalent, so no transformation needed"<<endl; }
				}
				
				//Otherwise it is a binding site
				else
				{
					const char* numOfBonds = pComp->Attribute("numberOfBonds");
					if(numOfBonds!=NULL)
					{
						char *p; double b = strtod(numOfBonds,&p);
						if(*p != '\n' && *p != '\0')
						{
							try {
								b = parameter.find(numOfBonds)->second;
							} catch(exception& e) {
								cout<<"!!!!Invalid numberOfBonds value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
								return false;
							}
						}
						
						
						//////////////////////////////////
						///  NEEDS WORK HERE....
						//Check here for a possible unbinding event
						if((int)b==0)
						{
							if(tm->isBonded(compName))
							{
								//ch
								//int bSiteIndex = tm->getMoleculeType()->getBindingSiteIndex(compName);
								//int templateBSiteIndex = tm->getTemplateBsiteIndexFromMoleculeBsiteIndex(bSiteIndex);
								//TemplateMolecule * tm2 = tm->getBondedTemplateMolecule(templateBSiteIndex);
								
								//alreadyUnboundMapping.push_back()
								//Transformation::genUnbindingTransform(tm,compName, r);
								//cout<<"site "<<compName<<" appears to be bonded, 
							}
							else
							{
								cout<<"site "<<compName<<" appears to be bonded, so we don't have to do much here"<<endl;
							}
							
						}
						
						//Check here for a possible binding event
						else if((int)b==1)
						{
							if(!tm->isBonded(compName))
							{
								cout<<"site "<<compName<<" appears to be open, so mark it for probable addition of a bond."<<endl; 
								bSiteSiteMapping[compId] = compName;
								bSiteMolMapping[compId] = molUid;
							}
							else
							{
								cout<<"site "<<compName<<"is already bonded, so we don't have to do much here"<<endl;
							}
							
						}
						
					}
				}
				
				
				
			}
		}
		
		
//		map<const char*, const char*, strCmp>::iterator it;
	//	for ( it=bSiteSiteMapping.begin() ; it != bSiteSiteMapping.end(); it++ )
//			cout << (*it).first << " => " << (*it).second << endl;
	}
	
	//Here is where we add the bond generation mappings to the template molecules
	TiXmlElement *pListOfBonds = pListOfProducts->NextSiblingElement("ListOfBonds");
	if(pListOfBonds)
	{
		//First get the information on the bonds in the complex
		TiXmlElement *pBond;
		for ( pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond")) 
		{
			const char *bondId = pBond->Attribute("id");
			const char *bSite1 = pBond->Attribute("site1");
			const char *bSite2 = pBond->Attribute("site2");
			if(!bondId || !bSite1 || !bSite2)
			{
				cout<<"!! Invalid Bond tag for product: "<<patternName<<".  Quitting."<<endl;
				return false;	
			}
			cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;
							
							
			//Get the information on this bond that tells us which molecules to connect
			try {
				const char *bSiteName1 = bSiteSiteMapping.find(bSite1)->second;
				const char * bSiteMolId1 = bSiteMolMapping.find(bSite1)->second;
				TemplateMolecule *t1 = reactants.find(bSiteMolId1)->second;
				
				const char *bSiteName2 = bSiteSiteMapping.find(bSite2)->second;
				const char * bSiteMolId2 = bSiteMolMapping.find(bSite2)->second;
				TemplateMolecule *t2 = reactants.find(bSiteMolId2)->second;
					
				Transformation::genBindingTransform(t1,t2, bSiteName1, bSiteName2, r);
			} catch (exception& e) {
				cout<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating product '"<<patternName<<"'. Quitting"<<endl;
				return false;
			}
		}
	}
		
//			//First get the type of molecule and retrieve the moleculeType object from the system
//			const char *molName = pMol->Attribute("name");
//			const char *molUid = pMol->Attribute("id");
//			if(!molName || ! molUid)
//			{
//				cout<<"!!!Error.  Invalid 'Molecule' tag found when creating pattern '"<<patternName<<"'. Quitting"<<endl;
//				return NULL;
//			}
//			
//			//Skip anything that is a null molecule
//			if(strcmp(molName,"Null")==0) continue;
//			
//			MoleculeType *mt = s->getMoleculeTypeByName(molName);
//			cout<<"\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;
//			
//			//Loop through the components of the molecule in order to set state values
//			TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
//			if(pListOfComp)
//			{
//				TiXmlElement *pComp;
//				for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component")) 
//				{
//					const char *compId = pComp->Attribute("id");
//					const char *compName = pComp->Attribute("name");
//					const char *compBondCount = pComp->Attribute("numberOfBonds");
//					const char *compStateValue = pComp->Attribute("state");
//					if(!compId || !compName || !compBondCount)
//					{
//						cout<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
//						return 0;
//					}
//					if(compStateValue)
//					{
//						double s = strtod(compStateValue,&p);
//						if(*p != '\n' && *p != '\0')
//						{
//							try {
//								s = parameter.at(compStateValue);
//							} catch(exception& e) {
//								cout<<"!!!!Invalid state value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
//								return NULL;
//							}
//						}
//						stateName.push_back(compName);
//						stateValue.push_back(s);
//					}
//					else
//					{
//						const char* numOfBonds = pComp->Attribute("numberOfBonds");
//						if(numOfBonds!=NULL)
//						{
//							double b = strtod(numOfBonds,&p);
//							if(*p != '\n' && *p != '\0')
//							{
//								try {
//									b = parameter.at(numOfBonds);
//								} catch(exception& e) {
//									cout<<"!!!!Invalid numberOfBonds value for: '"<<molUid<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
//									return NULL;
//								}
//							}
//							if(b==0) emptyBondSite.push_back(compName);
//							else
//							{
//								//Add the b site mapping that will let us later easily connect binding sites
//								//with the molecules involved
//								bSiteSiteMapping[compId] = compName;
//								bSiteMolMapping[compId] = tMolecules.size();
//							}
//						}
//					}
//				}
//			}
//			
//			//create the actual templateMolecule
//			TemplateMolecule *m = new TemplateMolecule(mt);
//					
//			//Loop through the states and set the constraints we need to set
//			int k=0;
//			for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
//				m->addStateValue((const char *)(*snIter),(int)stateValue.at(k));
//			for(snIter = emptyBondSite.begin(); snIter != emptyBondSite.end(); snIter++ )
//				m->addEmptyBindingSite((const char *)(*snIter));
//				
//			templates.insert(pair <const char *, TemplateMolecule *> (molUid,m));
//			tMolecules.push_back(m);	
//			stateName.clear();
//			stateValue.clear();
//			emptyBondSite.clear();
//		}
//		
//		
//		//Here is where we add the bonds to the template molecules in the pattern
//		TiXmlElement *pListOfBonds = pListOfMol->NextSiblingElement("ListOfBonds");
//		if(pListOfBonds)
//		{
//			//First get the information on the bonds in the complex
//			TiXmlElement *pBond;
//			for ( pBond = pListOfBonds->FirstChildElement("Bond"); pBond != 0; pBond = pBond->NextSiblingElement("Bond")) 
//			{
//				const char *bondId = pBond->Attribute("id");
//					const char *bSite1 = pBond->Attribute("site1");
//					const char *bSite2 = pBond->Attribute("site2");
//					if(!bondId || !bSite1 || !bSite2)
//					{
//						cout<<"!! Invalid Bond tag for pattern: "<<patternName<<".  Quitting."<<endl;
//						return false;	
//					}
//					cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;
//					
//					
//					//Get the information on this bond that tells us which molecules to connect
//					try {
//						const char *bSiteName1 = bSiteSiteMapping.at(bSite1);
//						int bSiteMolIndex1 = bSiteMolMapping.at(bSite1);
//						const char *bSiteName2 = bSiteSiteMapping.at(bSite2);
//						int bSiteMolIndex2 = bSiteMolMapping.at(bSite2);
//						TemplateMolecule::bind(tMolecules.at(bSiteMolIndex1),bSiteName1,tMolecules.at(bSiteMolIndex2),bSiteName2);
//					} catch (exception& e) {
//						cout<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating pattern '"<<patternName<<"'. Quitting"<<endl;
//						return false;
//					}
//				}
//			}
	
	
	
	
	
	
	
	return true;
	
}















