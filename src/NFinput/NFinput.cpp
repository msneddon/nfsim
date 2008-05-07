#include "NFinput.hh"









using namespace NFinput;
using namespace std;




System * NFinput::initializeFromXML(char * filename)
{
	cout<<"\tTrying to read xml model specification file: "<<filename<<endl;
	TiXmlDocument doc(filename);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		cout<<"\t\tread was successful... beginning parse..."<<endl;
		
		//First declare our system
		System *s;
		
		//Read in the root node, which should give us the system's name
		TiXmlHandle hDoc(&doc);
		TiXmlElement *pModel = hDoc.FirstChildElement().Node()->FirstChildElement("model");
		if(!pModel) { cout<<"\tNo 'model' tag found.  Quitting."; return NULL; }
		
		//Make sure the basics are there
		const char * modelName =  pModel->Attribute("id");
		if(modelName)
		{
			s=new System((char *)modelName);
			cout<<"\tCreating system: "<<s->getName()<<endl;
		}
		else
		{
			s=new System("noname");
			cout<<"\tNo System name given, creating system: "<<s->getName()<<endl;
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
		//and save the parameters in an map we call parameter
		cout<<"\n\tReading parameter list..."<<endl;
		map<const char*, double, strCmp> parameter;
		if(!initParameters(pListOfParameters, parameter))
		{
			delete s;
			return NULL;
		}
		
		
		
		
		cout<<"\n\tReading list of MoleculeTypes..."<<endl;
		if(!initMoleculeTypes(pListOfMoleculeTypes, s))
		{
			delete s;
			return NULL;
		}
		
		cout<<"\n\tReading list of Species..."<<endl;
		if(!initStartSpecies(pListOfSpecies, s, parameter))
		{
			delete s;
			return NULL;
		}
		
		
//		cout<<"\n\tReading list of Reaction Rules..."<<endl;
//		if(!initReactionRules(pListOfReactionRules, s, parameter))
//		{
//			delete s;
//			return NULL;
//		}
//		
		cout<<"\n\tReading list of Observables..."<<endl;
		if(!initObservables(pListOfObservables, s, parameter))
		{
			delete s;
			return NULL;
		}
		
		/////////////////////////////////////////
	
		cout<<"\n\n-------------------------\n";
		for(int m=0; m<s->getNumOfMoleculeTypes(); m++)
			s->getMoleculeType(m)->printDetails();
		
		s->prepareForSimulation();
		
		s->registerOutputFileLocation("/home/msneddon/Desktop/new_xml/exampleOut.txt");
		s->outputAllObservableNames();
		s->outputAllObservableCounts();
		return s;
	}
	else
	{
		cout<<"could not open.."<<endl;
	}
	
	
	return 0;
}




bool NFinput::initParameters(TiXmlElement *pListOfParameters, map <const char*,double, strCmp> &parameter)
{
	/*TiXmlElement *pParamElement;
	for ( pParamElement = pListOfParameters->FirstChildElement("Parameter"); pParamElement != 0; pParamElement = pParamElement->NextSiblingElement("Parameter")) 
	{
		const char * a = pParamElement->Attribute("id");
		if(!a)
		{
			cout<<"\t\t!!A Parameter is undefined! It is missing the 'id' attribute!  Quitting.\n";
			return false;
		}
			
		const char * v = pParamElement->Attribute("value");
		if(!v)
		{
			cout<<"\t\t!!A Parameter '"<<a<<"' does not have the 'value' attribute defined! Quitting.\n";
			return false;
		}
			
		char *p;
		double d = strtod(v,&p);
		if(*p != '\n' && *p != '\0')
		{
			cout<<"\t\t!!Parameter '"<<a<<"' has a 'value' attribute that is not well formatted!\n\t\t";
			cout<<"You gave '"<<v<<"', but that is not a valid number! Quitting.\n";
			return false;	
		}
			
		parameter[a]=d;
		cout<<"\t\t Identified parameter:\t"<<a<<"\tValue:"<<d<<endl;
	}
	*/
	return true;
}


bool NFinput::initMoleculeTypes(TiXmlElement * pListOfMoleculeTypes, System * s) 
{
/*	bool output = false;
	vector <const char*> stateLabels;
	vector <const char*> bSiteLabels;
	
	
	TiXmlElement *pMoTypeEl;
	for ( pMoTypeEl = pListOfMoleculeTypes->FirstChildElement("MoleculeType"); pMoTypeEl != 0; pMoTypeEl = pMoTypeEl->NextSiblingElement("MoleculeType")) 
	{
		const char *typeName = pMoTypeEl->Attribute("id");
		if(!typeName)
		{
			cout<<"!!!Error:  MoleculeType tag must contain the id attribute.  Quitting."<<endl;
			return false;	
		}
		if(strcmp(typeName,"Null")==0)
			continue;
			
		
		cout<<"\t\tReading and Creating Moleculetype: "<<typeName<<"(";
		TiXmlElement *pListOfComp = pMoTypeEl->FirstChildElement("ListOfComponentTypes");
		if(pListOfComp)
		{
		
			
			//Grab the needed info on the molecule
			TiXmlElement *pComp;
			for ( pComp = pListOfComp->FirstChildElement("ComponentType"); pComp != 0; pComp = pComp->NextSiblingElement("ComponentType")) 
			{
				const char *comp = pComp->Attribute("id");
				if(stateLabels.size()!=0 || bSiteLabels.size()!=0)cout<<",";	
				
				TiXmlElement *pListOfAllowedStates = pComp->FirstChildElement("ListOfAllowedStates");
				if(pListOfAllowedStates)
				{
					stateLabels.push_back(comp);
					cout<<comp;
					TiXmlElement *pAlStates;
					for ( pAlStates = pListOfAllowedStates->FirstChildElement("AllowedState"); pAlStates != 0; pAlStates = pAlStates->NextSiblingElement("AllowedState")) 
					{
						const char *aState = pAlStates->Attribute("id");
						cout<<"~"<<aState;
					}
				}
				else
				{
					bSiteLabels.push_back(comp);
					cout<<comp;
				}
			}
		}
		cout<<")"<<endl;
		
		//Now, actually create the molecule
		int numOfBsites = bSiteLabels.size();
		char ** bSiteNames = new char * [numOfBsites];
		for(int b=0;b<numOfBsites; b++)
		{
			bSiteNames[b] = (char *)bSiteLabels.at(b);	
			if(output) cout<<"bSiteNames["<<b<<"] = "<<bSiteNames[b]<<endl;
		}
		
		
		int numOfStates = stateLabels.size();
		char ** stateNames = new char * [numOfStates];
		int * stateValues = new int [numOfStates];
		for(int i=0;i<numOfStates; i++)
		{
			
			stateNames[i] = (char *)stateLabels.at(i);
			stateValues[i] = 0;
			if(output) cout<<"stateNames["<<i<<"] = "<<stateNames[i]<<endl;
		}
		
		
		new MoleculeType(typeName,stateNames,stateValues,numOfStates,bSiteNames,numOfBsites,s);

		bSiteLabels.clear();
		stateLabels.clear();
	}
			
			
	*/
	return true;
}





bool NFinput::initStartSpecies(TiXmlElement * pListOfSpecies, System * s, map <const char*,double, strCmp> &parameter) 
{
	
/*	//A vector to hold molecules as we are creating the species
	vector < vector <Molecule *> > molecules;
	
	//A vector that maps binding site ids into a molecule location in the molecules vector
	//and the name of the binding site
	map <const char*, const char *, strCmp> bSiteSiteMapping;
	map <const char*, int, strCmp> bSiteMolMapping;
	
	vector <const char *> stateName;
	vector <double> stateValue;
	
	vector<const char *>::iterator snIter;
	
	
	TiXmlElement *pSpec;
	for ( pSpec = pListOfSpecies->FirstChildElement("Species"); pSpec != 0; pSpec = pSpec->NextSiblingElement("Species")) 
	{
		//First get the species name 
		const char *speciesName = pSpec->Attribute("id");
		if(!speciesName)
		{
			cout<<"Species tag without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		//Get the number of molecules of this species to create
		const char *specCount = pSpec->Attribute("concentration");
		if(!specCount)
		{
			cout<<"Species "<<speciesName<<" does not have a 'concentration' attribute.  Quitting"<<endl;
			return false;
		}
		
		//Try to parse out the number of species, or look it up in the parameter map
		char *p; double c = strtod(specCount,&p);
		if(*p != '\n' && *p != '\0')
		{
			try {
				c = parameter.at(specCount);
			} catch(exception& e) {
				cout<<"Could not find parameter: "<<specCount<<" when creating species "<<speciesName<<". Quitting"<<endl;
				return false;
			}
		}
		
		cout<<"\t\tCreating "<<c<<" instances of the Species: "<<speciesName<<endl;
		TiXmlElement *pListOfMol = pSpec->FirstChildElement("ListOfMolecules");
		if(!pListOfMol)
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!! Warning:: Species "<<speciesName<<" contains no molecules!"<<endl;
			continue;
		}
		
			
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Now loop through the molecules
		TiXmlElement *pMol;
		for ( pMol = pListOfMol->FirstChildElement("Molecule"); pMol != 0; pMol = pMol->NextSiblingElement("Molecule")) 
		{
			//First get the type of molecule and retrieve the moleculeType object from the systme
			const char *molName = pMol->Attribute("name");
			const char *molUid = pMol->Attribute("id");
			if(!molName || ! molUid)
			{
				cout<<"!!!Error.  Invalid 'Molecule' tag found when creating species '"<<speciesName<<"'. Quitting"<<endl;
				return false;
			}
			
			if(strcmp(molName,"Null")==0) continue;
			
			MoleculeType *mt = s->getMoleculeTypeByName(molName);
			//cout<<"\t\t\tIncluding Molecule of type: "<<molName<<" with local id: " << molUid<<endl;
			
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
						cout<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of species '"<<speciesName<<"'. Quitting"<<endl;
						return false;
					}
					if(compStateValue)
					{
						double s = strtod(compStateValue,&p);
						if(*p != '\n' && *p != '\0')
						{
							try {
								s = parameter.at(compStateValue);
							} catch(exception& e) {
								cout<<"!!!!Invalid state value for: '"<<molUid<<"' when creating species '"<<speciesName<<"'. Quitting"<<endl;
								return false;
							}
						}
						stateName.push_back(compName);
						stateValue.push_back(s);
					}
					else
					{
						//Add the b site mapping that will let us later easily connect binding sites
						//with the molecules involved
						bSiteSiteMapping[compId] = compName;
						bSiteMolMapping[compId] = molecules.size();
					}
				}
				
				//loop to create the actual molecules
				vector <Molecule *> currentM;
				molecules.push_back(currentM);
				for(int m=0; m<c; m++)
				{
					Molecule *m = new Molecule(mt);
					
					//Loop through the states and set the ones we need to set
					int k=0;
					for(snIter = stateName.begin(); snIter != stateName.end(); k++, snIter++ )
						m->setState((const char *)(*snIter), (int)stateValue.at(k));
					molecules.at(molecules.size()-1).push_back(m);
				}
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
				const char *bondId = pBond->Attribute("id");
				const char *bSite1 = pBond->Attribute("site1");
				const char *bSite2 = pBond->Attribute("site2");
				if(!bondId || !bSite1 || !bSite2)
				{
					cout<<"!! Invalid Bond tag for species: "<<speciesName<<".  Quitting."<<endl;
					return false;	
				}
				//cout<<"reading bond "<<bondId<<" which connects "<<bSite1<<" to " <<bSite2<<endl;
				
				
				//Get the information on this bond that tells us which molecules to connect
				try {
					const char *bSiteName1 = bSiteSiteMapping.at(bSite1);
					int bSiteMolIndex1 = bSiteMolMapping.at(bSite1);
					const char *bSiteName2 = bSiteSiteMapping.at(bSite2);
					int bSiteMolIndex2 = bSiteMolMapping.at(bSite2);
				
					for(int j=0;j<c;j++)
						Molecule::bind(molecules.at(bSiteMolIndex1).at(j),bSiteName1,molecules.at(bSiteMolIndex2).at(j),bSiteName2);
				} catch (exception& e) {
					cout<<"!!!!Invalid site value for bond: '"<<bondId<<"' when creating species '"<<speciesName<<"'. Quitting"<<endl;
					return false;
				}
			}
		}
		
		//Tidy up and clear the lists for the next species
		for(unsigned int i=0; i<molecules.size(); i++)
				molecules.at(i).clear();
		molecules.clear();
		bSiteMolMapping.clear();
		bSiteSiteMapping.clear();
	}
	
	*/
	//If we got here, then we are indeed successful
	return true;
}


bool NFinput::initReactionRules(TiXmlElement * pListOfReactionRules, System * system, map <const char*,double, strCmp> &parameter)
{
	
	
	return true;
}

bool NFinput::initObservables(TiXmlElement * pListOfObservables, System * s, map <const char*,double, strCmp> &parameter)
{
/*	//actual molecules, just template molecules.
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
		TemplateMolecule *tm = readPattern(pListOfMol, s, parameter, observableName);
		if(tm==NULL) return false;
		
		MoleculeType *mt = tm->getMoleculeType();
		Observable *o  = new Observable(observableName,tm);
		mt->addObservable(o);
	}
	
	*/
	return true;
}
	






TemplateMolecule *NFinput::readPattern(TiXmlElement * pListOfMol, System * s, map <const char*,double, strCmp> &parameter, const char *patternName)
{
	//A vector to hold each template molecule as we compose the entire template molecule
/*	vector < TemplateMolecule * > tMolecules;
	
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
		
		MoleculeType *mt = s->getMoleculeTypeByName(molName);
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
							s = parameter.at(compStateValue);
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
								b = parameter.at(numOfBonds);
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
					const char *bSiteName1 = bSiteSiteMapping.at(bSite1);
					int bSiteMolIndex1 = bSiteMolMapping.at(bSite1);
					const char *bSiteName2 = bSiteSiteMapping.at(bSite2);
					int bSiteMolIndex2 = bSiteMolMapping.at(bSite2);
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
	
	return finalTemplate;*/return 0;
}





















