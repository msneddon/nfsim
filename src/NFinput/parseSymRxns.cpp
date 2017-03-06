#include "NFinput.hh"


using namespace NFinput;
using namespace std;


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

		//Grab the name of the rule
		string rxnName;
		if(!pRxnRule->Attribute("id")) {
			cerr<<"ReactionRule tag without a valid 'id' attribute.  Quiting"<<endl;
			return false;
		} else {
			rxnName = pRxnRule->Attribute("id");
		}
		if(verbose) cout<<"\n\t\tReading Reaction Rule: "<<rxnName<<" to find symmetries...  ";

		//Look at the list of patterns
		TiXmlElement *pListOfReactantPatterns = pRxnRule->FirstChildElement("ListOfReactantPatterns");
		if(!pListOfReactantPatterns) {
			cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!! Warning:: ReactionRule "<<rxnName<<" contains no reactant patterns!"<<endl;
			return true;
		}

		//Read the pattern for symmetry
		TiXmlElement *pReactant;
		for ( pReactant = pListOfReactantPatterns->FirstChildElement("ReactantPattern"); pReactant != 0; pReactant = pReactant->NextSiblingElement("ReactantPattern"))
		{
			const char *reactantName = pReactant->Attribute("id");
			if(!reactantName) {
				cerr<<"\nReactant tag in reaction "<<rxnName<<" without a valid 'id' attribute.  Quiting"<<endl;
				return false;
			}

			TiXmlElement *pListOfMols = pReactant->FirstChildElement("ListOfMolecules");
			if(pListOfMols) {
				if(!readPatternForSymmetry(pListOfMols, s, reactantName, comps, symComps, verbose)) return false;
			}
			else {
				cerr<<"\nReactant pattern "<<reactantName <<" in reaction "<<rxnName<<" without a valid 'ListOfMolecules'!  Quiting."<<endl;
				return false;
			}
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
					cerr<<"It seems that I couldn't find the states you are referring to."<<endl;
					cerr<<"Looking for site: "<<site<<endl;
					return false;
			}
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

				//Skip this if we are adding a bond in the product pattern....
				// @TODO:  FIX THIS!  should reject adds in molecule species that are newly added!
				//if(site1.find("RP")>=0 || site2.find("RP")>=0) continue;
			}

			// Handle site1 and site2 separately!
			//  If we can't find a site, look to see if it's on a product molecule.
			//  If so, then we'll move on and assume the missing site won't affect symmetry
			// TODO: think carefully w.r.t. symmetry and new Molecules.
			// --Justin
			if(comps.find(site1)!=comps.end() )
			{
				component c1 = comps.find(site1)->second;
				MoleculeType *mt1 = c1.mt;

				if( mt1->isEquivalentComponent(c1.name) )
				{
					symRxnCenter.insert(pair <string, component> (site1,c1));
					symComps.erase(site1);
				}
			}
		    else
		    {
		    	if ( site1.find("_PP") != string::npos )
		    	{
		    		if (verbose)
		    		{
			    		cout << "\n\t\t\tAddBond transform includes site '" << site1 << "' on a newly created molecule."
			    		     << "\n\t\t\t(I am ignoring this site with respect to symmetry.)"
			    		     <<	endl;
		    		}
		    	}
		    	else
		    	{
		    		cerr<<"\nError in adding a binding operation in ReactionClass: '"+rxnName+"'."<<endl;
					cerr<<"It seems that either I couldn't find the binding sites you are refering to."<<endl;
					return false;
		    	}
			}

			// Now handle site2.
			if( comps.find(site2)!=comps.end())
			{
				component c2 = comps.find(site2)->second;
				MoleculeType *mt2 = c2.mt;

				if( mt2->isEquivalentComponent(c2.name) )
				{
					symRxnCenter.insert(pair <string, component> (site2,c2));
					symComps.erase(site2);
				}
			}
		    else
		    {
		    	if ( site2.find("_PP") != string::npos )
		    	{
		    		if (verbose)
		    		{
			    		cout << "\n\t\t\tAddBond transform includes site '" << site1 << "' on a newly created molecule."
			    		     << "\n\t\t\t(I am ignoring this site with respect to symmetry.)"
			    		     <<	endl;
		    		}
		    	}
		    	else
		    	{
		    		cerr<<"\nError in adding a binding operation in ReactionClass: '"+rxnName+"'."<<endl;
					cerr<<"It seems that either I couldn't find the binding sites you are refering to."<<endl;
					return false;
		    	}
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

				//Skip this if we are messing with a bond in the product pattern....
				// @TODO:  FIX THIS!  should reject adds in molecule species that are newly added!
				//if(site1.find("RP")>=0 || site2.find("RP")>=0) continue;

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
				cout.flush();
				cerr<<"\nError in adding an unbinding operation in ReactionClass: '"+rxnName+"'."<<endl;
				cerr<<"It seems that I couldn't find the binding sites you are refering to."<<endl;
				cerr<<"Looking for site: "<<site1<<endl;
				cerr<<"Or site: "<<site2<<endl;
				return false;
			}
		}


		if(verbose) {
			if(symComps.size()>0 || symRxnCenter.size()>0) {
				cout<<"\n\t\t\tFound "<< symRxnCenter.size() <<" at rxn center, ";
				cout<<symComps.size()<<" outside rxn center."<<endl;

				cout<<"\t\t\t\tat the center: "<<endl;
				map <string,component>::iterator mapIter;
				for(mapIter=symRxnCenter.begin();mapIter!=symRxnCenter.end(); mapIter++) {
					cout<<"\t\t\t\t\t"<<mapIter->first<<"   "<<mapIter->second.name<<endl;

				}
				cout<<"\t\t\t\ton the side: "<<endl;
				for(mapIter=symComps.begin();mapIter!=symComps.end(); mapIter++) {
					cout<<"\t\t\t\t\t"<<mapIter->first<<"   "<<mapIter->second.name<<endl;
				}
			} else {
				cout<<"\t\t\tNo symmetry found.\n";
			}
		}

		return true;

	} catch (...) {
		cout<<"caught something.."<<endl;
		return false;
	}
}


void createMoleculeSymMap(
		map<string,component> &symMap,
		int mId,
		vector <vector <vector <component> > > &symmetries,
		vector <bool> &isRxnCenter,
		vector <vector <int> > &originalPosition,
		vector <int> &currentPosition
		)
{
	for(unsigned int k=0; k<currentPosition.size(); k++) {
		component *c;
		if(isRxnCenter.at(k)) {
			c = &symmetries.at(mId).at(0).at(originalPosition.at(k).at(currentPosition.at(k)));
		} else {
			c = &symmetries.at(mId).at(1).at(originalPosition.at(k).at(currentPosition.at(k)));
		}

		component compCopy(c->mt, c->name);
		compCopy.symPermutationName=c->symPermutationName;
		compCopy.uniqueId = c->uniqueId;
		compCopy.numOfBondsLabel=c->numOfBondsLabel;
		compCopy.stateConstraintLabel=c->stateConstraintLabel;

		symMap.insert(pair <string, component> (c->uniqueId,compCopy));
	}
}


void createFullSymMaps(
		vector<map<string,component> > &permutations, //The output will be stored here
		vector <vector <map <string,component> > > &symMaps, // the main input of the stored symMaps
		bool verbose
	)
{

	vector <int> currentPosition;
	for(unsigned int i=0; i<symMaps.size(); i++) currentPosition.push_back(0);

	int counter = 1;
	bool isDone = false;
	while(!isDone)
	{
		map <string, component> singlePermutation;
		for(unsigned int k=0; k<currentPosition.size(); k++) {
			map<string, component>::iterator it;
			if(verbose) {
				if(k==0) cout<<"\t\t\t\t"<<counter<<": [ ";
				else cout<<" [ ";
			}

			cout.flush();
			for ( it=symMaps.at(k).at(currentPosition.at(k)).begin() ; it != symMaps.at(k).at(currentPosition.at(k)).end(); it++)
			{
				component *c = &it->second;
				component compCopy(c->mt, c->name);
				compCopy.symPermutationName=c->symPermutationName;
				compCopy.uniqueId = c->uniqueId;
				compCopy.numOfBondsLabel=c->numOfBondsLabel;
				compCopy.stateConstraintLabel=c->stateConstraintLabel;
				singlePermutation.insert(pair <string, component> (c->uniqueId,compCopy));
				if(verbose) cout<<it->second.symPermutationName<<" ";
			}
			if(verbose) {
				if(k==currentPosition.size()-1) cout<<"]";
				else cout<<"] x";
			}
		}


		if(verbose) cout<<endl;
		counter++;
		permutations.push_back(singlePermutation);

		//Cycle to the next permutation, just as we did before
		int currentMap = symMaps.size()-1;
		do {
			currentPosition.at(currentMap)++;
			if(currentPosition.at(currentMap)>=(int)symMaps.at(currentMap).size()) {
				currentPosition.at(currentMap) = 0;
				currentMap--;
			} else { break; }
			if(currentMap<0) {
				isDone = true; break;
			}
		} while(true);


	}

	if(verbose) cout<<endl;
}


//Build an extended vector object called symmetries that
//saves all possible names for all possible components
void assembleFullSymmetryListOnRxnCenter(
	vector <vector <vector <component> > > &symmetries, //for output
	vector <string> &moleculeIds,    //also for output
	map<string,component> &symComps //the input of symmetric components
	)
{
	map<string, component>::iterator it;
	for ( it=symComps.begin() ; it != symComps.end(); it++)
	{
		//First, get the information about this symmetric component
		string id = it->first;
		component c = (*it).second;

		//Next, identify if this component is in a molecule we haven't yet considered...
		int length = id.find_last_not_of("_")-2;
		string thisMoleculeId = id.substr(0,length);
		int moleculeIndex = -1;

		for(unsigned int i=0; i<moleculeIds.size(); i++) {
			if(moleculeIds.at(i).compare(thisMoleculeId)==0) {
				moleculeIndex = i; break;
			}
		}
		if(moleculeIndex==-1) {
			moleculeIndex = moleculeIds.size();
			moleculeIds.push_back(thisMoleculeId);

			//Create the vector to store all of our potential permutations
			vector <vector <component> > v;
			vector <component> symRxnCenterComp;
			vector <component> symNonRxnCenterComp;
			v.push_back(symRxnCenterComp);
			v.push_back(symNonRxnCenterComp);
			symmetries.push_back(v);
		}

		//Get the list of equivalent components, and loop over them
		//to remember them in our new vector we are constructing
		int *eq; int n_eq; //here we get the number of equivalent sites
		c.mt->getEquivalencyClass(eq,n_eq, c.name);
		for(int e=0; e<n_eq; e++) {
			component newSymComp(c.mt, c.name);
			string name(c.mt->getComponentName(eq[e]));
			newSymComp.symPermutationName=name;
			newSymComp.uniqueId = id;
			symmetries.at(moleculeIndex).at(0).push_back(newSymComp);
		}
	}
}


//Checks only at reaction center (Assumes all symmetries are at reaction center!)
bool isMoleculePermuationValid(
	int mId,
	vector <vector <vector <component> > > &symmetries,
	vector <vector <int> > &originalPosition,
	vector <string> &uniqueComponents,
	vector <int> &currentPosition,
	vector <map <string,component> > &thisMoleculeSymMap
	)
{
	//First check if this permutation is self consistent (meaning each
	//unique component is mapped to only one symmetric component)
	vector <string> usedNames;
	for(unsigned int k=0; k<currentPosition.size(); k++) {
		component *c;
		c = &symmetries.at(mId).at(0).at(originalPosition.at(k).at(currentPosition.at(k)));

		for(unsigned int u=0; u<usedNames.size(); u++) {
			if(usedNames.at(u).compare(c->symPermutationName)==0) {
				return false;
			}
		}
		usedNames.push_back(c->symPermutationName);
	}

	return true;
}


bool NFinput::generateRxnPermutations(vector<map<string,component> > &permutations,
		map<string,component> &symComps,
		map<string,component> &symRxnCenter,
		bool verbose)
{
	//First, make sure we have some symmetric sites.  If not, just return and
	//carry on as normal...
	//if(symComps.size()==0 && symRxnCenter.size()==0) {
	if(symRxnCenter.size()==0) {
		map <string,component> m;
		permutations.push_back(m);
		return true;
	}

	if(verbose) cout<<"\t\t\tGenerating symmetric permutations..."<<endl;

	vector <vector <vector <component> > > symmetries;
	vector <string> moleculeIds;

	//Assemble the list of possible components for each symmetric class on a reaction center
	assembleFullSymmetryListOnRxnCenter(symmetries,moleculeIds,symRxnCenter);

	//Something to hold all of our intermediate symmetric site maps
	//which will be for each molecule, there is a vector of maps
	//that map component ids in the pattern to component names
	vector <vector <map <string,component> > > symMaps;

	// Now, from the full list, generate all permutations for each molecule
	for(unsigned int mId=0; mId<symmetries.size(); mId++)
	{
		//a vector containing the set of mappings for the molecule
		vector <map <string,component> > thisMoleculeSymMap;

		//extract out the necessary information
		vector <component> symRxnCenterComp = symmetries.at(mId).at(0);
		//vector <component> symNonRxnCenterComp = symmetries.at(mId).at(1);

		// Here we assemble the list of unique components.  For each permutation,
		// we want each unique component assigned only once.
		vector <bool> isRxnCenter;
		vector <vector <int> > originalPosition;
		vector <string> uniqueComponents;
		for(unsigned int k=0; k<symRxnCenterComp.size(); k++) {
			bool found = false;
			for(unsigned int i=0; i<uniqueComponents.size(); i++) {
				if(uniqueComponents.at(i).compare(symRxnCenterComp.at(k).uniqueId)==0) {
					originalPosition.at(i).push_back(k);
					found=true; break;
				}
			}
			if(!found) {
				uniqueComponents.push_back(symRxnCenterComp.at(k).uniqueId);
				isRxnCenter.push_back(true);
				vector <int> v; v.push_back(k);
				originalPosition.push_back(v);
			}
		}

		vector <int> currentPosition;
		for(unsigned int i=0; i<uniqueComponents.size(); i++) currentPosition.push_back(0);


		bool isDone = false;
		while(!isDone)
		{
			//Determine if the current permutation is valid...
			bool isValid = isMoleculePermuationValid(mId,
					symmetries,originalPosition,uniqueComponents,currentPosition,
					thisMoleculeSymMap);

			//If it is valid, then we have to add it to the list of permutations on this molecule
			if(isValid) {
				//Create the sym Map for this molecule
				map <string,component> moleculeSymMapForThisPermutation;
				createMoleculeSymMap(moleculeSymMapForThisPermutation,mId,
						symmetries,isRxnCenter,originalPosition,currentPosition);
				thisMoleculeSymMap.push_back(moleculeSymMapForThisPermutation);
			}


			// go to the next permutation, by ratcheting up like an odometer
			int currentComponent = uniqueComponents.size()-1;
			do {
				currentPosition.at(currentComponent)++;
				if(currentPosition.at(currentComponent)>=(int)originalPosition.at(currentComponent).size()) {
					currentPosition.at(currentComponent) = 0;
					currentComponent--;
				} else {
					break;
				}
				if(currentComponent<0) {
					isDone = true; break;
				}
			} while(true);

		}

		//Finally, add the symmetric maps of this molecule to the list...
		symMaps.push_back(thisMoleculeSymMap);
	}

	//and here we can create the symmetric permutation map
	createFullSymMaps(permutations, symMaps, verbose);

	return true;
}


bool NFinput::lookup(component *&c, string id, map<string,component> &comps, map<string,component> &symMap) {
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

		//Loop through the components of the molecule
		TiXmlElement *pListOfComp = pMol->FirstChildElement("ListOfComponents");
		if(pListOfComp)
		{
			TiXmlElement *pComp;
			for ( pComp = pListOfComp->FirstChildElement("Component"); pComp != 0; pComp = pComp->NextSiblingElement("Component"))
			{
				//Get the basic components of this molecule
				string compId, compName, compBondCount, compStateLabel;
				if(!pComp->Attribute("id") || !pComp->Attribute("name") || !pComp->Attribute("numberOfBonds")) {
					cerr<<"!!!Error.  Invalid 'Component' tag found when creating '"<<molUid<<"' of pattern '"<<patternName<<"'. Quitting"<<endl;
					return false;
				} else {
					compId = pComp->Attribute("id");
					compName = pComp->Attribute("name");
					compBondCount = pComp->Attribute("numberOfBonds");
					compStateLabel = "none";
					if(pComp->Attribute("state")) {
						compStateLabel = pComp->Attribute("state");
					}
				}


				//Declare and remember this component...
				component c(moltype, compName);
				c.numOfBondsLabel=compBondCount;
				c.stateConstraintLabel=compStateLabel;
				comps.insert(pair <string, component> (compId,c));

				if(moltype->isEquivalentComponent(compName)) {
						symComps.insert(pair <string, component> (compId,c));
				}

				//Make sure the number of binding sites makes sense here
				if(pComp->Attribute("numberOfBonds")) {
					string numOfBonds = pComp->Attribute("numberOfBonds");
					int numOfBondsInt = -1;
					if(numOfBonds!="+" && numOfBonds!="*" &&numOfBonds!="?") {
						try {
							numOfBondsInt = NFutil::convertToInt(numOfBonds);
						} catch (std::runtime_error e) {
							cerr<<"I couldn't parse the numberOfBonds value when creating pattern: "<<patternName<<endl;
							cerr<<e.what()<<endl;
							return false;
						}
					}
				}


			} //end loop over components
		} //end if statement for compenents to exist

	}
	return true;
}
