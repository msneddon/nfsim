
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
			//if(verbose) cout<<"\t\t\tReading Reactant Pattern: "<<reactantName<<endl;

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


//		//First extract out the state changes
//						TiXmlElement *pStateChange;
//						for ( pStateChange = pListOfOperations->FirstChildElement("StateChange"); pStateChange != 0; pStateChange = pStateChange->NextSiblingElement("StateChange"))
//						{
//							//Make sure all the information about the state change is here
//							string site, finalState;
//							if(!pStateChange->Attribute("site") || !pStateChange->Attribute("finalState")) {
//								cerr<<"A specified state change operation in ReactionClass: '"+rxnName+"' does not "<<endl;
//								cerr<<"have a valid site or finalState attribute.  Quitting."<<endl;
//								return false;
//							} else {
//								site = pStateChange->Attribute("site");
//								finalState = pStateChange->Attribute("finalState");
//							}
//
//							//First grab the component that is going to change...
//							component *c;
//							int finalStateInt = 0;
//							if(!lookup(c, site, comps, symMap)) return false;
//
//							//handle both increment and decrement states first...
//							if(finalState=="PLUS") {
//								if(!ts->addIncrementStateTransform(c->t,c->symPermutationName)) return false;
//
//								if(verbose) {
//									cout<<"\t\t\t***Identified increment state of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
//									cout<<") to new state value: oldStateValue+1"<<endl;
//								}
//							}
//							else if(finalState=="MINUS") {
//								if(!ts->addDecrementStateTransform(c->t,c->symPermutationName)) return false;
//
//
//								if(verbose) {
//									cout<<"\t\t\t***Identified decrement state of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
//									cout<<") to new state value: oldStateValue-1"<<endl;
//								}
//
//							}
//							else {
//
//								//Here, we handle your typical state change operation
//								try {
//									if(allowedStates.find(c->t->getMoleculeTypeName()+"_"+c->symPermutationName+"_"+finalState)==allowedStates.end()) {
//										cout<<"Error! in NFinput, when looking up state: "<<c->t->getMoleculeTypeName()+"_"+c->symPermutationName+"_"+finalState<<endl;
//										cout<<"Could not find this in the list of allowed states!  exiting!"<<endl;
//										exit(1);
//									}
//									finalStateInt = allowedStates.find(c->t->getMoleculeTypeName()+"_"+c->symPermutationName+"_"+finalState)->second;
//									//cout<<"found:"<<finalStateInt<<endl;
//								} catch (exception& e) {
//									cerr<<"Error in adding a state change operation in ReactionClass: '"+rxnName+"'."<<endl;
//									cerr<<"It seems that the final state is not valid."<<endl;
//									return false;
//								}
//								if(!ts->addStateChangeTransform(c->t,c->symPermutationName,finalStateInt)) return false;
//								if(verbose) {
//									cout<<"\t\t\t***Identified state change of site: "+c->t->getMoleculeTypeName()+"("+c->symPermutationName;
//									cout<<") to new state value: " + finalState<<endl;
//								}
//							}
//						}


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
				//if(site.find("RP")>=0) continue;
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



//////  DEPRECATED
//makes sure the component hasn't been used yet
//bool isValid(vector <vector <component> > &symRxnCenterComp, vector <int> &currentPos) {
//	for(unsigned int i=0; i<symRxnCenterComp.size(); i++) {
//		for(unsigned int j=i+1; j<symRxnCenterComp.size(); j++) {
//			if(symRxnCenterComp.at(i).at(currentPos.at(i)).symPermutationName
//					== symRxnCenterComp.at(j).at(currentPos.at(j)).symPermutationName) return false;
//		}
//	}
//	return true;
//}
//void dumpState(vector <vector <component> > &symRxnCenterComp, vector <int> &currentPos) {
//	cout<<"( ";
//	for(unsigned int s=0; s<symRxnCenterComp.size(); s++)
//		cout<<symRxnCenterComp.at(s).at(currentPos.at(s)).symPermutationName<<" ";
//	cout<<")"<<endl;
//}

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
		newComp.numOfBondsLabel=c.numOfBondsLabel;
		newComp.stateConstraintLabel=c.stateConstraintLabel;
		symMap.insert(pair <string, component> (uniqueId.at(s),newComp));
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
void assembleFullSymmetryList(
		vector <vector <vector <component> > > &symmetries, //for output
		vector <string> &moleculeIds,    //also for output
		map<string,component> &symComps, //the input of symmetric components
		bool isRxnCenter  //set to true if you are looking at reaction centers
		)
{
	//cout<<"\n\nreaction centers: "<<endl;
		map<string, component>::iterator it;
		for ( it=symComps.begin() ; it != symComps.end(); it++)
		{
			//cout<<it->first<<"   "<<it->second.name<<endl;

			//First, get the information about this symmetric component
			string id = it->first;
			component c = (*it).second;

			//Next, identify if this component is in a molecule we haven't yet considered...
			int length = id.find_last_not_of("_")-2;
			string thisMoleculeId = id.substr(0,length);
			int moleculeIndex = -1;

			for(unsigned int i=0; i<moleculeIds.size(); i++) {
				if(moleculeIds.at(i).compare(thisMoleculeId)==0)
					moleculeIndex = i; break;
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
				if(isRxnCenter) symmetries.at(moleculeIndex).at(0).push_back(newSymComp);
				else  symmetries.at(moleculeIndex).at(1).push_back(newSymComp);
			}
		}
}



//Build an extended vector object called symmetries that
//saves all possible names for all possible components
void assembleFullSymmetryListOnRxnCenter(
		vector <vector <vector <component> > > &symmetries, //for output
		vector <string> &moleculeIds,    //also for output
		map<string,component> &symComps //the input of symmetric components
		)
{
	//cout<<"\n\nreaction centers: "<<endl;
		map<string, component>::iterator it;
		for ( it=symComps.begin() ; it != symComps.end(); it++)
		{
			//cout<<it->first<<"   "<<it->second.name<<endl;

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



	///// Output for debugging
	//cout<<"currentPosition (originalPosition) Array"<<endl;
	//for(unsigned int k=0; k<currentPosition.size(); k++) {
	//	cout<<currentPosition.at(k)<<" (";
	//	cout<<originalPosition.at(k).at(currentPosition.at(k))<<") ";
	//}
	//cout<<endl;
}



//
void assembleOffRxnCenterSymClasses(
		vector <vector <vector <string> > > &offRxnCenterSymClasses, //the output
		vector <string> &moleculeIds, //input list of molecule names
		map<string,component> &symComps)  //input list of symmetric components off the rxn center
{
	offRxnCenterSymClasses.clear();
	for(unsigned int i=0; i<moleculeIds.size(); i++) {
		vector <vector <string> > v;
		offRxnCenterSymClasses.push_back(v);
	}

	//Loop through the symComponents
	map <string,component>::iterator mapIter;
	for(mapIter=symComps.begin();mapIter!=symComps.end(); mapIter++)
	{
		//cout<<mapIter->first<<"   "<<mapIter->second.name<<"   bond:"<<mapIter->second.numOfBondsLabel;
		//cout<<"  state:"<<mapIter->second.stateConstraintLabel<<endl;

		string id = mapIter->first;
		component c = mapIter->second;
		int length = id.find_last_not_of("_")-2;
		string thisMoleculeId = id.substr(0,length);

		int mIndex = -1;
		for(unsigned int i=0; i<moleculeIds.size(); i++) {
			if(thisMoleculeId.compare(moleculeIds.at(i))==0) {
				//cout<<" found at index: "<<i<<endl;
				mIndex=i;
				break;
			}
		}
		if(mIndex==-1) { cout<<"ERROR in parseSymRxns.cpp - in assebmly of off rxn center sym classes"<<endl; exit(1); }



		if(offRxnCenterSymClasses.at(mIndex).size()==0)
		{
			//Create a new class
			vector <string> newClass;
			newClass.push_back(id);
			offRxnCenterSymClasses.at(mIndex).push_back(newClass);
		} else {

			//Identify the sym class, or make your own
			bool hasMatchedExistingSymClass = false;

			//Loop through the existing sym classes
			for(unsigned int k=0; k<offRxnCenterSymClasses.at(mIndex).size(); k++) {

				//Look at the first element in the potential sym class (which must always exist)
				string potentialSymCompId = offRxnCenterSymClasses.at(mIndex).at(k).at(0);
				component potentialSymComp = symComps.find(potentialSymCompId)->second;

				//If the name matches, the state constraint label matches, and
				//the binding state (either not bound, or bound) matches
				if(c.name.compare(potentialSymComp.name)==0) {
					if(c.stateConstraintLabel.compare(potentialSymComp.stateConstraintLabel)==0) {

						int numOfBondsInt = -1;
						try {
							numOfBondsInt = NFutil::convertToInt(c.numOfBondsLabel);
						} catch (std::runtime_error e) {
							//This means the number of bonds label was not a number, but
							//make sure that it still matches
							if(c.numOfBondsLabel.compare(potentialSymComp.numOfBondsLabel)==0) {
								offRxnCenterSymClasses.at(mIndex).at(k).push_back(id);
								hasMatchedExistingSymClass = true;
							}
						}

						//Only add it if the number of bonds is zero.  If it is one, then
						//we are both bonded, but we don't know what we are bonded to, so
						//for now, we will create both permutations
						if(numOfBondsInt==0) {
							if(c.numOfBondsLabel.compare(potentialSymComp.numOfBondsLabel)==0) {
								offRxnCenterSymClasses.at(mIndex).at(k).push_back(id);
								hasMatchedExistingSymClass = true;
							}
						}
					}
				}
			}
			if(!hasMatchedExistingSymClass) {
				//then we make our own...
				vector <string> newClass;
				newClass.push_back(id);
				offRxnCenterSymClasses.at(mIndex).push_back(newClass);
			}

		}



	}


	///// Output for debugging
	//	cout<<endl<<endl;
	//	for(unsigned int k=0; k<offRxnCenterSymClasses.size(); k++) {
	//		cout<<"Molecule: "<<moleculeIds.at(k)<<endl;
	//		for(unsigned int j=0; j<offRxnCenterSymClasses.at(k).size(); j++) {
	//			cout<<"   ** Sym Class "<<j<<endl;
	//			for(unsigned int i=0; i<offRxnCenterSymClasses.at(k).at(j).size(); i++) {
	//				cout<<"      -Component: "<<offRxnCenterSymClasses.at(k).at(j).at(i)<<endl;
	//			}
	//		}
	//	}
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

	/// Output for debugging
	//map <string,component>::iterator mapIter;
	//for(mapIter=symComps.begin();mapIter!=symComps.end(); mapIter++) {
	//	cout<<mapIter->first<<"   "<<mapIter->second.name<<"   bond:"<<mapIter->second.numOfBondsLabel;
	//	cout<<"  state:"<<mapIter->second.stateConstraintLabel<<endl;
	//}


	if(verbose) cout<<"\t\t\tGenerating symmetric permutations..."<<endl;

	vector <vector <vector <component> > > symmetries;
	vector <string> moleculeIds;

	//Assemble the list of possible components for each symmetric class on a reaction center
	assembleFullSymmetryListOnRxnCenter(symmetries,moleculeIds,symRxnCenter);

	//Some more output for debugging
	//for(unsigned int i=0; i<symmetries.size(); i++) {
	//	cout<<">>>> "<<moleculeIds.at(i)<<endl;
	//	vector <component> symRxnCenterComp = symmetries.at(i).at(0);
	//	cout<<"  ** Rxn Center Components: "<<endl;
	//	for(unsigned int k=0; k<symRxnCenterComp.size(); k++) {
	//		cout<<"      -(" <<symRxnCenterComp.at(k).uniqueId<<"-"<<symRxnCenterComp.at(k).name<<")  "<<symRxnCenterComp.at(k).symPermutationName<<endl;
	//	}
	//}

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
		///// Output for debugging
		//cout<<"\n\nUnique Components:"<<endl;
		//for(unsigned int i=0; i<uniqueComponents.size(); i++) {
		//	cout<<uniqueComponents.at(i)<<"  count: " << originalPosition.at(i).size()<<endl;
		//	for(unsigned int j=0; j<originalPosition.at(i).size(); j++) {
		//		cout<<" index location: "<<originalPosition.at(i).at(j)<<endl;
		//	}
		//} cout<<endl;


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

				/// Output for debugging
				//cout<<"   [ ";
				//for(unsigned int k=0; k<currentPosition.size(); k++) {
				//	component *c;
				//	if(isRxnCenter.at(k)) {
				//		c = &symmetries.at(mId).at(0).at(originalPosition.at(k).at(currentPosition.at(k)));
				//		cout<<"*";
				//	} else {
				//		c = &symmetries.at(mId).at(1).at(originalPosition.at(k).at(currentPosition.at(k)));
				//	}
				//	cout<<c->symPermutationName<<" ";
				//}
				//cout<<" ]"<<endl;
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
				} //else {/*cout<<"no"<<endl;*/ }  //just a check for debugging

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















///////////////// OLD VERSION
//Checks if a given permutation on molecule components is valid (meaning
//we don't use the same name twice in any given permutation
//bool isMoleculePermuationValid(
//		int mId,
//		vector <vector <vector <component> > > &symmetries,
//		vector <bool> &isRxnCenter,
//		vector <vector <int> > &originalPosition,
//		vector <string> &uniqueComponents,
//		vector <int> &currentPosition,
//		vector <vector <vector <string> > > &offRxnCenterSymClasses,
//		vector <map <string,component> > &thisMoleculeSymMap
//		)
//{
//	//First check if this permutation is self consistent (meaning each
//	//unique component is mapped to only one symmetric component)
//	vector <string> usedNames;
//	for(unsigned int k=0; k<currentPosition.size(); k++) {
//		component *c;
//		if(isRxnCenter.at(k)) {
//			c = &symmetries.at(mId).at(0).at(originalPosition.at(k).at(currentPosition.at(k)));
//		} else {
//			c = &symmetries.at(mId).at(1).at(originalPosition.at(k).at(currentPosition.at(k)));
//		}
//
//		for(unsigned int u=0; u<usedNames.size(); u++) {
//			if(usedNames.at(u).compare(c->symPermutationName)==0) {
//				return false;
//			}
//		}
//		usedNames.push_back(c->symPermutationName);
//	}
//
//	///// Output for debugging
//	//cout<<"currentPosition (originalPosition) Array"<<endl;
//	//for(unsigned int k=0; k<currentPosition.size(); k++) {
//	//	cout<<currentPosition.at(k)<<" (";
//	//	cout<<originalPosition.at(k).at(currentPosition.at(k))<<") ";
//	//}
//	//cout<<endl;
//
//
//
//	///////////////////////////////////////////////////////
//	//Second, make sure we haven't already used a rearrangement of this permutation
//	//before by looking at the past maping sets
//
//	//if there is nothing to check, we can stop now
//	if(offRxnCenterSymClasses.at(mId).size()==0) return true;
//
//
//	//Init the mappedValue array
//	vector <vector <string> > mappedValue;
//	for(unsigned int j=0; j<offRxnCenterSymClasses.at(mId).size(); j++) {
//		vector <string> v;
//		for(unsigned int i=0; i<offRxnCenterSymClasses.at(mId).at(j).size(); i++) {
//			v.push_back("_");
//		}
//		mappedValue.push_back(v);
//	}
//
//    //Loop through all the past sym maps that we have accepted
//	for(unsigned int k=0; k<thisMoleculeSymMap.size(); k++) {
//
//		//Compare this map to the states to catch any duplicate rearrangements
//
//		//put the values into the mappedValue array
//		for(unsigned int j=0; j<offRxnCenterSymClasses.at(mId).size(); j++) {
//			for(unsigned int i=0; i<offRxnCenterSymClasses.at(mId).at(j).size(); i++) {
//				string compSymName = thisMoleculeSymMap.at(k).find(offRxnCenterSymClasses.at(mId).at(j).at(i))->second.symPermutationName;
//				mappedValue.at(j).at(i) = compSymName;
//			}
//		}
//
//		///// Output for debugging
//		//cout<<"\n\n  Mapped states... "<<endl;
//		//for(unsigned int j=0; j<offRxnCenterSymClasses.at(mId).size(); j++) {
//		//	cout<<"   ** Sym Class "<<j<<endl;
//		//	for(unsigned int i=0; i<offRxnCenterSymClasses.at(mId).at(j).size(); i++) {
//		//		cout<<"      -Component: "<<offRxnCenterSymClasses.at(mId).at(j).at(i);
//		//		cout<<"  "<<mappedValue.at(j).at(i)<<endl;
//		//	}
//		//}
//
//		//Because we know that each state value is only used once, for each symmetry class,
//		//we just count the number of checks.  If we checked off the same number as the size
//		//of the sym component vector, then we matched them all and so we just have a permutation
//		//on the classes.  Count those as a matched Sym class.  If we matched them all, then
//		//the permutation is not valid...
//
//		int matchedSymClasses = 0;
//		for(unsigned int j=0; j<offRxnCenterSymClasses.at(mId).size(); j++) { //loops through each sym class
//
//			if(offRxnCenterSymClasses.at(mId).at(j).size()==1) continue;
//
//
//			int checks = 0;
//			//loop through each component in the sym class
//			for(unsigned int i=0; i<offRxnCenterSymClasses.at(mId).at(j).size(); i++) {
//
//				///// Output for debugging
//				//cout<<" LOOKING AT     -Component: "<<offRxnCenterSymClasses.at(mId).at(j).at(i);
//				//cout<<"  "<<mappedValue.at(j).at(i)<<endl;
//
//
//				//Look up what we are trying to map here, from the set of input vectors
//				string componentToBeMapped = "";
//				for(unsigned int t=0; t<currentPosition.size(); t++) {
//
//					//Only look if we are NOT at the reaction center
//					if(isRxnCenter.at(t)) { continue; }
//
//					component *c = &symmetries.at(mId).at(1).at(originalPosition.at(t).at(currentPosition.at(t)));
//
//					//cout<<"    now checking out "<<currentPosition.at(t)<<" ("<<originalPosition.at(t).at(currentPosition.at(t))<<")"<<endl;
//					//cout<<"     "<<c->uniqueId<<endl;
//
//					if(c->uniqueId.compare(offRxnCenterSymClasses.at(mId).at(j).at(i))==0) {
//						componentToBeMapped=c->symPermutationName;
//						//cout<<" in new map, found :"<<componentToBeMapped<<"  "<<c->uniqueId<<endl;
//						break;
//					}
//				}
//				if(componentToBeMapped.empty()) {
//					cout<<"In parseSymRxns.cpp::Nothing mapped to this position, when checking if permutation is valid!!!"<<endl;
//					exit(1);
//				}
//
//
//
//				//See if that was mapped to any other location in the mapped value vector
//				for(unsigned int p=0; p<mappedValue.at(j).size(); p++) {
//					if(componentToBeMapped.compare(mappedValue.at(j).at(p))==0) {
//						checks++;
//					}
//				}
//			}
//			if(checks==(int)offRxnCenterSymClasses.at(mId).at(j).size()) {
//				matchedSymClasses++;
//			} else {
//				break;
//			}
//		}
//
//		//cout<<"Number of matched Sym Classes: "<<matchedSymClasses<<endl;
//
//		//If the number of matched sym classes equals the number of actual
//		//sym classes, then all of our sym classes are a permutation of this
//		//mapping, so this configuration is not valid.
//		if(matchedSymClasses==(int)offRxnCenterSymClasses.at(mId).size()) {
//			return false;
//		}
//	}
//
//	//well, if we got here, everything looks valid...
//	return true;
//}



//	bool NFinput::generateRxnPermutations(vector<map<string,component> > &permutations,
//			map<string,component> &symComps,
//			map<string,component> &symRxnCenter,
//			bool verbose)
//	{
//		//First, make sure we have some symmetric sites.  If not, just return and
//		//carry on as normal...
//		if(symComps.size()==0 && symRxnCenter.size()==0) {
//		//if(symRxnCenter.size()==0) {
//			map <string,component> m;
//			permutations.push_back(m);
//			return true;
//		}
//
//	///// Output for debugging
//	//	map <string,component>::iterator mapIter;
//	//	for(mapIter=symComps.begin();mapIter!=symComps.end(); mapIter++) {
//	//		cout<<mapIter->first<<"   "<<mapIter->second.name<<"   bond:"<<mapIter->second.numOfBondsLabel;
//	//		cout<<"  state:"<<mapIter->second.stateConstraintLabel<<endl;
//	//	}
//	//	exit(1);
//
//	/*For each molecule template in the pattern, we need to have the list of
//	 * all the symmetric
//
//	Strategy: generate every possible permutation, remove the ones that are clearly not
//	valid.  Then, when we create the template molecules from those permutations later
//
//	*/
//
//	//Some vectors to store the full symmetry list and the molecule names
//	vector <vector <vector <component> > > symmetries;
//	vector <string> moleculeIds;
//
//	//Assemble the full list of all possible states that each component,
//	//either at the reaction center or not, can be in.
//	assembleFullSymmetryList(symmetries,moleculeIds,symRxnCenter,true);
//	assembleFullSymmetryList(symmetries,moleculeIds,symComps,false);
//
//	//Now assemble the sym classes that are off the reaction center
//	//to tell us if ordering matters in the permutations...
//	vector <vector <vector <string> > > offRxnCenterSymClasses;
//	assembleOffRxnCenterSymClasses(
//			offRxnCenterSymClasses,
//			moleculeIds,
//			symComps);
//
//
//	///// Output for debugging
//		for(unsigned int i=0; i<symmetries.size(); i++) {
//			cout<<">>>> "<<moleculeIds.at(i)<<endl;
//			vector <component> symRxnCenterComp = symmetries.at(i).at(0);
//			vector <component> symNonRxnCenterComp = symmetries.at(i).at(1);
//
//			cout<<"  ** Rxn Center Components: "<<endl;
//			for(unsigned int k=0; k<symRxnCenterComp.size(); k++) {
//				cout<<"      -(" <<symRxnCenterComp.at(k).uniqueId<<"-"<<symRxnCenterComp.at(k).name<<")  "<<symRxnCenterComp.at(k).symPermutationName<<endl;
//			}
//
//			cout<<"  ** Other Symmetric Components: "<<endl;
//			for(unsigned int k=0; k<symNonRxnCenterComp.size(); k++) {
//				cout<<"      -(" <<symNonRxnCenterComp.at(k).uniqueId<<"-"<<symNonRxnCenterComp.at(k).name<<")  "<<symNonRxnCenterComp.at(k).symPermutationName<<endl;
//			}
//		}
//
//
//	//Something to hold all of our intermediate symmetric site maps
//	//which will be for each molecule, there is a vector of maps
//	//that map component ids in the pattern to component names
//	vector <vector <map <string,component> > > symMaps;
//
//
//	// Now, from the full list, generate all permutations for each molecule
//	for(unsigned int mId=0; mId<symmetries.size(); mId++)
//	{
//		//a vector containing the set of mappings for the molecule
//		vector <map <string,component> > thisMoleculeSymMap;
//
//		//extract out the necessary information
//		vector <component> symRxnCenterComp = symmetries.at(mId).at(0);
//		vector <component> symNonRxnCenterComp = symmetries.at(mId).at(1);
//
//
//		// Here we assemble the list of unique components.  For each permutation,
//		// we want each unique component assigned only once.
//		vector <bool> isRxnCenter;
//		vector <vector <int> > originalPosition;
//		vector <string> uniqueComponents;
//		for(unsigned int k=0; k<symRxnCenterComp.size(); k++) {
//			bool found = false;
//			for(unsigned int i=0; i<uniqueComponents.size(); i++) {
//				if(uniqueComponents.at(i).compare(symRxnCenterComp.at(k).uniqueId)==0) {
//					originalPosition.at(i).push_back(k);
//					found=true; break;
//				}
//			}
//			if(!found) {
//				uniqueComponents.push_back(symRxnCenterComp.at(k).uniqueId);
//				isRxnCenter.push_back(true);
//				vector <int> v; v.push_back(k);
//				originalPosition.push_back(v);
//			}
//		}
//		for(unsigned int k=0; k<symNonRxnCenterComp.size(); k++) {
//			bool found = false;
//			for(unsigned int i=0; i<uniqueComponents.size(); i++) {
//				if(uniqueComponents.at(i).compare(symNonRxnCenterComp.at(k).uniqueId)==0) {
//					originalPosition.at(i).push_back(k);
//					found=true; break;
//				}
//			}
//			if(!found) {
//				uniqueComponents.push_back(symNonRxnCenterComp.at(k).uniqueId);
//				isRxnCenter.push_back(false);
//				vector <int> v; v.push_back(k);
//				originalPosition.push_back(v);
//			}
//		}
//
//		///// Output for debugging
//		//		cout<<"\n\nUnique Components:"<<endl;
//		//		for(unsigned int i=0; i<uniqueComponents.size(); i++) {
//		//			cout<<uniqueComponents.at(i)<<"  count: " << originalPosition.at(i).size()<<endl;
//		//			for(unsigned int j=0; j<originalPosition.at(i).size(); j++) {
//		//				cout<<" index location: "<<originalPosition.at(i).at(j)<<endl;
//		//			}
//		//		} cout<<endl;
//
//
//		vector <int> currentPosition;
//		for(unsigned int i=0; i<uniqueComponents.size(); i++) currentPosition.push_back(0);
//
//
//		bool isDone = false;
//		while(!isDone)
//		{
//			//Determine if the current permutation is valid...
//			bool isValid = isMoleculePermuationValid(mId,
//					symmetries,isRxnCenter,originalPosition,uniqueComponents,currentPosition,
//					offRxnCenterSymClasses,
//					thisMoleculeSymMap);
//
//			//If it is valid, then we have to add it to the list of permutations on this molecule
//			if(isValid) {
//
//				//Create the sym Map for this molecule
//				map <string,component> moleculeSymMapForThisPermutation;
//				createMoleculeSymMap(moleculeSymMapForThisPermutation,mId,
//						symmetries,isRxnCenter,originalPosition,currentPosition);
//				thisMoleculeSymMap.push_back(moleculeSymMapForThisPermutation);
//
//				///// Output for debugging
////				cout<<"   [ ";
////				for(unsigned int k=0; k<currentPosition.size(); k++) {
////					component *c;
////					if(isRxnCenter.at(k)) {
////						c = &symmetries.at(mId).at(0).at(originalPosition.at(k).at(currentPosition.at(k)));
////						cout<<"*";
////					} else {
////						c = &symmetries.at(mId).at(1).at(originalPosition.at(k).at(currentPosition.at(k)));
////					}
////					cout<<c->symPermutationName<<" ";
////				}
////				cout<<" ]"<<endl;
//			}
//
//
//			// go to the next permutation, by ratcheting up like an odometer
//			int currentComponent = uniqueComponents.size()-1;
//			do {
//				currentPosition.at(currentComponent)++;
//				if(currentPosition.at(currentComponent)>=(int)originalPosition.at(currentComponent).size()) {
//					currentPosition.at(currentComponent) = 0;
//					currentComponent--;
//				} else {
//					break;
//				}
//				if(currentComponent<0) {
//					isDone = true; break;
//				}
//			} while(true);
//
//		}
//
//
//		//Finally, add the symmetric maps of this molecule to the list...
//		symMaps.push_back(thisMoleculeSymMap);
//	}
//
//
//	//With our sym maps for each molecule assembled, we can arrange them into permutations
//	createFullSymMaps(permutations, symMaps, verbose);
//	return true;

/////////////////////////// END OF NEW CODE

//
//	//Vectors to hold the info we need as we are generating things...
//	vector <vector <component> > symRxnCenterComp;
//	vector <string> uniqueId;
//	vector <int> currentPos;
//
//
//	//First, translate the map into vectors that have the list
//	//of possible state values
//	map<string, component>::iterator it;
//	for ( it=symRxnCenter.begin() ; it != symRxnCenter.end(); it++)
//	{
//		//Create a new set to hold all the names of this equivalent site.
//		vector <component> v;
//
//		//Get the generalized component for this site or state
//		component c = (*it).second;
//		int *eq; int n_eq; //here we get the number of equivalent sites
//		//if(c.type==component::BSITE) {
//		c.mt->getEquivalencyClass(eq,n_eq, c.name);
//		//} else if(c.type==component::STATE) {
//		//	c.mt->getEquivalencyStateClass(eq,n_eq, c.name);
//		//}
//
//
//		//Loop through the equivalent sites or states and add it to the vector
//		for(int e=0; e<n_eq; e++) {
//			component newComp(c.mt, c.name);
//			//if(c.type==component::BSITE) {
//			string name(c.mt->getComponentName(eq[e]));
//			newComp.symPermutationName=name;
//			//}
//			//else if(c.type==component::STATE) {
//			//	string name(c.mt->getStateName(eq[e]));
//			//	newComp.symPermutationName=name;
//			//}
//			v.push_back(newComp);
//		}
//
//		//finally put that list of equivalent sites on the main vector
//		symRxnCenterComp.push_back(v);
//		currentPos.push_back(0);
//		uniqueId.push_back((*it).first);
//	}
//
//	//Do the same for the nonRxnCenters ...
//	vector <vector <component> > symNotRxnCenterComp;
//	vector <string> uniqueIdNotRxnCenter;
//	vector <int> currentPosNotRxnCenter;
//	for ( it=symComps.begin() ; it != symComps.end(); it++)
//	{
//		//Create a new set to hold all the names of this equivalent site.
//		vector <component> v;
//
//		//Get the generalized component for this site or state
//		component c = (*it).second;
//		int *eq; int n_eq; //here we get the number of equivalent sites
//		c.mt->getEquivalencyClass(eq,n_eq, c.name);
//
//		//Loop through the equivalent sites or states and add it to the vector
//		for(int e=0; e<n_eq; e++) {
//			component newComp(c.mt, c.name);
//			string name(c.mt->getComponentName(eq[e]));
//			newComp.symPermutationName=name;
//			v.push_back(newComp);
//		}
//
//		//finally put that list of equivalent sites on the main vector
//		symNotRxnCenterComp.push_back(v);
//		currentPosNotRxnCenter.push_back(0);
//		uniqueIdNotRxnCenter.push_back((*it).first);
//	}
//
//	// Output for debugging...
//	for(unsigned int i=0; i<symRxnCenterComp.size(); i++) {
//		cout<<"sym class "<<i<<": ";
//		for(unsigned int j=0; j<symRxnCenterComp.at(i).size(); j++)
//			cout<<symRxnCenterComp.at(i).at(j).symPermutationName<<" ";
//		cout<<endl;
//	}
//
//	// Output for debugging...
//	for(unsigned int i=0; i<symNotRxnCenterComp.size(); i++) {
//		cout<<"sym class "<<i<<": ";
//		for(unsigned int j=0; j<symNotRxnCenterComp.at(i).size(); j++)
//			cout<<symNotRxnCenterComp.at(i).at(j).symPermutationName<<" ";
//		cout<<endl;
//	}
//
//
//	// This creates all valid permutations on the reaction center
//	int dumpCounter = 0;
//	int activeIndex = 0; bool finished=false;
//	if(isValid(symRxnCenterComp, currentPos)) {
//		dumpCounter++;
//		cout<<dumpCounter<<": ";
//		dumpState(symRxnCenterComp, currentPos);
//
//		map<string,component> symMap;
//		createSymMap(symMap,uniqueId,symRxnCenterComp, currentPos);
//		permutations.push_back(symMap);
//	} else {
//		cout<<"Invalid Configuration! : ";
//		dumpState(symRxnCenterComp, currentPos);
//	}
//	while(true) {
//		bool looped = false;
//		while(((unsigned int)currentPos.at(activeIndex)+1)>=symRxnCenterComp.at(activeIndex).size()) {
//			looped=true;
//			currentPos.at(activeIndex)=0;
//			activeIndex++;
//			if((unsigned int)activeIndex>=symRxnCenterComp.size()) { finished=true; break; }
//		}
//		if(finished) break;
//		currentPos.at(activeIndex) = currentPos.at(activeIndex)+1;
//		if(looped) { activeIndex=0; } //Reset to the beginning of the list
//		if(isValid(symRxnCenterComp, currentPos)) {
//			dumpCounter++;
//			cout<<dumpCounter<<": ";
//			dumpState(symRxnCenterComp, currentPos);
//
//			map<string,component> symMap;
//			createSymMap(symMap,uniqueId,symRxnCenterComp, currentPos);
//			permutations.push_back(symMap);
//		} else {
//			cout<<"Invalid Configuration! : ";
//			dumpState(symRxnCenterComp, currentPos);
//		}
//	}
//
//
//
//	// This will create all valid permutations on the non reaction centers
//	vector<map<string,component> > permutationsNotCenter;
//	dumpCounter = 0;
//	activeIndex = 0; finished=false;
//	if(isValid(symNotRxnCenterComp, currentPosNotRxnCenter)) {
//		dumpCounter++;
//		cout<<dumpCounter<<": ";
//		dumpState(symNotRxnCenterComp, currentPosNotRxnCenter);
//
//		map<string,component> symMapNotCenter;
//		createSymMap(symMapNotCenter,uniqueIdNotRxnCenter,symNotRxnCenterComp, currentPosNotRxnCenter);
//		permutationsNotCenter.push_back(symMapNotCenter);
//	} else {
//		cout<<"Invalid Configuration! : ";
//		dumpState(symNotRxnCenterComp, currentPosNotRxnCenter);
//	}
//	while(true) {
//		bool looped = false;
//		while(((unsigned int)currentPosNotRxnCenter.at(activeIndex)+1)>=symNotRxnCenterComp.at(activeIndex).size()) {
//			looped=true;
//			currentPosNotRxnCenter.at(activeIndex)=0;
//			activeIndex++;
//			if((unsigned int)activeIndex>=symNotRxnCenterComp.size()) { finished=true; break; }
//		}
//		if(finished) break;
//		currentPosNotRxnCenter.at(activeIndex) = currentPosNotRxnCenter.at(activeIndex)+1;
//		if(looped) { activeIndex=0; } //Reset to the beginning of the list
//		if(isValid(symNotRxnCenterComp, currentPosNotRxnCenter)) {
//			dumpCounter++;
//			cout<<dumpCounter<<": ";
//			dumpState(symNotRxnCenterComp, currentPosNotRxnCenter);
//
//			map<string,component> symMapNotCenter;
//			createSymMap(symMapNotCenter,uniqueIdNotRxnCenter,symNotRxnCenterComp, currentPos);
//			permutationsNotCenter.push_back(symMapNotCenter);
//		} else {
//			cout<<"Invalid Configuration! : ";
//			dumpState(symNotRxnCenterComp, currentPosNotRxnCenter);
//		}
//	}
//
//
//	exit(1);
//	return true;
//}
