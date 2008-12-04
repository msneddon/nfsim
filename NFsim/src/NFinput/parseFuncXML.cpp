#include "NFinput.hh"





using namespace NFinput;
using namespace std;



bool createLocalFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			System *s,
			map <string,double> &parameter,
			bool verbose);



bool createFunction(string name,
		string expression,
		vector <string> &argNames,
		vector <string> &refNames,
		vector <string> &refTypes,
		System *s,
		map <string,double> &parameter,
		bool verbose)
{



	cout<<endl<<endl;
	cout<<"creating the functions."<<endl;

	//bool referencesFunction = false;
	//bool hasLocal

	//Step 1: determine if it is local or global
	if(argNames.size()==0) { //No arguments, so must be a global function...
		cout<<"must be a global function..."<<endl;
		vector <string> varRefNames;
		vector <string> varRefTypes;
		vector <string> paramNames;


		int otherFuncRefCounter=0;
		for(unsigned int rn=0; rn<refNames.size(); rn++) {
			if(refTypes.at(rn)=="Function") {
				cout<<"identified function reference."<<endl;
				otherFuncRefCounter++;
			} else if(refTypes.at(rn)=="Constant") {
				paramNames.push_back(refNames.at(rn));
			} else {
				varRefNames.push_back(refNames.at(rn));
				varRefTypes.push_back(refTypes.at(rn));
			}
		}

		//Treat as usual global function
		if(otherFuncRefCounter==0) {

			//FuncFactory::create();
			GlobalFunction *gf = new GlobalFunction(name, expression,
					varRefNames, varRefTypes, paramNames, s);
			if(!s->addGlobalFunction(gf)) {
				cerr<<"!!!Error:  Function name '"<<name<<"' has already been used.  You can't have two\n";
				cerr<<"functions with the same name, so I'll just stop now."<<endl;
				return false;
			}
		}
		//Treat as special rate law global function
		else if (otherFuncRefCounter==1) {
			for(unsigned int rn=0; rn<refNames.size(); rn++) {
				if(refTypes.at(rn)=="Function") {
					FunctionReference *fr = new FunctionReference(name,expression,refNames.at(rn));
					s->addFunctionReference(fr);
				}
			}

		} else {
			cerr<<"!!!Error:  Functions can reference at most one other function!  Quitting."<<endl;
			return false;
		}




		return true;
	}


	//if we got here, we are creating a local function, so call the create local function function.
	return createLocalFunction(name, expression, argNames, refNames, refTypes, s, parameter, verbose);
}




bool createLocalFunction(string name,
			string expression,
			vector <string> &argNames,
			vector <string> &refNames,
			vector <string> &refTypes,
			System *s,
			map <string,double> &parameter,
			bool verbose)
	{





	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//@todo create a new function to handle these local functions
	//      handle local functions that reference another function (as in rateLaw1(x))


	cout<<"must be a local function..."<<endl;

	//Remember the original expression (useful for outputting)
	string originalExpression = expression;

	// First, extract out the parameter constants we will use, and make sure
	// that this isn't just a function reference to a local function...
	vector <string> paramNames;
	int otherFuncRefCounter=0;
	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Function") {
			otherFuncRefCounter++;
		} else if(refTypes.at(rn)=="Constant") {
			paramNames.push_back(refNames.at(rn));
		}
	}

	if(otherFuncRefCounter==1) {
		//handle the function reference...
		cout<<"handling local function reference!  (not yet handled)..."<<endl;
		return true;
	}
	if(otherFuncRefCounter!=0) {
		cerr<<"!!!Error:  Functions can reference at most one other function!  Quitting."<<endl;
		return false;
	}


	//some vectors to store results
	vector <string> obsUsedName;
	vector <string> obsUsedExpressionRef;
	vector <int> obsUsedScope;
	vector <bool> markForRemoval;

	//First, figure out the scope of the observables that we
	//referenced in the function.  This requires the following annoying
	//set of nested loops:
	for(unsigned int rn=0; rn<refNames.size(); rn++) {
		if(refTypes.at(rn)=="Observable") {
			string::size_type sPos=expression.find(refNames.at(rn));
			for( ; sPos!=string::npos; sPos=sPos=expression.find(refNames.at(rn),sPos+1)) {

				//first, add this reference to our vector lists
				obsUsedName.push_back(refNames.at(rn));
				obsUsedExpressionRef.push_back(refNames.at(rn));
				obsUsedScope.push_back(-1); //assume global scope (with -1) until we discover otherwise
				markForRemoval.push_back(false); //assume that we aren't removing this guy, yet

				//Now, go about finding the scope by finding the next enclosing parentheses...
				string::size_type openPar = expression.find_first_of('(',sPos);
				string::size_type closePar = expression.find_first_of(')',sPos);
				if(openPar!=string::npos && closePar!=string::npos) {
					if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at

						string possibleArg = expression.substr(openPar+1,closePar-openPar-1);
						NFutil::trim(possibleArg);

						for(unsigned int aIndex=0; aIndex<argNames.size(); aIndex++) {
							if(argNames.at(aIndex)==possibleArg) {
								//hurah!  we found the local scope of this guy

								//Now we have to gently excise the reference, and replace it with
								//a marker so we can uniquely identify it later...
								int lastIndex = obsUsedExpressionRef.size();
								string identifier = "_"+NFutil::toString(lastIndex);
								expression.replace(openPar,closePar-openPar+1,identifier);
								obsUsedExpressionRef.pop_back();
								obsUsedExpressionRef.push_back(obsUsedName.at(lastIndex-1)+identifier);

								obsUsedScope.pop_back();
								obsUsedScope.push_back(aIndex);
								break; //break cause we're done with this scope...
							}
						}
					}
				} //end if statement to find open and closed parentheses
			}//end for loop over sPos
		}
	} //end loop over the possible references


	//as an optimization, certain observables may be used in more than one position
    //with the same scope.  The loops above uniquely identified each instance of each
	//reference - but that means that the function will evaluate separately each
	//instance found in the expression.  If we have the same reference with the same
	//scope, then we can reduce by renaming them as a single reference.

	//so we loop over each used reference...
	for(unsigned int i=0; i<obsUsedExpressionRef.size(); i++) {

		//compare this reference to each of the other used references...
		for(unsigned int j=i+1; j<obsUsedExpressionRef.size(); j++) {

			//If we are already removing element j, then we already found it
			//and replaced it, so we can stop.
			if(markForRemoval.at(j)) continue;

			//cout<<"comparing: "<<obsUsedExpressionRef.at(i)<<" to "<< obsUsedExpressionRef.at(j)<<endl;
			if(obsUsedName.at(i)==obsUsedName.at(j)) {
				if(obsUsedScope.at(i)==obsUsedScope.at(j)) {
					//find it in the expression, replace it with the original, and mark it for removal
					string::size_type refPos =expression.find(obsUsedExpressionRef.at(j));
					expression.replace(refPos, obsUsedExpressionRef.at(j).size(),obsUsedExpressionRef.at(i));
					markForRemoval.at(j)=true;
				}
			}
		}
	}

	////////////////////
	// almost there...

	//Final vectors for storing the list after the removal (it is easier to just copy
	//than risk errors on indexing when we are removing from the original...)
	vector <string> finalObsUsedExpressionRef;
	vector <string> finalObsUsedName;
	vector <int> finalObsUsedScope;

	//Fill the final vectors
	for(unsigned int i=0; i<obsUsedExpressionRef.size(); i++) {
		if(!markForRemoval.at(i)) {
			finalObsUsedExpressionRef.push_back(obsUsedExpressionRef.at(i));
			finalObsUsedName.push_back(obsUsedName.at(i));
			finalObsUsedScope.push_back(obsUsedScope.at(i));
		}
	}


	//Ah.  and so we get here.  We now have:
	// 1) original functional expression (originalExpression)
	// 2) new revised functional expression with scope notation removed (expression)
	// 3) a list of each observable referenced and its referenced name (finalObsUsedExpressionRef)
	// 4) a list of what observables those references depend on (finalObsUsedName)
	// 5) a list of the scope of each of those references (finalObsUsedScope)
	// 6) the original list of arguments (argNames)
	// 7) the reduced list of parameter constants (paramNames)

	//so we can finally create our local function creator...



	// Functions needed: 1) LocalFuncCreator class
	//                   2) NFinput::CreateFreeCopyOfObs(obsName)
	//                   3) update system to handle refs to global or local functions




	cout<<endl<<endl<<".."<<endl;
	//exit(1);
}















////  New Function Parser
bool NFinput::initGlobalFunctions(
	TiXmlElement * pListOfFunctions,
	System * system,
	map <string,double> &parameter,
	bool verbose)
{
	try {
		vector <string> argNames;
		vector <string> refNames;
		vector <string> refTypes;

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
			if(verbose) cout<<"\t\tReading Function: "+funcName+"(";


			//Get the list of arguments for this function
			TiXmlElement *pListOfArgs = pFunction->FirstChildElement("ListOfArguments");
			if(pListOfArgs) {
				//Loop through each arguement
				TiXmlElement *pArg; bool firstArg = true;
				for ( pArg = pListOfArgs->FirstChildElement("Argument"); pArg != 0; pArg = pArg->NextSiblingElement("Argument"))
				{
					if(verbose && !firstArg) cout<<", ";
					firstArg=false;
					//Check again for errors by making sure the argument has a name
					string argName;
					if(!pArg->Attribute("id")) {
						if(verbose) cout<<" ?? ...\n";
						cerr<<"!!!Error:  Argument tag in Function: '" + funcName + "' must contain the id attribute.  Quitting."<<endl;
						return false;
					}

					argName = pArg->Attribute("id");

					argNames.push_back(argName);

					if(verbose) cout<<argName;
				}
			}
			if(verbose) cout<<")"<<endl;


			//Get the list of References
			TiXmlElement *pListOfRefs = pFunction->FirstChildElement("ListOfReferences");
			if(pListOfRefs) {

				//Loop through each reference
				TiXmlElement *pRef;
				for ( pRef = pListOfRefs->FirstChildElement("Reference"); pRef != 0; pRef = pRef->NextSiblingElement("Reference"))
				{
					//Check again for errors by making sure the parameter has a name
					if(!pRef->Attribute("name")) {
						cerr<<"!!!Error:  Reference tag in Function: '" + funcName + "' must contain a proper id.  Quitting."<<endl;
						return false;
					} else if(!pRef->Attribute("type")) {
						cerr<<"!!!Error:  Reference tag in Function: '" + funcName + "' must contain a proper type.  Quitting."<<endl;
						return false;
					}
//
					string refName = pRef->Attribute("name");
					string refType = pRef->Attribute("type");

					if(verbose) cout<<"Reference: "+refType+" "+refName<<endl;

					refNames.push_back(refName);
					refTypes.push_back(refType);
				}
			}


			//Read in the actual function definition
			string funcExpression = "";
			TiXmlElement *pExpression = pFunction->FirstChildElement("Expression");
			if(pExpression) {
				if(!pExpression->GetText()) {
					cerr<<"!!!Error:  Expression tag must actually contain a string for the function.  Quitting."<<endl;
										return false;
				}
				funcExpression = pExpression->GetText();
				//GlobalFunction *gf = new GlobalFunction(funcName, functionDefintion, argNames, argTypes, paramNames, paramValues);
				//if(!system->addGlobalFunction(gf)) {
				//	cerr<<"!!!Error:  Function name '"<<funcName<<"' has already been used.  You can't have two\n";
				//	cerr<<"functions with the same name, so I'll just stop now."<<endl;
				//	return false;
				//}
				if(verbose) cout<<"\t\t\t = "<<funcExpression<<endl;
			} else {
				cerr<<"!!!Error:  Expression tag for a function must exist!  Quitting."<<endl;
				return false;
			}

			//Here we actually generate the function or the function generator
			if(!createFunction(funcName,
					funcExpression,
					argNames,
					refNames,
					refTypes,
					system,
					parameter,
					verbose)) {
				return false;
			}

			//And here we clear our arrays
			argNames.clear();
			refNames.clear();
			refTypes.clear();
		}




		cout<<"done reading functions!"<<endl;
	//	exit(0);
		//Getting here means we read everything we could successfully
		return true;
	} catch (...) {
		//Uh oh! we got some unknown exception thrown, so we must abort!
		cerr<<"I caught some unknown error when I was trying to parse out a Global Function.\n";
		cerr<<"I'm at a loss for words right now, so you're on you're own."<<endl;
		return false;
	}
}










//////  Original Function Parser
//bool NFinput::initGlobalFunctions(
//	TiXmlElement * pListOfFunctions,
//	System * system,
//	map <string,double> &parameter,
//	bool verbose)
//{
//	try {
//		vector <string> argNames;
//		vector <string> argTypes;
//		vector <string> paramNames;
//		vector <double> paramValues;
//
//		//Loop through the Function tags...
//		TiXmlElement *pFunction;
//		for ( pFunction = pListOfFunctions->FirstChildElement("Function"); pFunction != 0; pFunction = pFunction->NextSiblingElement("Function"))
//		{
//			//Check if MoleculeType tag has a name...
//			if(!pFunction->Attribute("id")) {
//				cerr<<"!!!Error:  Function tag must contain the id attribute.  Quitting."<<endl;
//				return false;
//			}
//
//			//Read in the Function Name
//			string funcName = pFunction->Attribute("id");
//			if(verbose) cout<<"\t\tReading and Creating Function: "+funcName+"(";
//
//
//			//Get the list of arguments for this function
//			TiXmlElement *pListOfArgs = pFunction->FirstChildElement("ListOfArgs");
//			if(pListOfArgs) {
//				//Loop through each arguement
//				TiXmlElement *pArg; bool firstArg = true;
//				for ( pArg = pListOfArgs->FirstChildElement("Arg"); pArg != 0; pArg = pArg->NextSiblingElement("Arg"))
//				{
//					if(verbose && !firstArg) cout<<", ";
//					firstArg=false;
//					//Check again for errors by making sure the argument has a name
//					string argName, argType;
//					if(!pArg->Attribute("id") || !pArg->Attribute("name") || !pArg->Attribute("type")) {
//						if(verbose) cout<<" ?? ...\n";
//						cerr<<"!!!Error:  Arg tag in Function: '" + funcName + "' must contain the all attributes including id, name, and type.  Quitting."<<endl;
//						return false;
//					}
//
//					argName = pArg->Attribute("name");
//					argType = pArg->Attribute("type");
//					argNames.push_back(argName);
//					argTypes.push_back(argType);
//
//					if(verbose) cout<<argName;
//				}
//			}
//			if(verbose) cout<<")"<<endl;
//
//
//			//Get the list of Parameter Constants for this function
//			TiXmlElement *pListOfParamConst = pFunction->FirstChildElement("ListOfParameterConstants");
//			if(pListOfParamConst) {
//
//				//Loop through each arguement
//				TiXmlElement *pParamConst;
//				for ( pParamConst = pListOfParamConst->FirstChildElement("ParameterConstant"); pParamConst != 0; pParamConst = pParamConst->NextSiblingElement("ParameterConstant"))
//				{
//					//Check again for errors by making sure the parameter has a name
//					if(!pParamConst->Attribute("id")) {
//						cerr<<"!!!Error:  ParameterConstant tag in Function: '" + funcName + "' must contain a proper id.  Quitting."<<endl;
//						return false;
//					}
//
//					string paramName = pParamConst->Attribute("id");
//					if(parameter.find(paramName)==parameter.end()) {
//						cerr<<"Could not find parameter constant: "<<paramName<<" when creating function "<<funcName<<". Quitting"<<endl;
//						return false;
//					}
//					paramNames.push_back(paramName);
//					paramValues.push_back(parameter.find(paramName)->second);
//				}
//			}
//
//
//			//Read in the actual function definition
//			TiXmlElement *pDefinition = pFunction->FirstChildElement("Definition");
//			if(pDefinition) {
//				if(!pDefinition->GetText()) {
//					cerr<<"!!!Error:  Definition tag must actually contain a string for the function.  Quitting."<<endl;
//										return false;
//				}
//				string functionDefintion = pDefinition->GetText();
//				GlobalFunction *gf = new GlobalFunction(funcName, functionDefintion, argNames, argTypes, paramNames, paramValues);
//				if(!system->addGlobalFunction(gf)) {
//					cerr<<"!!!Error:  Function name '"<<funcName<<"' has already been used.  You can't have two\n";
//					cerr<<"functions with the same name, so I'll just stop now."<<endl;
//					return false;
//				}
//				if(verbose) cout<<"\t\t\t= "<<functionDefintion<<endl;
//			} else {
//				cerr<<"!!!Error:  Definition tag for a function must exist!  Quitting."<<endl;
//				return false;
//			}
//
//			argNames.clear();
//			argTypes.clear();
//			paramNames.clear();
//			paramValues.clear();
//		}
//
//
//
//		cout<<"done reading functions!"<<endl;
//		exit(0);
//		//Getting here means we read everything we could successfully
//		return true;
//	} catch (...) {
//		//Uh oh! we got some unknown exception thrown, so we must abort!
//		cerr<<"I caught some unknown error when I was trying to parse out a Global Function.\n";
//		cerr<<"I'm at a loss for words right now, so you're on you're own."<<endl;
//		return false;
//	}
//}
