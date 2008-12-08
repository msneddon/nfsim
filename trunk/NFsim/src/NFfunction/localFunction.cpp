/*
 * localFunction.cpp
 *
 *  Created on: Dec 4, 2008
 *      Author: msneddon
 */


#include "NFfunction.hh"



using namespace std;
using namespace NFcore;
using namespace mu;











string LocalFunction::getName() const {
	return this->name;
}
string LocalFunction::getNiceName() const {
	return nicename;
}
string LocalFunction::getExpression() const {
	return originalExpression;
}
string LocalFunction::getParsedExpression() const {
	return parsedExpression;
}









LocalFunction::LocalFunction(System *s,
					string name,
					string originalExpression,
					string parsedExpression,
					vector <string> &args,
					vector <string> &varRefNames,
					vector <string> &varObservableNames,
					vector <Observable *> & varObservables,
					vector <int> &varRefScope,
					vector <string> paramNames)
{
	cout<<"Attempting to create local function: "<<name<<endl;

	if(args.size()>1) {
		cerr<<"For efficiency, local functions currently support a maximum of 1 argument."<<endl;
		cerr<<"Quitting now."<<endl;
		exit(1);
	}
	if(args.size()<1) {
		cerr<<"When creating a local Function, ERROR!! there were no args, so it is a global function."<<endl;
		cerr<<"Quitting now."<<endl;
		exit(1);
	}

	//Do the basics first...
	this->name = name;
	this->originalExpression=originalExpression;
	this->parsedExpression=parsedExpression;


	//Move the vectors into our neat little arrays
	this->n_args=args.size();
	this->argNames = new string[n_args];
	for(unsigned int i=0; i<n_args; i++) {
		this->argNames[i]=args.at(i);
	}

	this->n_varRefs=varRefNames.size();
	this->varRefNames = new string[n_varRefs];
	this->varObservableNames = new string[n_varRefs];
	this->varLocalObservables = new Observable *[n_varRefs];
	this->varRefScope = new int[n_varRefs];
	for(unsigned int i=0; i<n_varRefs; i++) {
		this->varRefNames[i] = varRefNames.at(i);
		this->varObservableNames[i] = varObservableNames.at(i);
		this->varLocalObservables[i] = varObservables.at(i);
		this->varRefScope[i] = varRefScope.at(i);
	}

	this->n_params=paramNames.size();
	this->paramNames = new string[n_params];
	for(unsigned int i=0; i<n_params; i++) {
		this->paramNames[i] = paramNames.at(i);
	}

	//now assemble the nicename
	nicename = this->name + "(";
	for(unsigned int i=0;i<n_args; i++) {
		if(i==0) nicename+=argNames[i];
		else nicename+=","+argNames[i];
	}
	nicename+=")";


	cout<<nicename<<endl;



	p=0;


/*
		//Grab the observables that the function requires
		n_obs=observables.size();
		this->obs = new Observable * [n_obs];
		this->obsVal = new int [n_obs];
		for(unsigned int i=0; i<n_obs; i++)
			this->obs[i]=observables.at(i);

		//Also grab the stateCounters
		this->n_sc=stateCounters.size();
		this->sc = new StateCounter * [n_sc];
		for(unsigned int i=0; i<n_sc; i++)
			this->sc[i]=stateCounters.at(i);


		//Grab the parameters that the function requires
		this->n_paramConst=paramConstNames.size();
		this->paramNames = new string[n_paramConst];
		this->paramValues = new double[n_paramConst];
		for(unsigned int i=0; i<n_paramConst; i++)
		{
			this->paramNames[i]=paramConstNames.at(i);
			this->paramValues[i]=paramConstValues.at(i);
		}


		//Finally, we can create the local function
		try {
			p=FuncFactory::create();

			//Point the function to the observable
			for(unsigned int i=0; i<n_obs; i++) {
				obs[i]->addReferenceToMyself(p);
			}

			//Point the function to the counter
			for(unsigned int i=0; i<n_sc; i++) {
				p->DefineVar(sc[i]->name,&(sc[i]->value));
			}

			//Set the constant values
			for(unsigned int i=0; i<n_paramConst; i++) {
				p->DefineConst(paramNames[i],paramValues[i]);
			}

			//Finally, we can set the expression
			p->SetExpr(this->funcString);

		//Catch anything that goes astray
		} catch (mu::Parser::exception_type &e) {
			cout<<"Error creating local function "<<name<<" in class LocalFunction!!  This is what happened:"<<endl;
			cout<< "  "<<e.GetMsg() << endl;
			cout<<"Quitting."<<endl;
			exit(1);
		}


		////////////////////////////
		//If everything works, we can now identify the type II molecule types, because we have all the observables
		//used to define this local function.  So let's add them now..
		vector <TemplateMolecule *> tmList;
		vector <MoleculeType *> addedMoleculeTypes;

		cout<<"Now remembering type II molecules..."<<endl;
		bool hasAdded = false;
		for(unsigned int i=0; i<n_obs; i++) {
			TemplateMolecule::traverse(this->obs[i]->getTemplateMolecule(),tmList);
			cout<<"traversed obs "<<i<<" and found: "<<tmList.size()<<" templates\n";

			for(unsigned int t=0; t<tmList.size(); t++) {
				MoleculeType *mt = tmList.at(t)->getMoleculeType();

				//Make sure we haven't added this molecule type before
				hasAdded = false;
				for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
					if(addedMoleculeTypes.at(m)==mt) {
						hasAdded=true; break;
					}
				}
				if(!hasAdded) {
					addedMoleculeTypes.push_back(mt);
					cout<<"remembering: "<<mt->getName()<<endl;
				} else {
					cout<<"ignoring: "<<mt->getName()<<endl;
				}
			}
			tmList.clear();
		}

		for(unsigned int i=0; i<n_sc; i++) {
			MoleculeType *mt = this->sc[i]->mt;

			//Make sure we haven't added this molecule type before
			hasAdded = false;
			for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
				if(addedMoleculeTypes.at(m)==mt) {
					hasAdded=true; break;
				}
			}
			if(!hasAdded) {
				addedMoleculeTypes.push_back(mt);
				cout<<"*remembering: "<<mt->getName()<<endl;
			} else {
				cout<<"*ignoring: "<<mt->getName()<<endl;
			}
		}

		for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
			int index = addedMoleculeTypes.at(m)->addLocalFunc_TypeII(this);
			this->typeII_mol.push_back(addedMoleculeTypes.at(m));
			this->typeII_localFunctionIndex.push_back(index);
		}



		//evaluation level is the degree to which the local function is evaluated
		//An evaluation level of 0 evaluates over the entire species.  A level of
		// 1 searches the immediate molecule only.  All other levels are currently
		//not supported....  default evaluates over the entire species
		evaluationLevel = 0;


		// finally add the function to the system so it is remembered...
		s->addLocalFunction(this);

*/

}





void LocalFunction::printDetails(System *s)
{
	cout<<"Local Function: "+this->nicename+"\n";
	cout<<" = "<<this->originalExpression<<endl;
	cout<<" parsed expression = "<<this->parsedExpression<<endl;

	cout<<"   -Variable References:"<<endl;
	for(unsigned int i=0; i<n_varRefs; i++) {
		if(varRefScope[i]==-1) {
			cout<<"         "<<varObservableNames[i]<<" (scope=global): ";
			cout<<s->getObservableByName(varRefNames[i])->getCount()<<endl;
		} else {
			cout<<"         "<<varObservableNames[i]<<" (scope=";
			cout<<argNames[varRefScope[i]]<<") last evaluated to: ";
			cout<<varLocalObservables[i]->getCount()<<endl;
		}
	}

	if(n_params>0) {
		cout<<"   -Constant Parameters:"<<endl;
		for(unsigned int i=0; i<n_params; i++) {
			cout<<"         "<<paramNames[i]<<" = " << s->getParameter(paramNames[i])<<endl;
		}
	}

	if(p!=0)
		cout<<"   Function last evaluated to: "<<FuncFactory::Eval(p)<<endl;
}






void LocalFunction::prepareForSimulation(System *s) {

	cout<<"preparing local function: "<<this->nicename<<endl;

	//Finally, we can create the local function
		try {
			p=FuncFactory::create();

			//Give the local observable to the function so it can be used
			for(unsigned int i=0; i<n_varRefs; i++) {
				if(this->varRefScope[i]==-1) { //for global variables, use the global observable
					s->getObservableByName(this->varObservableNames[i])->addReferenceToMyself(this->varRefNames[i],p);
				} else { //for local observables, use this function's observable
					this->varLocalObservables[i]->addReferenceToMyself(this->varRefNames[i],p);
				}
			}
	//
			//Point the function to the counter
	//		for(unsigned int i=0; i<n_sc; i++) {
	//			p->DefineVar(sc[i]->name,&(sc[i]->value));
	//		}
	//
			//Set the constant values
			for(unsigned int i=0; i<this->n_params; i++) {
				p->DefineConst(this->paramNames[i],s->getParameter(paramNames[i]));
			}

			//Finally, we can set the expression
			p->SetExpr(this->parsedExpression);

		//Catch anything that goes astray
		} catch (mu::Parser::exception_type &e) {
			cout<<"Error creating local function "<<name<<" in class LocalFunction!!  This is what happened:"<<endl;
			cout<< "  "<<e.GetMsg() << endl;
			cout<<"Quitting."<<endl;
			exit(1);
		}
}

double LocalFunction::evaluateOn(Molecule *m, int scope) {

	cout<<"evaluating local function: "<<this->nicename<<endl;


//	if(evaluationLevel==0) {
//		list <Molecule *> molList;
//		list <Molecule *>::iterator molIter;
//
//		//Get the species.  If we are using complex bookkeeping, we should take advantage
//		//of that here, although I don't here yet.
//		m->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);
//		for(unsigned int i=0; i<n_obs; i++) obs[i]->clear();
//		for(unsigned int i=0; i<n_sc; i++) sc[i]->reset();
//
//		//Update the observables and counters, as something has changed
//		for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
//			for(unsigned int i=0; i<n_obs; i++)
//				if(obs[i]->isObservable(*molIter)) obs[i]->straightAdd();
//			for(unsigned int i=0; i<n_sc; i++) {
//				//cout<<"adding to sc"<<endl;
//				sc[i]->add(*molIter);
//			}
//		}
//
//		//evaluate the function
//		double newValue = FuncFactory::Eval(p);
//
//		//Update the molecules (Type I) that needed this function evaluated...
//		for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
//			for(unsigned int ti=0; ti<typeI_mol.size(); ti++) {
//				if((*molIter)->getMoleculeType()==typeI_mol.at(ti)) {
//					(*molIter)->setLocalFunctionValue(newValue,this->typeI_localFunctionIndex.at(ti));
//					(*molIter)->updateDORRxnValues();
//				}
//			}
//		}
//		return newValue;
//
//	//Evaluation level of 1 means that we search only this molecule when
//	//evaluating the function - the steps are the same as for the entire species, except we evaluate
//	//only on the molecule
//	} else if(evaluationLevel==1) {
//
//
//		for(unsigned int i=0; i<n_obs; i++) obs[i]->clear();
//		for(unsigned int i=0; i<n_sc; i++) sc[i]->reset();
//
//		for(unsigned int i=0; i<n_obs; i++) {
//			if(obs[i]->isObservable(m)) obs[i]->straightAdd();
//		}
//
//
//
//		for(unsigned int i=0; i<n_sc; i++) {
//			sc[i]->add(m);
//		}
//
//		double newValue = FuncFactory::Eval(p);
//		for(unsigned int ti=0; ti<typeI_mol.size(); ti++) {
//			if(m->getMoleculeType()==typeI_mol.at(ti)) {
//				m->setLocalFunctionValue(newValue,this->typeI_localFunctionIndex.at(ti));
//				m->updateDORRxnValues();
//			}
//		}
//		return newValue;
//
//
//	} else {
//		cout<<"Internal error in LocalFunction::evaluateOn()! trying to evaluate a function with a bad evaluation level."<<endl;
//		exit(1);
//	}
//
//	//Should never get here...
//	return 0;



	return -1;
}


//
//
//LocalFunction::LocalFunction(System *s,
//					string name,
//					string funcString,
//					vector <Observable *> &observables,
//					vector <StateCounter *> &stateCounters,
//					vector <string> &paramConstNames,
//					vector <double> &paramConstValues) {
//
//	cout<<"Attempting to create local function: "<<name<<endl;
//
//	if(paramConstNames.size()!=paramConstValues.size()) {
//		cerr<<"Trying to create a local function, but your parameter vectors don't match up!"<<endl;
//		cerr<<"Quitting!"<<endl;
//		exit(1);
//	}
//
//	//Extract out the basic information about this function
//	this->name = name;
//	this->funcString = funcString;
//
//	//Grab the observables that the function requires
//	n_obs=observables.size();
//	this->obs = new Observable * [n_obs];
//	this->obsVal = new int [n_obs];
//	for(unsigned int i=0; i<n_obs; i++)
//		this->obs[i]=observables.at(i);
//
//	//Also grab the stateCounters
//	this->n_sc=stateCounters.size();
//	this->sc = new StateCounter * [n_sc];
//	for(unsigned int i=0; i<n_sc; i++)
//		this->sc[i]=stateCounters.at(i);
//
//
//	//Grab the parameters that the function requires
//	this->n_paramConst=paramConstNames.size();
//	this->paramNames = new string[n_paramConst];
//	this->paramValues = new double[n_paramConst];
//	for(unsigned int i=0; i<n_paramConst; i++)
//	{
//		this->paramNames[i]=paramConstNames.at(i);
//		this->paramValues[i]=paramConstValues.at(i);
//	}
//
//
//	//Finally, we can create the local function
//	try {
//		p=FuncFactory::create();
//
//		//Point the function to the observable
//		for(unsigned int i=0; i<n_obs; i++) {
//			obs[i]->addReferenceToMyself(p);
//		}
//
//		//Point the function to the counter
//		for(unsigned int i=0; i<n_sc; i++) {
//			p->DefineVar(sc[i]->name,&(sc[i]->value));
//		}
//
//		//Set the constant values
//		for(unsigned int i=0; i<n_paramConst; i++) {
//			p->DefineConst(paramNames[i],paramValues[i]);
//		}
//
//		//Finally, we can set the expression
//		p->SetExpr(this->funcString);
//
//	//Catch anything that goes astray
//	} catch (mu::Parser::exception_type &e) {
//		cout<<"Error creating local function "<<name<<" in class LocalFunction!!  This is what happened:"<<endl;
//		cout<< "  "<<e.GetMsg() << endl;
//		cout<<"Quitting."<<endl;
//		exit(1);
//	}
//
//
//	////////////////////////////
//	//If everything works, we can now identify the type II molecule types, because we have all the observables
//	//used to define this local function.  So let's add them now..
//	vector <TemplateMolecule *> tmList;
//	vector <MoleculeType *> addedMoleculeTypes;
//
//	cout<<"Now remembering type II molecules..."<<endl;
//	bool hasAdded = false;
//	for(unsigned int i=0; i<n_obs; i++) {
//		TemplateMolecule::traverse(this->obs[i]->getTemplateMolecule(),tmList);
//		cout<<"traversed obs "<<i<<" and found: "<<tmList.size()<<" templates\n";
//
//		for(unsigned int t=0; t<tmList.size(); t++) {
//			MoleculeType *mt = tmList.at(t)->getMoleculeType();
//
//			//Make sure we haven't added this molecule type before
//			hasAdded = false;
//			for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
//				if(addedMoleculeTypes.at(m)==mt) {
//					hasAdded=true; break;
//				}
//			}
//			if(!hasAdded) {
//				addedMoleculeTypes.push_back(mt);
//				cout<<"remembering: "<<mt->getName()<<endl;
//			} else {
//				cout<<"ignoring: "<<mt->getName()<<endl;
//			}
//		}
//		tmList.clear();
//	}
//
//	for(unsigned int i=0; i<n_sc; i++) {
//		MoleculeType *mt = this->sc[i]->mt;
//
//		//Make sure we haven't added this molecule type before
//		hasAdded = false;
//		for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
//			if(addedMoleculeTypes.at(m)==mt) {
//				hasAdded=true; break;
//			}
//		}
//		if(!hasAdded) {
//			addedMoleculeTypes.push_back(mt);
//			cout<<"*remembering: "<<mt->getName()<<endl;
//		} else {
//			cout<<"*ignoring: "<<mt->getName()<<endl;
//		}
//	}
//
//	for(unsigned int m=0; m<addedMoleculeTypes.size(); m++) {
//		int index = addedMoleculeTypes.at(m)->addLocalFunc_TypeII(this);
//		this->typeII_mol.push_back(addedMoleculeTypes.at(m));
//		this->typeII_localFunctionIndex.push_back(index);
//	}
//
//
//
//	//evaluation level is the degree to which the local function is evaluated
//	//An evaluation level of 0 evaluates over the entire species.  A level of
//	// 1 searches the immediate molecule only.  All other levels are currently
//	//not supported....  default evaluates over the entire species
//	evaluationLevel = 0;
//
//
//	// finally add the function to the system so it is remembered...
//	s->addLocalFunction(this);
//}

LocalFunction::~LocalFunction() {
/*	for(unsigned int i=0; i<n_obs; i++)
		delete obs[i];
	delete [] obs;

	for(unsigned int i=0; i<n_sc; i++)
		delete sc[i];
	delete [] sc;

	delete [] paramNames;
	delete [] paramValues;*/
	if(p!=NULL) delete p;
}


/*void LocalFunction::setEvaluationLevel(int eLevel) {

	if(eLevel<0 || eLevel>1) {
		cout<<"Error when setting evaluation level of function: "<<getNiceName();
		cout<<"\nEvaluation level given was:"<<eLevel<<" but currently only supports levels of 0 or 1."<<endl;
		exit(1);
	}

	this->evaluationLevel = eLevel;
}*/


//This function is generally called by a DOR reaction class once the
//DOR reaction class has established that the value of this function
//is required for some moleculetype...
/*void LocalFunction::addTypeIMoleculeDependency(MoleculeType *mt) {

	//First, make sure we haven't added this bad boy yet
	for(unsigned int i=0; i<this->typeI_mol.size(); i++) {
		if(typeI_mol.at(i)==mt) return;
	}

	//First, add myself to the moleculeType
	int index = mt->addLocalFunc_TypeI(this);
	this->typeI_mol.push_back(mt);
	this->typeI_localFunctionIndex.push_back(index);
}*/

/*
int LocalFunction::getIndexOfTypeIFunctionValue(Molecule *m) {

	for(unsigned int i=0; i<this->typeI_mol.size(); i++) {
		if(typeI_mol.at(i)==m->getMoleculeType()) return this->typeI_localFunctionIndex.at(i);
	}
	cout<<"Error when getting the index of a Type I function value in LocalFunction:"<<endl;
	cout<<"Could not find the molecule type: '"<<m->getMoleculeType()->getName()<<"' as a type I molecule of this function: "<<this->getNiceName()<<endl;
	exit(1);
}
*/

//void LocalFunction::printDetails() {
/*
	cout<<"Local Function: "+name+"()="<<funcString<<"\n";
	for(unsigned int i=0; i<n_obs; i++)
		cout<<" -Observable '"<<obs[i]->getAliasName()<<"' has local value: "<<obs[i]->getCount()<<"\n";
	for(unsigned int i=0; i<n_sc; i++)
		cout<<" -StateCounter '"<<sc[i]->name<<"' has local value: "<<sc[i]->value<<"\n";

	for(unsigned int i=0; i<typeII_mol.size(); i++) {
		cout<<" -TypeII: "<<typeII_mol.at(i)->getName()<<" @ index "<<typeII_localFunctionIndex.at(i)<<"\n";
	}

	double val =FuncFactory::Eval(p);
	cout<<"  Last Evaluated to: "<<val<<endl<<endl;
	*/
//}


//Note: not the most effecient function in the world, but it does what it has to for now
//we should make this faster in the future...
/*double LocalFunction::evaluateOn(Molecule *m) {

	if(evaluationLevel==0) {
		list <Molecule *> molList;
		list <Molecule *>::iterator molIter;

		//Get the species.  If we are using complex bookkeeping, we should take advantage
		//of that here, although I don't here yet.
		m->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);
		for(unsigned int i=0; i<n_obs; i++) obs[i]->clear();
		for(unsigned int i=0; i<n_sc; i++) sc[i]->reset();

		//Update the observables and counters, as something has changed
		for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
			for(unsigned int i=0; i<n_obs; i++)
				if(obs[i]->isObservable(*molIter)) obs[i]->straightAdd();
			for(unsigned int i=0; i<n_sc; i++) {
				//cout<<"adding to sc"<<endl;
				sc[i]->add(*molIter);
			}
		}

		//evaluate the function
		double newValue = FuncFactory::Eval(p);

		//Update the molecules (Type I) that needed this function evaluated...
		for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
			for(unsigned int ti=0; ti<typeI_mol.size(); ti++) {
				if((*molIter)->getMoleculeType()==typeI_mol.at(ti)) {
					(*molIter)->setLocalFunctionValue(newValue,this->typeI_localFunctionIndex.at(ti));
					(*molIter)->updateDORRxnValues();
				}
			}
		}
		return newValue;

	//Evaluation level of 1 means that we search only this molecule when
	//evaluating the function - the steps are the same as for the entire species, except we evaluate
	//only on the molecule
	} else if(evaluationLevel==1) {


		for(unsigned int i=0; i<n_obs; i++) obs[i]->clear();
		for(unsigned int i=0; i<n_sc; i++) sc[i]->reset();

		for(unsigned int i=0; i<n_obs; i++) {
			if(obs[i]->isObservable(m)) obs[i]->straightAdd();
		}



		for(unsigned int i=0; i<n_sc; i++) {
			sc[i]->add(m);
		}

		double newValue = FuncFactory::Eval(p);
		for(unsigned int ti=0; ti<typeI_mol.size(); ti++) {
			if(m->getMoleculeType()==typeI_mol.at(ti)) {
				m->setLocalFunctionValue(newValue,this->typeI_localFunctionIndex.at(ti));
				m->updateDORRxnValues();
			}
		}
		return newValue;


	} else {
		cout<<"Internal error in LocalFunction::evaluateOn()! trying to evaluate a function with a bad evaluation level."<<endl;
		exit(1);
	}

	//Should never get here...
	return 0;
}*/
