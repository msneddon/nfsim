#include "NFfunction.hh"



using namespace std;
using namespace NFcore;
using namespace mu;


GlobalFunction::GlobalFunction(string name,
		string funcString,
		vector <string> &argNames,
		vector <string> &argTypes,
		vector <string> &paramConstNames,
		vector <double> &paramConstValues)
{
	if(argNames.size()!=argTypes.size()) {
		cerr<<"Trying to create a global function, but your argument vectors don't match up!"<<endl;
		cerr<<"Quitting!"<<endl;
		exit(1);
	}

	if(paramConstNames.size()!=paramConstValues.size()) {
		cerr<<"Trying to create a global function, but your parameter vectors don't match up!"<<endl;
		cerr<<"Quitting!"<<endl;
		exit(1);
	}

	this->name = name;
	this->funcString = funcString;

	this->n_args=argNames.size();
	this->argNames = new string[n_args];
	this->argTypes = new string[n_args];
	for(unsigned int a=0; a<n_args; a++)
	{
		this->argNames[a]=argNames.at(a);
		this->argTypes[a]=argTypes.at(a);
	}

	this->n_paramConst=paramConstNames.size();
	this->paramNames = new string[n_paramConst];
	this->paramValues = new double[n_paramConst];
	for(unsigned int i=0; i<n_paramConst; i++)
	{
		this->paramNames[i]=paramConstNames.at(i);
		this->paramValues[i]=paramConstValues.at(i);
	}
	p=0;
}



GlobalFunction::~GlobalFunction()
{
	delete [] argNames;
	delete [] argTypes;
	delete [] paramNames;
	delete [] paramValues;
	if(p!=NULL) delete p;
}




void GlobalFunction::prepareForSimulation(System *s)
{
	try {
		p=FuncFactory::create();
		for(unsigned int a=0; a<n_args; a++)
		{
			if(argTypes[a]=="MoleculeObservable") {
				Observable *obs = s->getObservableByName(argNames[a]);
				if(obs==NULL) {
					cout<<"When creating global function: "<<this->name<<endl<<" could not find the observable: ";
					cout<<argNames[a]<<" of type "<<argTypes[a]<<endl;
					cout<<"Quitting."<<endl;
					exit(1);
				}
				obs->addReferenceToMyself(p);
			} else {
				cout<<"Uh oh, an unrecognized argType ("<<argTypes[a]<<") for a function! "<<argNames[a]<<endl;
				cout<<"Try using the type: \"MoleculeObservable\""<<endl;
				cout<<"Quitting because this will give unpredicatable results, or just crash."<<endl;
				exit(1);
			}
		}

		for(unsigned int i=0; i<n_paramConst; i++)
		{
			p->DefineConst(paramNames[i],paramValues[i]);
		}
		p->SetExpr(this->funcString);
	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error preparing function "<<name<<" in class GlobalFunction!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}
}


void GlobalFunction::attatchRxn(ReactionClass *r)
{
	//unsigned int n_rxns;
	//ReactionClass *rxns;

}





void GlobalFunction::printDetails()
{
	cout<<"--------"<<endl;
	cout<<"Details of Function: '"<< this->name << "'"<<endl;
	cout<<" ="<<funcString<<endl;
	cout<<"   -Arguments:"<<endl;
	for(unsigned int a=0; a<n_args; a++)
		cout<<"     "<<argTypes[a]<<":  "<<argNames[a]<<endl;
	cout<<"   -Constant Parameters:"<<endl;
	for(unsigned int i=0; i<n_args; i++)
		cout<<"     "<<paramNames[i]<<" with value of  "<<paramValues[i]<<endl;





//	// Get the map with the variables
//	Parser::varmap_type variables = p->GetVar();
//	cout << (int)variables.size() << " variables."<<endl;
//	mu::Parser::varmap_type::const_iterator item = variables.begin();
//
//	// Query the variables
//	for (; item!=variables.end(); ++item)
//	{
//	  cout << "  Name: " << item->first << " Address: [0x" << item->second << "]  Value: "<< *(item->second)<<"\n";
//	}
//

	cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;

}







StateCounter::StateCounter(string name, MoleculeType *mt, string stateName) {
	this->name=name;
	this->mt = mt;
	this->stateIndex = mt->getCompIndexFromName(stateName);
	this->value=0;

	//Make sure this is a state we can count on!
	if(!mt->isIntegerComponent(stateName)) {
		//if it is not an integer component state, we must abort because
		//we can not evaluate the sum of a non integer component
		cerr<<"Trying to create a stateCounter: '"<<name<<"' on the state: '"<<stateName<<"'\n";
		cerr<<"of MoleculeType: '"<<mt->getName()<<"', but the state you have selected cannot\n";
		cerr<<"have integer values, so summations on this state are undefined.  I am quitting."<<endl;
		exit(1);
	}
}
StateCounter::~StateCounter() {
	mt=0;
}

void StateCounter::add(Molecule *m) {
	if(m->getMoleculeType()==mt) {

		cout<<"matched moleculeType"<<endl;
		value+=m->getComponentState(stateIndex);
		cout<<"found component state: "<<m->getComponentState(stateIndex)<<endl;
		cout<<"updating value to: "<< value<<endl;
	}
}





/////////////////////////////////////////

LocalFunction::LocalFunction(System *s,
					string name,
					string funcString,
					vector <Observable *> &observables,
					vector <StateCounter *> &stateCounters,
					vector <string> &paramConstNames,
					vector <double> &paramConstValues) {

	cout<<"Attempting to create local function: "<<name<<endl;

	if(paramConstNames.size()!=paramConstValues.size()) {
		cerr<<"Trying to create a local function, but your parameter vectors don't match up!"<<endl;
		cerr<<"Quitting!"<<endl;
		exit(1);
	}

	//Extract out the basic information about this function
	this->name = name;
	this->funcString = funcString;

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




	// finally add the function to the system so it is remembered...
	s->addLocalFunction(this);
}

LocalFunction::~LocalFunction() {
	for(unsigned int i=0; i<n_obs; i++)
		delete obs[i];
	delete [] obs;

	for(unsigned int i=0; i<n_sc; i++)
		delete sc[i];
	delete [] sc;

	delete [] paramNames;
	delete [] paramValues;
	if(p!=NULL) delete p;
}



//This function is generally called by a DOR reaction class once the
//DOR reaction class has established that the value of this function
//is required for some moleculetype...
void LocalFunction::addTypeIMoleculeDependency(MoleculeType *mt) {

	//First, make sure we haven't added this bad boy yet
	for(unsigned int i=0; i<this->typeI_mol.size(); i++) {
		if(typeI_mol.at(i)==mt) return;
	}

	//First, add myself to the moleculeType
	int index = mt->addLocalFunc_TypeI(this);
	this->typeI_mol.push_back(mt);
	this->typeI_localFunctionIndex.push_back(index);
}


void LocalFunction::printDetails() {

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
}


//Note: not the most effecient function in the world, but it does what it has to for now
//we should make this faster in the future...
double LocalFunction::evaluateOn(Molecule *m) {

	list <Molecule *> molList;
	list <Molecule *>::iterator molIter;

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

	//Update the molecules that needed this function evaluated...
	for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
		for(unsigned int ti=0; ti<typeI_mol.size(); ti++) {
			if((*molIter)->getMoleculeType()==typeI_mol.at(ti)) {
				(*molIter)->setLocalFunctionValue(newValue,this->typeI_localFunctionIndex.at(ti));
			}
		}
	}



	return FuncFactory::Eval(p);
}


//			//List of observables that this local function depends on
//			Observable ** observables;
//			int * observableValues;
//
//
//			string name;
//			string funcString;
//
//			unsigned int n_args;
//			string *argNames;
//			string *argTypes;
//
//			unsigned int n_paramConst;
//			string *paramNames;
//			double *paramValues;















