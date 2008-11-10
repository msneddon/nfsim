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

LocalFunction::LocalFunction(string name,
					string funcString,
					vector <Observable *> &observables,
					vector <StateCounter *> &stateCounters,
					vector <string> &paramConstNames,
					vector <double> &paramConstValues) {

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


void LocalFunction::printDetails() {

	cout<<"Local Function: "+name+"()="<<funcString<<"\n";
	for(unsigned int i=0; i<n_obs; i++)
		cout<<" -Observable '"<<obs[i]->getAliasName()<<"' has local value: "<<obs[i]->getCount()<<"\n";
	for(unsigned int i=0; i<n_sc; i++)
		cout<<" -StateCounter '"<<sc[i]->name<<"' has local value: "<<sc[i]->value<<"\n";

	double val =FuncFactory::Eval(p);
	cout<<"Evaluates to: "<<val<<endl<<endl;
}


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
			cout<<"adding to sc"<<endl;
			sc[i]->add(*molIter);
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















