#include "NFfunction.hh"



using namespace std;
using namespace NFcore;
using namespace mu;


GlobalFunction::GlobalFunction(string name,
		string funcExpression,
		vector <string> &varRefNames,
		vector <string> &varRefTypes,
		vector <string> &paramNames,
		System *s)
{
	if(varRefNames.size()!=varRefTypes.size()) {
		cerr<<"Trying to create a global function, but your variable reference vectors don't match up in size!"<<endl;
		cerr<<"Quitting!"<<endl;
		exit(1);
	}

	this->name = name;
	this->funcExpression = funcExpression;

	this->n_varRefs=varRefNames.size();
	this->varRefNames = new string[n_varRefs];
	this->varRefTypes = new string[n_varRefs];
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		this->varRefNames[vr]=varRefNames.at(vr);
		this->varRefTypes[vr]=varRefTypes.at(vr);
	}

	this->n_params=paramNames.size();
	this->paramNames = new string[n_params];
	for(unsigned int i=0; i<n_params; i++) {
		this->paramNames[i]=paramNames.at(i);
	}
	p=0;
}



GlobalFunction::~GlobalFunction()
{
	delete [] varRefNames;
	delete [] varRefTypes;
	delete [] paramNames;
	if(p!=NULL) delete p;
}




void GlobalFunction::prepareForSimulation(System *s)
{
	try {
		p=FuncFactory::create();
		for(unsigned int vr=0; vr<n_varRefs; vr++)
		{
			if(varRefTypes[vr]=="Observable") {
				Observable *obs = s->getObservableByName(varRefNames[vr]);
				if(obs==NULL) {
					cout<<"When creating global function: "<<this->name<<endl<<" could not find the observable: ";
					cout<<varRefNames[vr]<<" of type "<<varRefTypes[vr]<<endl;
					cout<<"Quitting."<<endl;
					exit(1);
				}
				obs->addReferenceToMyself(p);
			} else {
				cout<<"here"<<endl;
				cout<<"Uh oh, an unrecognized argType ("<<varRefTypes[vr]<<") for a function! "<<varRefNames[vr]<<endl;
				cout<<"Try using the type: \"MoleculeObservable\""<<endl;
				cout<<"Quitting because this will give unpredicatable results, or just crash."<<endl;
				exit(1);
			}
		}

		for(unsigned int i=0; i<n_params; i++) {
			p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
		}
		p->SetExpr(this->funcExpression);

	}
	catch (mu::Parser::exception_type &e)
	{
		cout<<"Error preparing function "<<name<<" in class GlobalFunction!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}
}

void GlobalFunction::updateParameters(System *s) {
	//cout<<"Updating parameters for function: "<<name<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
	}

}




void GlobalFunction::attatchRxn(ReactionClass *r)
{
	//unsigned int n_rxns;
	//ReactionClass *rxns;

}





void GlobalFunction::printDetails()
{
	cout<<"Global Function: '"<< this->name << "()'"<<endl;
	cout<<" ="<<funcExpression<<endl;
	cout<<"   -Variable References:"<<endl;
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		cout<<"         "<<varRefTypes[vr]<<":  "<<varRefNames[vr]<<" = " << ""<<endl;
	}
	cout<<"   -Constant Parameters:"<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		cout<<"         "<<paramNames[i]<<endl;
	}





//	// Get the map with the variables
//	mu::Parser::varmap_type variables = p->GetVar();
//	cout << (int)variables.size() << " variables."<<endl;
//	mu::Parser::varmap_type::const_iterator item = variables.begin();
//	// Query the variables
//	for (; item!=variables.end(); ++item)
//	{
//	  cout << "  Name: " << item->first << " Address: [0x" << item->second << "]  Value: "<< *(item->second)<<"\n";
//	}


	if(p!=0) {
		// AS-2021
		if (this->fileFunc) {
			cout<<"   Function currently evaluates to: "<<this->fileEval()<<endl;
		} else {
			cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
		}
		// cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
		// AS-2021
	}
}

// AS-2021
void GlobalFunction::loadParamFile(string filePath) 
{
	// cout<<"######## loading file: "<<filePath<<endl;
	// cout<<"######## for function: "<<paramName<<endl;
	// setup our vectors
	vector <double> time;
	vector <double> values;
	// open file for reading
	ifstream file(filePath.c_str());
	// strings for looping over the file
	string line, word, content;
	string a,b;
	while (file >> a >> b) {
		// convert a to double
		istringstream aos(a);
		double d;
		aos >> d;
		// add it to time
		time.push_back(d);
		// convert b to double
		istringstream bos(b);
		bos >> d;
		// add it to values
		values.push_back(d);
	}
	// put the vectors into data vector
	this->data.push_back(time);
	this->data.push_back(values);
	return;
};

void GlobalFunction::addCounterPointer(double *counter){
	this->counter = counter;
}

void GlobalFunction::enableFileDependency(string filePath) {
	// load file
	this->loadParamFile(filePath);
	// we just want to keep a record of this
	this->filePath = filePath;
	// this sets it up so that this function knows it's supposed
	// to be pulling values from a file
	this->fileFunc = true;
}

double GlobalFunction::fileEval() {
	// TODO: Error checking and reporting
	// initialize index
	int ctrInd = 0;
	// distance from current value
	double cdist;
	// distance from new value to compare to current
	double ndist;
	// counter val
	double ctrVal = (*this->counter);
	// find the index closest in time
	for (int i=0;i<this->data[0].size();i++) {
		// if it's the first value, calculate and move on
		if(i==0) {
			cdist = abs(ctrVal-this->data[0][i]);
		} else {
			// calculate the new distance
			ndist = abs(ctrVal-this->data[0][i]);
			// if the distance is increasing, we are done
			// this relies on the assumption that the original
			// time series is ordered
			if(ndist>cdist) {
				break;
			} else {
				// we are getting closer to the value we want
				// so we keep going
				cdist = ndist;
				ctrInd += 1;
			}
		}
	}
	// index can't be larger than the array size
	if (ctrInd>=this->data[0].size()) {
		ctrInd = this->data[0].size()-1;
	}
	// debug stuff
	// cout<<"counter value was: "<<ctrVal<<endl;
	// cout<<"ctr array result was: "<<data[0][ctrInd]<<endl;
	// cout<<"value array result was: "<<data[1][ctrInd]<<endl;
	// return value from the value array
	return this->data[1][ctrInd];
}
// AS-2021

void GlobalFunction::printDetails(System *s)
{
	cout<<"Global Function: '"<< this->name << "()'"<<endl;
	cout<<" ="<<funcExpression<<endl;
	cout<<"   -Variable References:"<<endl;
	for(unsigned int vr=0; vr<n_varRefs; vr++) {
		cout<<"         "<<varRefTypes[vr]<<":  "<<varRefNames[vr]<<" = " << s->getObservableByName(varRefNames[vr])->getCount()<<endl;
	}
	cout<<"   -Constant Parameters:"<<endl;
	for(unsigned int i=0; i<n_params; i++) {
		cout<<"         "<<paramNames[i]<<" = " << s->getParameter(paramNames[i])<<endl;
	}

	
	if(p!=0) {
		// AS-2021
		if (this->fileFunc) {
			cout<<"   Function relies on file: "<<this->filePath<<endl;
			cout<<"   Function currently evaluates to: "<<this->fileEval()<<endl;
		} else {
			cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
		}
		// cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
		// AS-2021
	}
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

	//	cout<<"matched moleculeType"<<endl;
		value+=m->getComponentState(stateIndex);
	//	cout<<"found component state: "<<m->getComponentState(stateIndex)<<endl;
	//	cout<<"updating v`alue to: "<< value<<endl;
	}
}
