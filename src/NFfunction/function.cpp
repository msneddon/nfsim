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

	// AS-2021
	this->fileFunc = false;
	// AS-2021
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
		if (this->fileFunc==true) {
			this->fileUpdate();
		}
		// AS-2021
		cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
	}
}

// AS-2021
void GlobalFunction::loadParamFile(string filePath) 
{
	// setup our vectors
	vector <double> time;
	vector <double> values;
	// open file for reading
	ifstream file(filePath.c_str());
	// Report if file doesn't exist
	if(!file.good()){
		cout<<"Error preparing function "<<this->name<<" in class GlobalFunction!!"<<endl;
		cout<<"File doesn't look like it exists"<<endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}
	// TODO: Err out this doesn't work
	try {
		// strings for looping over the file
		string line, word, content;
		string a,b;
		// TODO: Err out if the format is wrong
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
	} catch (exception const & e) {
		cout<<"Error preparing function "<<this->name<<" in class GlobalFunction!!"<<endl;
		cout<<"Failed to either open or read the file."<<endl;
		cout<<"Quitting."<<endl;
		exit(1);
	};
	return;
};

void GlobalFunction::addCounterPointer(double *counter){
	this->ctrType = "Observable";
	this->counter = counter;
}

void GlobalFunction::setCtrName(string name) {
	this->ctrName = name;
}

// unhooking system timer option for now
// void GlobalFunction::addSystemPointer(System *s) {
// 	this->ctrType = "System";
// 	this->sysPtr = s;
// }

void GlobalFunction::enableFileDependency(string filePath) {
	// load file
	// TODO: Err out if this fails
	try {
		this->loadParamFile(filePath);
	} catch (exception const & e) {
		cout<<"Error preparing function "<<name<<" in class GlobalFunction!!"<<endl;
		cout<<"Quitting."<<endl;
		exit(1);
	};
	// we just want to keep a record of this
	this->filePath = filePath;
	// this sets it up so that this function knows it's supposed
	// to be pulling values from a file
	this->fileFunc = true;
	// initialize internal index
	this->currInd = 0;
	// pull data lenght so we can reuse it
	this->dataLen = data[0].size();
}

double GlobalFunction::getCounterValue() {
	// depending on the type of the observable counter
	// get the actual value
	double ctrVal;
	if (ctrType == "Observable") {
		ctrVal = (*counter);
	} 
	// unhooking system timer option for now
	// else {
	// 	// not sure but this is likely slower
	// 	ctrVal = this->sysPtr->getCurrentTime();
	// }
	return ctrVal;
}
void GlobalFunction::fileUpdate() {
	// TODO: Error checking and reporting
	
	// get counter val
	double ctrVal = this->getCounterValue();

	// basic step function implementation
	// if we got past the last point, keep returning
	// the last point
	if (currInd>dataLen-1) {
		currInd = dataLen-1;
		p->DefineConst(ctrName,data[1][currInd]);
		return;
	} else if (currInd==dataLen-1) {
		p->DefineConst(ctrName,data[1][currInd]);
		return;
	}
	// a simple way to do interval locating 
	if (data[0][currInd] < data[0][currInd+1]) {
		// next point is higher than the current point, we
		// are waiting for the counter value to be higher 
		// than our current point
		
		// return 0 if we don't have data yet
		if(data[0][0]>=ctrVal) {
			// we haven't gotten to the point where
			// we can get a value out, return 0
			// cout<<"not there yet, returning 0"<<endl;
			p->DefineConst(ctrName,0);
			return;
		} 
		// go up by one if the counter value got past 
		// the next value in the array
		if (ctrVal>=data[0][currInd+1]) {
			currInd += 1;
		}
	// note that this makes no sense if they are equal
	// TODO: Raise error if they are equal. Better yet, parse 
	// it ahead of time and make sure that doesn't happen
	} else {
		// next point is lower than the current point, we
		// are waiting for the counter value to be lower 
		// than our current point

		// return 0 if we don't have data yet
		if(data[0][0]<=ctrVal) {
			// we haven't gotten to the point where
			// we can get a value out, return 0
			// cout<<"not there yet, returning 0"<<endl;
			p->DefineConst(ctrName,0);
			return;
		}
		// go up by one if the counter value got past 
		// the next value in the array
		if (ctrVal<=data[0][currInd+1]) {
			currInd += 1;
		}
	}
	// // return value from the value array
	p->DefineConst(ctrName,data[1][currInd]);
	return;
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
		if (this->fileFunc==true) {
			cout<<"   Function relies on file: "<<this->filePath<<endl;
			this->fileUpdate();
		}
		// AS-2021
		cout<<"   Function currently evaluates to: "<<FuncFactory::Eval(p)<<endl;
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
