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
}



GlobalFunction::~GlobalFunction() 
{
	delete [] argNames;
	delete [] argTypes;
	delete [] paramNames;
	delete [] paramValues;
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


void GlobalFunction::printDetails()
{
	cout<<"--------"<<endl;
	cout<<"Details of Function: '"<< this->name << "'"<<endl;
	cout<<" ="<<funcString<<endl;
	cout<<"   -Arguements:"<<endl;
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
	
