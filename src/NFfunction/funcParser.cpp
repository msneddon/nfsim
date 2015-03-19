#include "NFfunction.hh"

#include <math.h>

using namespace NFcore;
using namespace mu;




Parser * FuncFactory::create(string functionString, vector <string> & variableNames, vector <double *> & variablePtrs)
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 6.02214179e23;

	Parser *p = new Parser();
	try
	{
		p->DefineConst("_PI",PI);
		p->DefineConst("_e",E);
		p->DefineConst("_Na",NA);

		if(variableNames.size()!=variablePtrs.size())  {
			cout<<"Error parsing function in FuncFactory!!  Your variableNames vector and "<<endl;
			cout<<"you variablePtrs vector do not appear to match up!"<<endl;
			cout<<"The function you gave me was: "<<functionString<<endl;
			cout<<"For that, I will quit"<<endl;
			exit(1);
		} else {
			for(unsigned int v=0;v<variableNames.size(); v++) {
				p->DefineVar(variableNames.at(v), variablePtrs.at(v));
			}
		}

		p->SetExpr(functionString);
	}
	catch (Parser::exception_type &e)
	{
		cout<<"Error parsing function in FuncFactory!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}

	return p;
}

Parser * FuncFactory::create()
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 6.02214179e23;
	Parser *p = new Parser();
	try
	{
		p->DefineConst("_PI",PI);
		p->DefineConst("_e",E);
		p->DefineConst("_Na",NA);
	}
	catch (Parser::exception_type &e)
	{
		cout<<"Error creating function in FuncFactory!!  This is what happened:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Quitting."<<endl;
		exit(1);
	}
	return p;
}




double FuncFactory::Eval(mu::Parser *p)
{
	if(p==NULL) {
		cout<<"In FuncFactory: Trying to evaluate a null Parser! You are probably trying\n";
		cout<<"to use a GlobalFunction before it has been prepared! Preparing a GlobalFunction\n";
		cout<<"connects it to Observables, so it must be done before you can use it!\n";
		cout<<"  we've all made this mistake before, but now I'm exiting..."<<endl;
		exit(1);
	}
	try {
		return p->Eval();
	} catch (Parser::exception_type &e) {
		cout<<"Error evaluating function in FuncFactory!!  "<<endl;
		cout<<"The function was: "<<p->GetExpr()<<endl;
		cout<<"And this is what went wrong:"<<endl;
		cout<< "  "<<e.GetMsg() << endl;
		cout<<"Terminating your simulation. Better luck next time."<<endl;
		exit(1);
	}
	return 0;
}


void FuncFactory::test()
{
	double PI = 3.14159265358979323846;
	double NA= 6.02214179e23;
	double E = 6.02214179e23;
	cout<<"Beginning diagnostic tests..."<<endl;

	{
	//Test 1: check constants are correct and basic math functions seem right
	cout<<" 1) simple test of constants and predefined functions: ";
	string functionString("sin(_e*cos(3.2/_PI))+ln(_Na*1.14e-11)");//sin(_e*cos(3.2/_PI)+0.11*(1-_Na*1.14)");
	vector <string> variableNames;
	vector <double *> variablePtrs;

	Parser *p = FuncFactory::create(functionString,variableNames,variablePtrs);
	double result = sin(E*cos(3.2/PI))+log(NA*1.14e-11);
	double funcResult = p->Eval();
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;
	delete p;
	}


	{
	//Test 2a: check that variables are working properly
	cout<<" 2a) test that variable input works properly: ";

	string functionString = "1-(d1/d2)*sin(d1*d2)+1.3*d1^2";

	string d1_name("d1");
	string d2_name("d2");
	vector <string> variableNames;
	variableNames.push_back(d1_name);
	variableNames.push_back(d2_name);

	double d1 = -1.412;
	double d2 = 2.01e1;
	vector <double *> variablePtrs;
	variablePtrs.push_back(&d1);
	variablePtrs.push_back(&d2);

	Parser *p = FuncFactory::create(functionString,variableNames,variablePtrs);
	double result = 1-(d1/d2)*sin(d1*d2)+1.3*pow(d1,2.0);
	double funcResult = FuncFactory::Eval(p);
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;


	//Test 2b: check that variables are working properly
	cout<<" 2b) another check of the variables: ";
	d1 = 50.1;
	d2 = 3.1;
	result = 1-(d1/d2)*sin(d1*d2)+1.3*pow(d1,2.0);
	funcResult = FuncFactory::Eval(p);
	if(abs(funcResult - result)<0.0001)
		cout<<"pass."<<endl;
	else
		cout<<"fail! p->Eval() = "<<funcResult<<"  but should be: "<<result<<endl;

	delete p;
	}

	//Thats all the test I can think of!
	cout<<endl<<"Testing complete."<<endl;
}

