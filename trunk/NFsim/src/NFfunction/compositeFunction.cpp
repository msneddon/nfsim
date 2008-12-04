/*
 * compositeFunction.cpp
 *
 *  Created on: Dec 4, 2008
 *      Author: msneddon
 */

#include "NFfunction.hh"



using namespace std;
using namespace NFcore;
using namespace mu;


CompositeFunction::CompositeFunction(System *s,
					string name,
					string expression,
					vector <string> &functions,
					vector <string> &paramNames)
{

	cout<<"creating composite function"<<endl;
	this->name = name;
	this->originalExpression=expression;
	this->parsedExpression="";


	this->n_params=paramNames.size();
	this->paramNames = new string[n_params];
	for(unsigned int i=0; i<n_params; i++) {
		this->paramNames[i]=paramNames.at(i);
	}


	this->n_allFuncs=functions.size();
	allFuncNames = new string[n_allFuncs];
	for(unsigned int i=0; i<n_allFuncs; i++) {
		this->allFuncNames[i] = functions.at(i);
	}

	p=0;
}
CompositeFunction::~CompositeFunction()
{


}


void CompositeFunction::setGlobalObservableDependency(ReactionClass *r, System *s) {

	for(int i=0; i<this->n_gfs; i++) {
		GlobalFunction *gf=gfs[i];

		for(int vr=0; vr<gf->getNumOfVarRefs(); vr++) {
			if(gf->getVarRefType(vr)=="Observable") {
				Observable *obs = s->getObservableByName(gf->getVarRefName(vr));
				obs->addDependentRxn(r);
			} else {
				cerr<<"When creating a FunctionalRxnClass of name: "+r->getName()+" you provided a function that\n";
				cerr<<"depends on an observable type that I can't yet handle! (which is "+gf->getVarRefType(vr)+"\n";
				cerr<<"try using type: 'MoleculeObservable' for now.\n";
				cerr<<"quiting..."<<endl; exit(1);
			}
		}
	}
}

//call this immediately after you have read in all the functions, but before preparing for the simulation
//or reading in reactions
void CompositeFunction::finalizeInitialization(System *s)
{
	//first, find the global functions by name
	vector <GlobalFunction *> gf_tempVector;
	for(unsigned int i=0; i<n_allFuncs; i++) {
		GlobalFunction *gf=s->getGlobalFunctionByName(allFuncNames[i]);
		if(gf!=0) gf_tempVector.push_back(gf);
	}

	this->n_gfs = gf_tempVector.size();
	this->gfNames = new string[n_gfs];
	this->gfValues = new double[n_gfs];
	this->gfs = new GlobalFunction * [n_gfs];
	for(int i=0; i<n_gfs; i++) {
		gfNames[i] = gf_tempVector.at(i)->getName();
		gfValues[i] = 0;
		gfs[i] = gf_tempVector.at(i);
	}

	//now the local functions...
	vector <LocalFunction *> lf_tempVector;
	for(unsigned int i=0; i<n_allFuncs; i++) {
		LocalFunction *lf=s->getLocalFunctionByName(allFuncNames[i]);
		if(lf!=0) lf_tempVector.push_back(lf);
	}

	this->n_lfs = lf_tempVector.size();
	this->lfNames = new string[n_lfs];
	this->lfValues = new double[n_lfs];
	this->lfs = new LocalFunction * [n_lfs];
	for(int i=0; i<n_lfs; i++) {
		lfNames[i] = lf_tempVector.at(i)->getName();
		lfValues[i] = 0;
		lfs[i] = lf_tempVector.at(i);
	}



	//Now we have to do the dirty work of parsing out the function expression
	//so that we can properly get at the variables and functions we need
	this->parsedExpression=originalExpression;

	//First
	for(int f=0; f<n_gfs; f++) {
		string::size_type sPos=parsedExpression.find(gfNames[f]);
		for( ; sPos!=string::npos; sPos=parsedExpression.find(gfNames[f],sPos+1)) {

			string::size_type openPar = parsedExpression.find_first_of('(',sPos);
			string::size_type closePar = parsedExpression.find_first_of(')',sPos);
			if(openPar!=string::npos && closePar!=string::npos) {
				if(closePar>openPar) { //if we got here, we found a valid parenthesis to look at
					string inBetween = parsedExpression.substr(openPar+1,closePar-openPar-1);
					NFutil::trim(inBetween);

					if(inBetween.size()==0) {
						parsedExpression.replace(openPar,closePar-openPar+1,"");
					}
				}
			}
		}
	}

	cout<<"now the expression is: "<<parsedExpression<<endl;


	///////// do the same for local functions here (can be a bit tricky, because different
	///////// composite functions will have different numbers of arguments



}

void CompositeFunction::updateParameters(System *s)
{




}

void CompositeFunction::prepareForSimulation(System *s)
{
	try {
			p=FuncFactory::create();
			for(int f=0; f<n_gfs; f++) {
				p->DefineVar(gfNames[f],&gfValues[f]);
			}

			//Define local function variables here...

			for(unsigned int i=0; i<n_params; i++) {
				p->DefineConst(paramNames[i],s->getParameter(paramNames[i]));
			}
			p->SetExpr(this->parsedExpression);
		}
		catch (mu::Parser::exception_type &e)
		{
			cout<<"Error preparing function "<<name<<" in class CompositeFunction!!  This is what happened:"<<endl;
			cout<< "  "<<e.GetMsg() << endl;
			cout<<"Quitting."<<endl;
			exit(1);
		}


	cout<<"preparing composite function.."<<endl;
//	exit(0);
}


double CompositeFunction::evaluateOn(Molecule **molList) {

	//1 evaluate all global functions
	for(int f=0; f<n_gfs; f++) {
		gfValues[f]=FuncFactory::Eval(gfs[f]->p);
	}


	//2 evaluate all local functions

	//evaluate this function
	return FuncFactory::Eval(p);
}
