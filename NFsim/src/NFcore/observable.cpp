#include <iostream>
#include "observable.hh"


using namespace std;
using namespace NFcore;





Observable::Observable(string name)
{
	this->obsName = name;
	this->n_templates= 0;
	this->templateMolecules = 0;
	this->n_dependentRxns=0;
	this->dependentRxns= new ReactionClass *[n_dependentRxns];
	this->count=0;
	this->type=Observable::NO_TYPE;
}

Observable::~Observable()
{
	delete [] dependentRxns;
	delete [] templateMolecules;

	this->n_templates = 0;
	this->templateMolecules = 0;
	this->n_dependentRxns = 0;
	this->dependentRxns = 0;
}

void Observable::add()
{
	//First, we add to our observable count
	count++;

	//Next, we update our dependent reactions, if there are any
	for(int r=0; r<n_dependentRxns; r++) {
		double old_a = dependentRxns[r]->get_a();
		double new_a = dependentRxns[r]->update_a();
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot(old_a,new_a);
	}
}


void Observable::straightAdd()
{
	count++;
}

void Observable::straightSubtract()
{
	count--;
}

void Observable::subtract()
{
	if(count==0){
		cout<<"Error in observable count!! Removing from an empty observable!"<<endl;
		cout<<"Observable named: "<<obsName<<endl;
		exit(1);
	}

	count--;

	//Next, we update our dependent reactions, if there are any
	for(int r=0; r<n_dependentRxns; r++) {
		double old_a = dependentRxns[r]->get_a();
		double new_a = dependentRxns[r]->update_a();
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot(old_a,new_a);
	}
}


void Observable::getTemplateMoleculeList(int &n_templates, TemplateMolecule **&tmList)
{
	n_templates = this->n_templates;
	tmList = this->templateMolecules;
}


void Observable::addReferenceToMyself(mu::Parser *p)
{
	p->DefineVar(obsName,&count);
}
void Observable::addReferenceToMyself(string referenceName, mu::Parser *p)
{
	p->DefineVar(referenceName,&count);
}
void Observable::addDependentRxn(ReactionClass *r)
{
	cout<<"Observable: "<<this->obsName<<" adding dependent rxn: "<<r->getName()<<endl;
	cout<<"n dependent rxns: "<<n_dependentRxns<<endl;
	ReactionClass ** newDepRxns = new ReactionClass * [n_dependentRxns+1];
	for(int i=0; i<n_dependentRxns; i++) {
		newDepRxns[i] = dependentRxns[i];
	}
	newDepRxns[n_dependentRxns] = r;
	delete [] dependentRxns;
	dependentRxns = newDepRxns;
	n_dependentRxns++;
}




///////////////////////////////////////////////////////////////////////////////////////////////



MoleculesObservable::MoleculesObservable(string name, TemplateMolecule *tm) :
	Observable(name)
{
	n_templates=1;
	templateMolecules = new TemplateMolecule * [n_templates];
	templateMolecules[0]=tm;

	this->type=Observable::MOLECULES;
}
MoleculesObservable::MoleculesObservable(string name, vector <TemplateMolecule *> &tmList) :
	Observable(name)
{
	n_templates=tmList.size();
	templateMolecules = new TemplateMolecule * [n_templates];
	for(int t=0; t<n_templates; t++) {
		templateMolecules[t] = tmList.at(t);
	}

	this->type=Observable::MOLECULES;
}

MoleculesObservable::~MoleculesObservable()
{

}


//NOTE:  cloning a molecule DOES NOT clone any dependent reactions.  You will have to
//add that again to this observable, if there are any.  This is the behaviour, because
//observables are generally only cloned in this way for local functions
Observable * MoleculesObservable::clone() {
	vector <TemplateMolecule *> tmList;
	for(int t=0; t<n_templates; t++)
		tmList.push_back(templateMolecules[t]);
	return new MoleculesObservable(obsName+"_clone",tmList);
}

int MoleculesObservable::isObservable(Molecule *m) const
{
	int matches = 0;
	for(int t=0; t<n_templates; t++) {
		if(templateMolecules[t]->compare(m)) matches++;
	}
	return matches;
}

int MoleculesObservable::isObservable(Complex *c) const
{
	cerr<<"Comparing a Molecules observable '"<<obsName<<"' to a complex!"<<endl;
	cerr<<"You can only compare Species observable to a complexes!  Quitting."<<endl;
	exit(1);
}




///////////////////////////////////////////////////////////////////


SpeciesObservable::SpeciesObservable(string name, vector <TemplateMolecule *> &tmList) :
	Observable(name)
{
	n_templates=tmList.size();
	templateMolecules = new TemplateMolecule * [n_templates];
	for(int t=0; t<n_templates; t++) {
		templateMolecules[t] = tmList.at(t);
	}

	this->type=Observable::SPECIES;
}

SpeciesObservable::~SpeciesObservable()
{

}


//NOTE:  cloning a molecule DOES NOT clone any dependent reactions.  You will have to
//add that again to this observable, if there are any.  This is the behaviour, because
//observables are generally only cloned in this way for local functions
Observable * SpeciesObservable::clone() {
	vector <TemplateMolecule *> tmList;
	for(int t=0; t<n_templates; t++)
		tmList.push_back(templateMolecules[t]);
	return new SpeciesObservable(obsName+"_clone",tmList);
}

int SpeciesObservable::isObservable(Molecule *m) const
{
	cerr<<"Comparing a Species observable '"<<obsName<<"' to a molecule!"<<endl;
	cerr<<"You can only compare Species observable to a complexes!  Quitting."<<endl;
	exit(1);

}

int SpeciesObservable::isObservable(Complex *c) const
{
	int matches = 0;
	for(int t=0; t<n_templates; t++) {
		for(c->molIter=c->complexMembers.begin(); c->molIter!=c->complexMembers.end();c->molIter++) {

			//For each template, we only have to find one match, then we match for sure.
			if(templateMolecules[t]->compare((*c->molIter))) {
				matches++;
				break;
			}
		}
	}

	return matches;
}





