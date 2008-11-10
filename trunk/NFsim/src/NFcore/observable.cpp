#include <iostream>
#include "NFcore.hh"


using namespace std;
using namespace NFcore;


Observable::Observable(string aliasName, TemplateMolecule * templateMolecule)
{
	this->aliasName = aliasName;
	this->templateMolecule = templateMolecule;
	this->count = 0;
}

Observable::~Observable()
{
	templateMolecule = 0;
}


void Observable::add() {
	//cout<<"adding to observable: "<<this->aliasName<<endl;
	count++;
	for(rxnIter = dependentRxns.begin(); rxnIter != dependentRxns.end(); rxnIter++ ) {
		double old_a = (*rxnIter)->get_a();
		templateMolecule->getMoleculeType()->getSystem()->update_A_tot(old_a,(*rxnIter)->update_a());
	}
};


void Observable::subtract() {
	if(count==0){ cout<<"Error in observable count!!"<<endl; exit(1); }
	count--;
	//cout<<"subtracting from observable: "<<this->aliasName<<endl;
	for(rxnIter = dependentRxns.begin(); rxnIter != dependentRxns.end(); rxnIter++ ) {
		double old_a = (*rxnIter)->get_a();
		templateMolecule->getMoleculeType()->getSystem()->update_A_tot(old_a,(*rxnIter)->update_a());
	}
};

void Observable::addDependentRxn(ReactionClass *r) {
	this->dependentRxns.push_back(r);
}



bool Observable::isObservable(Molecule * m) const {

	//cout<<"comparing obs: "<<this->aliasName<<" to: "<<endl;
	//m->printDetails();
	//templateMolecule->printDetails();
	bool answer = templateMolecule->compare(m);
	//if(answer) cout<<"Match."<<endl;
	//else cout<<"nope"<<endl;

	return answer;


}


void Observable::addReferenceToMyself(mu::Parser * p) {
	p->DefineVar(this->aliasName,&count);
}

