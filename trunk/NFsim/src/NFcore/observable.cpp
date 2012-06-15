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
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot(dependentRxns[r],old_a,new_a);
	}
}

/* add multiple new matches to an observable (rather than call 'add' a bunch of times --justin */
/* useful for counters! */
void Observable::add( int n_matches )
{
	//debug
	//cout << "Observable::add( " << n_matches << " )" << endl;

	// First, we add to our observable count
	count += n_matches;

	// Next, we update our dependent reactions, if there are any
	for (int r=0; r<n_dependentRxns; r++)
	{
		double old_a = dependentRxns[r]->get_a();
		double new_a = dependentRxns[r]->update_a();
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot( dependentRxns[r], old_a, new_a);
	}
}


void Observable::straightAdd()
{
	count++;
}

void Observable::straightAdd(int n_matches)
{
	count += n_matches;
}

void Observable::subtract()
{
	if(count==0){
		cerr << "Error in observable count!! Removing from an empty observable!"
		     << "Observable named: " << obsName << endl;
		exit(1);
	}

	count--;

	//Next, we update our dependent reactions, if there are any
	for(int r=0; r<n_dependentRxns; r++) {
		double old_a = dependentRxns[r]->get_a();
		double new_a = dependentRxns[r]->update_a();
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot(dependentRxns[r],old_a,new_a);
	}
}

/* Remove multiple matches fron an observable (rather than call 'subtract' a bunch of times --justin */
/* Necessary to make counters fast! */
void Observable::subtract( int n_matches )
{
	//debug
	//cout << "Observable::subtract( " << n_matches << " )" << endl;

	if (count - n_matches < 0)
	{
		cerr << "Error in observable count!! Removing " << n_matches << " matches will result in a negative match count!"
		     << "Observable named: " << obsName << endl;
		exit(1);
	}

	// First, we subtract from our observable count
	count -= n_matches;

	// Next, we update our dependent reactions, if there are any
	for (int r=0; r<n_dependentRxns; r++)
	{
		double old_a = dependentRxns[r]->get_a();
		double new_a = dependentRxns[r]->update_a();
		templateMolecules[0]->getMoleculeType()->getSystem()->update_A_tot( dependentRxns[r], old_a, new_a);
	}
}

void Observable::straightSubtract()
{
	count--;
}

void Observable::straightSubtract(int n_matches)
{
	count -= n_matches;
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
	//cout<<"Observable: "<<this->obsName<<" adding dependent rxn: "<<r->getName()<<endl;
	//cout<<"n dependent rxns: "<<n_dependentRxns<<endl;
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

//	cout<<" creating observable "<< name <<endl;
//	if(this->getName()=="CbpTot")
//	{
//		vector <TemplateMolecule *> tms;
//		TemplateMolecule::traverse(tmList.at(0),tms,false);
//		for(int k=0; k<tms.size(); k++)
//		{
//			tms.at(k)->printDetails(cout);
//		}
//		exit(1);
//	}

	//tmList.at(0)->printDetails(cout);
	//tmList.at(0)->printPattern(cout);
	//cout<<"-------------\n"<<endl;

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
	// DEBUG code
	//cout << "MoleculesObservable::isObservable( " << m->getUniqueID() << " )" << endl;
	//m->printDetails();
	//cout << "observable name: " << obsName << endl;
	//cout << "number of templates: " << n_templates << endl;

	int matches = 0;
	for(int t=0; t<n_templates; t++) {
		//cout << "moleculeType: " << (templateMolecules[t]->getMoleculeTypeName()) << endl;
		//cout << "connected to: " << (templateMolecules[t]->getN_connectedTo()) << endl;
		//templateMolecules[t]->printDetails();
		// try to get match counts, rather than just a boolean
		//cout<<endl<<endl<<endl;
		//cout<<"starting!"<<endl;

		if ( templateMolecules[t]->compare(m) ) {
			//cout<<"  adding one"<<endl;
			matches += m->getPopulation();
			//return 1;
		}
		//else { cout<<"  nothing."<<endl; }
	}
	//cout << "total_matches: " << matches << endl;
	return matches;
}

int MoleculesObservable::isObservable(Complex *c) const
{
	cerr<<"Comparing a Molecules observable '"<<obsName<<"' to a complex!"<<endl;
	cerr<<"You can only compare Species observable to a complexes!  Quitting."<<endl;
	exit(1);
}




///////////////////////////////////////////////////////////////////


SpeciesObservable::SpeciesObservable(string name, vector <TemplateMolecule *> &tmList, vector <string> &stochRelation, vector <int> &stochQuantity) :
	Observable(name)
{
	n_templates=tmList.size();
	templateMolecules = new TemplateMolecule * [n_templates];
	for(int t=0; t<n_templates; t++) {
		templateMolecules[t] = tmList.at(t);
	}


	relation = new int [n_templates];
	quantity = new int [n_templates];
	for(int t=0; t<n_templates; t++) {
		if(stochRelation.at(t).empty()) {
			relation[t]=NO_RELATION;
			quantity[t]=-1;

		} else if (stochRelation.at(t).compare("==")==0) {
			relation[t]=EQUALS;
			quantity[t]=stochQuantity.at(t);

		} else if (stochRelation.at(t).compare("!=")==0) {
			relation[t]=NOT_EQUALS;
			quantity[t]=stochQuantity.at(t);

		} else if (stochRelation.at(t).compare(">")==0) {
			relation[t]=GREATER_THAN;
			quantity[t]=stochQuantity.at(t);

		} else if (stochRelation.at(t).compare("<")==0) {
			relation[t]=LESS_THAN;
			quantity[t]=stochQuantity.at(t);

		} else if (stochRelation.at(t).compare(">=")==0) {
			relation[t]=GREATOR_OR_EQUAL_TO;
			quantity[t]=stochQuantity.at(t);

		} else if (stochRelation.at(t).compare("<=")==0) {
			relation[t]=LESS_THAN_OR_EQUAL_TO;
			quantity[t]=stochQuantity.at(t);

		} else {
			relation[t]=NO_RELATION;
			quantity[t]=-1;
		}
	}

	this->type=Observable::SPECIES;
}

SpeciesObservable::~SpeciesObservable()
{
	delete [] relation;
	delete [] quantity;
}


//NOTE:  cloning a molecule DOES NOT clone any dependent reactions.  You will have to
//add that again to this observable, if there are any.  This is the behaviour, because
//observables are generally only cloned in this way for local functions
Observable * SpeciesObservable::clone() {
	vector <TemplateMolecule *> tmList;
	cout<<"in clone species observable, this is not yet updated to handle stoch observables.  fix me."<<endl;
	exit(1);
	for(int t=0; t<n_templates; t++)
		tmList.push_back(templateMolecules[t]);
	return 0; //new SpeciesObservable(obsName+"_clone",tmList);
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
		if(relation[t]==NO_RELATION) {
			for(c->molIter=c->complexMembers.begin(); c->molIter!=c->complexMembers.end();c->molIter++) {
				//For each template, we only have to find one match, then we match for sure.
				if ( templateMolecules[t]->compare(*(c->molIter)) ) {
					matches += (*(c->molIter))->getPopulation();
					break;
				}
			}
		}

		//Handle stoichiometric observables
		else {
			int localMatches = 0;
			for(c->molIter=c->complexMembers.begin(); c->molIter!=c->complexMembers.end();c->molIter++) {
				//For each template, we only have to find one match, then we match for sure.
				if(templateMolecules[t]->compare(*(c->molIter)) ) {
					localMatches++;
				}
			}

			// reset iterator to beginning
			c->molIter = c->complexMembers.begin();

			if(relation[t]==EQUALS) {
				if(localMatches==quantity[t]) matches += (*(c->molIter))->getPopulation();

			} else if(relation[t]==NOT_EQUALS) {
				if(localMatches!=quantity[t]) matches += (*(c->molIter))->getPopulation();

			} else if(relation[t]==GREATER_THAN) {
				if(localMatches>quantity[t])  matches += (*(c->molIter))->getPopulation();

			} else if(relation[t]==LESS_THAN) {
				if(localMatches<quantity[t])  matches += (*(c->molIter))->getPopulation();

			} else if(relation[t]==GREATOR_OR_EQUAL_TO) {
				if(localMatches>=quantity[t]) matches += (*(c->molIter))->getPopulation();

			} else if(relation[t]==LESS_THAN_OR_EQUAL_TO) {
				if(localMatches<=quantity[t]) matches += (*(c->molIter))->getPopulation();

			}

		}
	}
	return matches;
}





