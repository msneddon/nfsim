



#include "reaction.hh"


using namespace std;
using namespace NFcore;




FunctionalRxnClass::FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1, transformationSet)
{
	this->cf=0;
	this->gf=gf;
	for(int vr=0; vr<gf->getNumOfVarRefs(); vr++) {
		if(gf->getVarRefType(vr)=="Observable") {
			Observable *obs = s->getObservableByName(gf->getVarRefName(vr));
			obs->addDependentRxn(this);
		} else {
			cerr<<"When creating a FunctionalRxnClass of name: "+name+" you provided a function that\n";
			cerr<<"depends on an observable type that I can't yet handle! (which is "+gf->getVarRefType(vr)+"\n";
			cerr<<"try using type: 'MoleculeObservable' for now.\n";
			cerr<<"quiting..."<<endl; exit(1);
		}
	}
}

FunctionalRxnClass::FunctionalRxnClass(string name, CompositeFunction *cf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1, transformationSet)
{
	this->gf=0;
	this->cf=cf;
	this->cf->setGlobalObservableDependency(this,s);
}


FunctionalRxnClass::~FunctionalRxnClass() {};

double FunctionalRxnClass::update_a() {
	//cout<<"udpating a"<<endl;
	if(this->onTheFlyObservables==false) {
		cerr<<"Warning!!  You have on the fly observables turned off, but you are using functional\n";
		cerr<<"reactions which depend on observables.  Therefore, you cannot turn off onTheFlyObservables!\n";
		cerr<<"exiting now."<<endl;
		exit(1);
	}



	//	cout<<"here"<<endl;
	if(gf!=0) {
	//	cout<<"in here"<<endl;
		a=FuncFactory::Eval(gf->p);
	} else if(cf!=0) {
		a=cf->evaluateOn(0,0);
	//	cout<<"and here"<<endl;
	} else {
		cout<<"Error!  Functional rxn is not properly initialized, but is being used!"<<endl;
		exit(1);
	}

	if(a<0) {
		cout<<"Warning!!  The function you provided for functional rxn: '"<<name<<"' evaluates\n";
		cout<<"to a value less than zero!  You cannot have a negative propensity!";
		cout<<"here is the offending function: \n";
		gf->printDetails();
		cout<<"\nhere is the offending reaction: \n";
		this->printDetails();
		cout<<"\n\nquitting."<<endl;
		exit(1);
	}

	//Finally, use a check to
	for(unsigned int i=0; i<n_reactants; i++) {
		if(reactantLists[i]->size()==0 && a>0.000000001) {
			a=0;
	//		cout<<"Warning!  Function evaluates to positive rate for a reaction, but"<<endl;
//			cout<<"one of the reactant lists is empty!"<<endl;
//			this->printDetails();
//			cf->printDetails(reactantTemplates[0]->getMoleculeType()->getSystem());
			//exit(1);
		}
	}

	return a;
}

void FunctionalRxnClass::printDetails() const {

	if(gf!=0)
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<gf->getNiceName()<<"="<<FuncFactory::Eval(gf->p)<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	else if(cf!=0)
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<cf->getName()<<"="<<cf->evaluateOn(0,0)<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









MMRxnClass::MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet ) :
	BasicRxnClass(name,1,transformationSet)
{
	this->Km = Km;
	this->kcat = kcat;
	this->sFree=0;
	if(n_reactants!=2) {
		cerr<<"You have tried to create a reaction with a Michaelis-Menten rate law (named: '"+name+"'\n')";
		cerr<<"but you don't have the correct number of reactants!  Michaelis-Menten reactions require\n";
		cerr<<"exactly 2 reactants.  A substrate (always given first) and an enzyme (always given second)\n";
		cerr<<"Read your tutorial next time... now I will quit."<<endl;
		exit(1);
	}
}
MMRxnClass::~MMRxnClass() {};

double MMRxnClass::update_a()
{
	double S = (double)reactantLists[0]->size();
	double E = (double)reactantLists[1]->size();
	sFree=0.5*( (S-Km-E) + pow((pow( (S-Km-E),2.0) + 4.0*Km*S),  0.5) );
	a=kcat*sFree*E/(Km+sFree);
	return a;
}

void MMRxnClass::printDetails() const {
	cout<<"ReactionClass: " << name <<"  ( Km="<<Km<<", kcat="<<kcat<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}









BasicRxnClass::BasicRxnClass(string name, double baseRate, TransformationSet *transformationSet) :
	ReactionClass(name,baseRate,transformationSet)
{
	this->reactionType = BASIC_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	reactantLists = new ReactantList *[n_reactants];
	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
}


BasicRxnClass::~BasicRxnClass()
{

//	this->reactantLists.at(0)->printDetails();

	if(DEBUG) cout<<"Destorying rxn: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
		//delete reactantTemplates[r]; DO NOT DELETE HERE (MoleculeType has responsibility of
		//deleting all template molecules of its type now.
		delete reactantLists[r];
	}
	delete [] reactantLists;
}

void BasicRxnClass::init()
{
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}


void BasicRxnClass::prepareForSimulation()
{

}






void BasicRxnClass::remove(Molecule *m, unsigned int reactantPos)
{

	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}


	//Get the specified reactantList
	ReactantList *rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);


	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}







void BasicRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex)
{
	cerr<<"You are trying to notify a Basic Reaction of a rate Factor Change!!! You should only use this"<<endl;
	cerr<<"function for DORrxnClass rules!  For this offense, I must abort now."<<endl;
	exit(1);
}
unsigned int BasicRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return reactantLists[reactantIndex]->size();
}
void BasicRxnClass::printFullDetails() const
{
	cout<<"BasicRxnClass: "<<name<<endl;
	for(unsigned int i=0; i<n_reactants; i++)
		reactantLists[i]->printDetails();
}
void BasicRxnClass::pickMappingSets(double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//So, we shall loop through the lists and extract out the MappingSets and
	//the molecule ID that the mappingSet was generated from
	//unsigned int *moleculeIDs = new unsigned int [n_reactants];
	for(unsigned int i=0; i<n_reactants; i++)
	{
		reactantLists[i]->pickRandom(mappingSet[i]);
		//mappingSet[i]->get(0)
	}

	//delete [] moleculeIDs;
}
































