#include "reaction.hh"


using namespace std;
using namespace NFcore;


FunctionalRxnClass::FunctionalRxnClass(string name, GlobalFunction *gf, TransformationSet *transformationSet, System *s) :
	BasicRxnClass(name,1,"",transformationSet,s)
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
	BasicRxnClass(name,1, "", transformationSet,s)
{
	this->gf=0;
	this->cf=cf;
	this->cf->setGlobalObservableDependency(this,s);
}


FunctionalRxnClass::~FunctionalRxnClass() {};


double FunctionalRxnClass::update_a() {
	if(this->onTheFlyObservables==false) {
		cerr<<"Warning!!  You have on the fly observables turned off, but you are using functional\n";
		cerr<<"reactions which depend on observables.  Therefore, you cannot turn off onTheFlyObservables!\n";
		cerr<<"exiting now."<<endl;
		exit(1);
	}

	if(gf!=0) {
		a=FuncFactory::Eval(gf->p);
	} else if(cf!=0) {
		int * reactantCounts = new int[this->n_reactants];
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCounts[r] = (int)getReactantCount(r);
		}
		a=cf->evaluateOn(0,0, reactantCounts, n_reactants);
		delete [] reactantCounts;
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

	// check here for the total rate flag - if this is set to true, then
	// use the rate exactly as given by the function, but if it is false,
	// then we have to multiply here by the reactant counts
	if(!this->totalRateFlag) {
		for(unsigned int i=0; i<n_reactants; i++)
			a*=(double)getCorrectedReactantCount(i);
	}
	else
	{
		// Check that we have at least one set of reactants!
		for(unsigned int i=0; i<n_reactants; i++) {
			if(getCorrectedReactantCount(i)==0) {
				a=0.0;
				break;
			}
		}
	}
	
	return a;
}


void FunctionalRxnClass::printDetails() const {

	string trate = "off";
	if(this->totalRateFlag) trate = "on";

	if(gf!=0)
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<gf->getNiceName()<<"="<<FuncFactory::Eval(gf->p)<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;
	else if(cf!=0) {
		int * reactantCounts = new int[this->n_reactants];
		for(unsigned int r=0; r<n_reactants; r++) {
			reactantCounts[r]=getReactantCount(r);
		}
		double value=cf->evaluateOn(0,0, reactantCounts, n_reactants);
		delete [] reactantCounts;
		cout<<"ReactionClass: " << name <<"  ( baseFunction="<<cf->getName()<<"="<<value<<",  a="<<a<<", fired="<<fireCounter<<" times, TotalRate="<<trate<<" )"<<endl;

	}

	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
		cout<<"	(count="<< this->getReactantCount(r) <<")."<<endl;
	}
	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}


MMRxnClass::MMRxnClass(string name, double kcat, double Km, TransformationSet *transformationSet,System *s) :
	BasicRxnClass(name,1,"",transformationSet,s)
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
	double S = (double)getCorrectedReactantCount(0);
	double E = (double)getCorrectedReactantCount(1);
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


BasicRxnClass::BasicRxnClass(string name, double baseRate, string baseRateName, TransformationSet *transformationSet, System *s) :
	ReactionClass(name,baseRate,baseRateName,transformationSet,s)
{
	this->reactionType = BASIC_RXN;  //set as normal reaction here, but deriving reaction classes can change this
	reactantLists = new ReactantList *[n_reactants];
	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
}


BasicRxnClass::~BasicRxnClass()
{
	if(DEBUG) cout<<"Destroying rxn: "<<name<<endl;

	for(unsigned int r=0; r<n_reactants; r++)
	{
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


bool BasicRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos)
{
	//Get the specified reactantList
	rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);

	//If this reaction has multiple instances, we always remove them all!
	// then we remap because other mappings may have changed.  Yes, this may
	// be more ineffecient, but it is the fast implementation
	if(rl->getHasClonedMappings()) {
		while(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
		}
	}

	//Here we get the standard update...
	while(m->getRxnListMappingId(rxnIndex)>=0) {
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
		//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}

	//Try to map it!
	ms = rl->pushNextAvailableMappingSet();
	symmetricMappingSet.clear();
	comparisonResult = reactantTemplates[reactantPos]->compare(m,rl,ms,false,&symmetricMappingSet);
	if(!comparisonResult) {
		//we must remove, if we did not match.  This will also remove
		//everything that was cloned off of the mapping set
		rl->removeMappingSet(ms->getId());
	} else {
		//TODO: it is necessary to remove elements that are not used anymore from the rl as well as from the m
		//for that
		//m->setRxnListMappingId(rxnIndex,-1);

		if (symmetricMappingSet.size() > 0){
            rl->removeMappingSet(ms->getId());
			for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
					m->setRxnListMappingId(rxnIndex,(*it)->getId());
            }
		}
		else{
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}

	return true;
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


double BasicRxnClass::update_a()
{
	// Use the total rate law convention (macroscopic rate)
	if(this->totalRateFlag) {
		a=baseRate;
		for(unsigned int i=0; i<n_reactants; i++)
			if(getCorrectedReactantCount(i)==0) a=0.0;

	// Use the standard microscopic rate
	} else {
		a = 1.0;
		for(unsigned int i=0; i<n_reactants; i++) {
			a*=getCorrectedReactantCount(i);
		}
		a*=baseRate;
	}
	return a;
}


int BasicRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	return isPopulationType[reactantIndex] ?
			   reactantLists[reactantIndex]->getPopulation()
			 : reactantLists[reactantIndex]->size();
}


int BasicRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}


void BasicRxnClass::pickMappingSets(double random_A_number) const
{
	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire

	//Select a reactant from each list
	for(unsigned int i=0; i<n_reactants; i++)
	{
		if ( isPopulationType[i] ) {
			reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
		} else {
			reactantLists[i]->pickRandom(mappingSet[i]);
		}
	}
}
