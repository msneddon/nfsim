//
//
//
//
//
//
//
//ObsDrxn::ObsDrxn( char * name, 
//				int n_reactants, 
//				TemplateMolecule ** reactantTemplates,
//				int n_observables,
//				Observable ** obs) : ReactionClass(name, n_reactants, reactantTemplates, 0)
//{
//	if(DEBUG) cout<<"Creating generalized Observable Dependent Reaction..."<<endl;
//	this->obs = obs;
//}
//
//
//ObsDrxn::~ObsDrxn()
//{
//	delete [] obs;
//}
//
//double ObsDrxn::update_a()
//{
//	cerr<<" !!! Calling update_a from a generalized ObsDrxn!!!  ("<<name<<")"<<endl;
//	
//	a = 1;
//	for(int r=0; r<n_reactants; r++)
//		a*= (lastIndex[r]+1);
//		
//	a*=rate;
//	return a;
//}
//
//
//void ObsDrxn::transformReactants(Molecule ** reactants, int nReactants)
//{
//	cerr<<" !!! Calling transformReactants from a generalized ObsDrxn!!!  ("<<name<<")"<<endl;
//}
//
//
//void ObsDrxn::printDetails() const
//{
//	cout<<"** ObsDrxn: " << name <<"  ( rate="<<rate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
//	for(int r=0; r<n_reactants; r++)
//	{
//		cout<<"      -"<< this->reactantTemplates[r]->getMoleculeTypeName();
//		cout<<"	(count="<< lastIndex[r]+1 <<")."<<endl;
//	}
//	
//}

