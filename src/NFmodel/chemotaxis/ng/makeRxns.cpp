#include "ng.hh"


using namespace NG;




void NG::rxn_freeRB_bind_unbind_tether(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB)
{
	//Free R binding tether
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("tether");
	reactantTemplates[1] = new TemplateMolecule(cheR);
	reactantTemplates[1]->addEmptyBindingSite("te");
	reactantTemplates[1]->addEmptyBindingSite("av");
	ReactionClass * r = new ReactionSimpleBinding("FREE_R_bind_TETHER", 
		n_reactants, reactantTemplates, Param::FREE_CHER_bind_TETHER, "tether", "te");
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	//Free B binding tether
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("tether");
	reactantTemplates[1] = new TemplateMolecule(cheB);
	reactantTemplates[1]->addEmptyBindingSite("te");
	reactantTemplates[1]->addEmptyBindingSite("av");
	r = new ReactionSimpleBinding("FREE_B_bind_TETHER", 
		n_reactants, reactantTemplates, Param::FREE_CHEB_bind_TETHER, "tether", "te");
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	//R unbind the tether
	n_reactants = 1;
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheR);
	TemplateMolecule *receptor = new TemplateMolecule(recDimer);
	TemplateMolecule::bind(reactantTemplates[0], "te", receptor, "tether");
	r = new ReactionUnbinding("R_unbind_TETHER",
		n_reactants,reactantTemplates,Param::CHER_unbind_TETHER,"te");
	r->setTraversalLimit(3);
	s->addReaction(r);
	
	//B unbind the tether
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheB);
	receptor = new TemplateMolecule(recDimer);
	TemplateMolecule::bind(reactantTemplates[0], "te", receptor, "tether");
	r = new ReactionUnbinding("B_unbind_TETHER",
		n_reactants,reactantTemplates,Param::CHEB_unbind_TETHER,"te");
	r->setTraversalLimit(3);
	s->addReaction(r);
}

void NG::rxn_freeRB_bind_unbind_active(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB)
{
	
	//First, free binding of CheR to the active site depends on the number of active
	//sites available.  Thus, we need to use a DOR rxn to keep track of this
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
	/////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	reactantTemplates[0]->addNotStateValue("m",1);  // can't bind if we are fully methylated
	
	reactantTemplates[1] = new TemplateMolecule(cheR);
	reactantTemplates[1]->addEmptyBindingSite("av");  //has to have an empty site to bind the active site
	reactantTemplates[1]->addEmptyBindingSite("te"); //this is the free rxn, so tether must be free too
	
	//Only the receptor is the DOR reactant	
	int DORreactantIndex = 0;
	
	//Only tar needs the DOR group name
	char *DORgroupName = DIMER_GROUP_NAME;
	
	//We can only bind to an active site based on the number of active sites that are
	//available for us to bind (in this case, the free site level)
	int DORgroupValueIndex = DimerGroup::FREE_SITE_LEVEL;

	ReactionClass * r = new NG_FreeActiveSiteBinding("FREE_R_bind_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::FREE_CHER_bind_ACTIVE, "av", "asite");
	s->addReaction(r);
	r->setTraversalLimit(2);
	
	
	
	//Now do the same for free cheB binding to receptor
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
	reactantTemplates[0]->addNotStateValue("m",0);  // can't bind if we are fully demethylated
	reactantTemplates[1] = new TemplateMolecule(cheB);
	if(Param::useCheBFeedbackLoop) reactantTemplates[1]->addStateValue("p",PHOS);
	reactantTemplates[1]->addEmptyBindingSite("av");  //has to have an empty site to bind the active site
	reactantTemplates[1]->addEmptyBindingSite("te"); //this is the free rxn, so tether must be free too
	DORreactantIndex = 0;
	DORgroupName = DIMER_GROUP_NAME;
	DORgroupValueIndex = DimerGroup::METH_SITE_LEVEL;
	r = new NG_FreeActiveSiteBinding("FREE_B_bind_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::FREE_CHEB_bind_ACTIVE, "av", "asite");
	s->addReaction(r);
	r->setTraversalLimit(2);
	
	
	
	//R unbind the tether
	n_reactants = 1;
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheR);
	TemplateMolecule *receptor = new TemplateMolecule(recDimer);
	TemplateMolecule::bind(reactantTemplates[0], "av", receptor, "asite");
	r = new ReactionUnbinding("R_unbind_ASITE",
		n_reactants,reactantTemplates,Param::CHER_unbind_ACTIVE,"av");
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	
	//B bind the tether
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheB);
	receptor = new TemplateMolecule(recDimer);
	TemplateMolecule::bind(reactantTemplates[0], "av", receptor, "asite");
	r = new ReactionUnbinding("B_unbind_ASITE",
		n_reactants,reactantTemplates,Param::CHEB_unbind_ACTIVE,"av");
	r->setTraversalLimit(2);
	s->addReaction(r);
}


void NG::rxn_tetheredRB_bind_activeNeighbor(System *s, MoleculeType *recDimer, MoleculeType *che, bool isCheR,
	char *neighborSiteName, char *reciprocalNeighborSiteName)
{
	int neighborSite = recDimer->getBindingSiteIndex(neighborSiteName);
	int reciprocalNeighborSite= recDimer->getBindingSiteIndex(reciprocalNeighborSiteName);
	
	if(isCheR)
	{
		int n_reactants = 1;
		TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
		//First, create this receptor
		reactantTemplates[0] = new TemplateMolecule(recDimer);
		reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
		reactantTemplates[0]->addNotStateValue("m",8);  // can't bind if we are fully methylated
		//Create the neighbor that is tethered
		TemplateMolecule *recNeighbor = new TemplateMolecule(recDimer);
		TemplateMolecule::bind(reactantTemplates[0], neighborSite, recNeighbor, reciprocalNeighborSite);
		
		//Create the tethered cheR
		TemplateMolecule * cheTemp = new TemplateMolecule(che);
		cheTemp->addEmptyBindingSite("av");
		TemplateMolecule::bind(recNeighbor,"tether",cheTemp, "te");
	
		
		int DORreactantIndex = 0;
		char *DORgroupName = DIMER_GROUP_NAME;
		int DORgroupValueIndex = DimerGroup::FREE_SITE_LEVEL;
	
		ReactionClass * r = new NG_TetheredActiveNeighborSiteBinding("TETHERED_R_bind_NEIGHBOR_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
			Param::TETHERED_CHER_bind_ACTIVE, 
			neighborSiteName,
			"tether", "asite", che->getBindingSiteIndex("av"));
		s->addReaction(r);
		r->setTraversalLimit(3);
	}
	else
	{
		int n_reactants = 1;
		TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
		//First, create this receptor
		reactantTemplates[0] = new TemplateMolecule(recDimer);
		reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
		reactantTemplates[0]->addNotStateValue("m",0);  // can't bind if we are fully demethylated
		
		//Create the neighbor that is tethered
		TemplateMolecule *recNeighbor = new TemplateMolecule(recDimer);
		TemplateMolecule::bind(reactantTemplates[0], neighborSite, recNeighbor, reciprocalNeighborSite);
		
		//Create the tethered cheR
		TemplateMolecule * cheTemp = new TemplateMolecule(che);
		if(Param::useCheBFeedbackLoop) cheTemp->addStateValue("p",PHOS);
		cheTemp->addEmptyBindingSite("av");
		TemplateMolecule::bind(recNeighbor,"tether",cheTemp, "te");
	
		
		int DORreactantIndex = 0;
		char *DORgroupName = DIMER_GROUP_NAME;
		int DORgroupValueIndex = DimerGroup::METH_SITE_LEVEL;
	
		ReactionClass * r = new NG_TetheredActiveNeighborSiteBinding("TETHERED_B_bind_NEIGHBOR_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
			Param::TETHERED_CHEB_bind_ACTIVE, 
			neighborSiteName,
			"tether", "asite", che->getBindingSiteIndex("av"));
		s->addReaction(r);
		r->setTraversalLimit(3);
	}
	
	
}

void NG::rxn_tetheredRB_bind_active(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB)
{
	
	//First, free binding of CheR to the active site depends on the number of active
	//sites available.  Thus, we need to use a DOR rxn to keep track of this
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
	reactantTemplates[0]->addNotStateValue("m",8);  // can't bind if we are fully methylated
	
	TemplateMolecule * cheRtemp = new TemplateMolecule(cheR);
	cheRtemp->addEmptyBindingSite("av");  //has to have an empty site to bind the active site
	TemplateMolecule::bind(reactantTemplates[0],"tether",cheRtemp,"te");

	int DORreactantIndex = 0;
	char *DORgroupName = DIMER_GROUP_NAME;
	
	//We can only bind to an active site based on the number of active sites that are
	//available for us to bind (in this case, the free site level)
	int DORgroupValueIndex = DimerGroup::FREE_SITE_LEVEL;

	ReactionClass * r = new NG_TetheredActiveSiteBinding("TETHERED_R_bind_THIS_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::TETHERED_CHER_bind_ACTIVE, cheR->getBindingSiteIndex("av"), "asite","tether");
	s->addReaction(r);
	r->setTraversalLimit(2);
	
	
	
	
	
	
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addEmptyBindingSite("asite");  // the active site must be available
	reactantTemplates[0]->addNotStateValue("m",0);  // can't bind if we are fully demethylated
	
	TemplateMolecule * cheBtemp = new TemplateMolecule(cheB);
	if(Param::useCheBFeedbackLoop) cheBtemp->addStateValue("p",PHOS);
	cheBtemp->addEmptyBindingSite("av");  //has to have an empty site to bind the active site
	TemplateMolecule::bind(reactantTemplates[0],"tether",cheBtemp,"te");

	DORreactantIndex = 0;
	DORgroupName = DIMER_GROUP_NAME;
	
	//We can only bind to an active site based on the number of active sites that are
	//available for us to bind (in this case, the free site level)
	DORgroupValueIndex = DimerGroup::METH_SITE_LEVEL;

	r = new NG_TetheredActiveSiteBinding("TETHERED_B_bind_THIS_ASITE",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::TETHERED_CHEB_bind_ACTIVE, cheB->getBindingSiteIndex("av"), "asite","tether");
	s->addReaction(r);
	r->setTraversalLimit(2);
	
	
	
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t0", "t3");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t1", "t4");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t2", "t5");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t3", "t0");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t4", "t1");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheR, true,"t5", "t2");
	
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t0", "t3");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t1", "t4");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t2", "t5");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t3", "t0");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t4", "t1");
	rxn_tetheredRB_bind_activeNeighbor(s, recDimer, cheB, false,"t5", "t2");
}



void NG::rxn_RB_meth_demeth(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB)
{
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addNotStateValue("m",8);  // can't meth if we are already fully methylated
	
	TemplateMolecule * cheRtemp = new TemplateMolecule(cheR);
	TemplateMolecule::bind(reactantTemplates[0],"asite",cheRtemp,"av");

	int DORreactantIndex = 0;
	char *DORgroupName = CLUSTER_NAME;
	
	int DORgroupValueIndex = P_OFF_INDEX;

	ReactionClass * r = new NG_BR_MethDemeth("R_meth_RECEPTOR",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::CHER_meth_RECEPTOR, "m", "asite",1);
	s->addReaction(r);
	r->setTraversalLimit(2);
	
	
	
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(recDimer);
	reactantTemplates[0]->addNotStateValue("m",0);  // can't demeth if we are already fully methylated
	TemplateMolecule * cheBtemp = new TemplateMolecule(cheB);
	TemplateMolecule::bind(reactantTemplates[0],"asite",cheBtemp,"av");

	DORreactantIndex = 0;
	DORgroupName = CLUSTER_NAME;
	DORgroupValueIndex = P_ON_INDEX;

	r = new NG_BR_MethDemeth("B_demeth_RECEPTOR",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::CHEB_demeth_RECEPTOR, "m", "asite",-1);
	s->addReaction(r);
	r->setTraversalLimit(2);
}



void NG::rxn_CheA_phos_CheY(System * s, MoleculeType * cheA, MoleculeType *cheY)
{
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(cheA);
	reactantTemplates[0]->addStateValue("p",PHOS);  //Set constraint that cheA is at p=1
	
	reactantTemplates[1] = new TemplateMolecule(cheY);
	reactantTemplates[1]->addStateValue("p",UNPHOS);  //Set constraint that cheY is at p=0
	
	ReactionClass * r = new ReactionChangeStateOfTwoMolecules("A_PHOS_Y", n_reactants, reactantTemplates,
		Param::CHEA_phos_CHEY,"p",UNPHOS, "p", PHOS);
	r->setTraversalLimit(1);
	s->addReaction(r);
}

void NG::ANrxn_add_CheY_autodephosFull(System *s, MoleculeType * cheY)
{
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(cheY);
	reactantTemplates[0]->addStateValue("p",PHOS);  //Set constraint that we are at p=1
	
	ReactionClass * r = new ReactionChangeState("Y_AUTO_DEPHOS", n_reactants, reactantTemplates, Param::AUTO_dephos_CHEY,"p",UNPHOS);
	r->setTraversalLimit(1);
	s->addReaction(r);
}



void NG::rxn_CheA_auto_phos(System * s, MoleculeType * cheA)
{
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheA);
	reactantTemplates[0]->addStateValue("p",UNPHOS);  //Set constraint that p is unphos	
	
	int DORreactantIndex = 0;
	char *DORgroupName = CLUSTER_NAME;
	int DORgroupValueIndex = P_ON_INDEX; //Auto Phos of CheA depends on Prob of cluster to be on
	
	ReactionClass * autoPhosA = new NG_AutoPhosA("AUTO_PHOS_A",n_reactants,reactantTemplates,DORreactantIndex,DORgroupName,DORgroupValueIndex,
		Param::AUTO_phos_CHEA,"p",PHOS);
	autoPhosA->setTraversalLimit(1);
	s->addReaction(autoPhosA);
}

void NG::rxn_CheA_phos_CheBFull(System * s, MoleculeType * cheA, MoleculeType *cheB)
{
	int n_reactants = 2;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(cheA);
	reactantTemplates[0]->addStateValue("p",PHOS);  //Set constraint that cheA is at p=1
	
	reactantTemplates[1] = new TemplateMolecule(cheB);
	reactantTemplates[1]->addStateValue("p",UNPHOS); 
	
	ReactionClass * r = new ReactionChangeStateOfTwoMolecules("A_PHOS_B", n_reactants, reactantTemplates, Param::CHEA_phos_CHEB,"p",UNPHOS, "p", PHOS);
	r->setTraversalLimit(2);
	s->addReaction(r);
}

void NG::rxn_add_CheB_autodephosFull(System *s, MoleculeType * cheB)
{
	int n_reactants = 1;
	TemplateMolecule ** reactantTemplates = new TemplateMolecule *[n_reactants];
	
	reactantTemplates[0] = new TemplateMolecule(cheB);
	reactantTemplates[0]->addStateValue("p",PHOS);  //Set constraint that we are at p=1
	reactantTemplates[0]->addEmptyBindingSite("av");  //Set constraint that we are not interacting with an active site
	
	ReactionClass * r = new ReactionChangeState("B_AUTO_DEPHOS", n_reactants, reactantTemplates, Param::AUTO_dephos_CHEB,"p",UNPHOS);
	r->setTraversalLimit(2);
	s->addReaction(r);
}




















////////////////// 
NG_FreeActiveSiteBinding::NG_FreeActiveSiteBinding( char * name, 
							int n_reactants, 
							TemplateMolecule ** reactantTemplates,
							int DORreactantIndex,
							char * DORgroupName,
							int DORgroupValueIndex,
							double baseRate,
							const char *cheBindingSite,
							const char *recActiveSite): DORrxn(name,n_reactants, reactantTemplates, DORreactantIndex, DORgroupName, DORgroupValueIndex, baseRate)
{
	this->recActiveSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(recActiveSite);
	this->cheBindingSiteIndex = reactantTemplates[1]->getMoleculeType()->getBindingSiteIndex(cheBindingSite);	
}

void NG_FreeActiveSiteBinding::transformReactants(Molecule ** reactants, int nReactants)
{
	Molecule::bind(reactants[0], recActiveSiteIndex, reactants[1], cheBindingSiteIndex);
}






/* DOR rxn constructor and deconstructor */
NG_TetheredActiveSiteBinding::NG_TetheredActiveSiteBinding( char * name, 
					int n_reactants, 
					TemplateMolecule ** reactantTemplates,
					int DORreactantIndex,
					char * DORgroupName,
					int DORgroupValueIndex,
					double baseRate,
					int cheBindingSiteIndex,
					const char *recActiveSite,
					const char *tetherSite): DORrxn(name,n_reactants, reactantTemplates, DORreactantIndex, DORgroupName, DORgroupValueIndex, baseRate)
{
	this->recActiveSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(recActiveSite);
	this->tetherSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(tetherSite);
	this->cheBindingSiteIndex = cheBindingSiteIndex;
}
								
void NG_TetheredActiveSiteBinding::transformReactants(Molecule ** reactants, int nReactants)
{
	Molecule *che = reactants[0]->getBondedMolecule(this->tetherSiteIndex);
	if(che==NULL) {
		cerr<<"tethered binding to active site failed!  cheR was not there!"<<endl;
		this->printDetails();
		exit(1);
	}
	Molecule::bind(reactants[0], this->recActiveSiteIndex, che, this->cheBindingSiteIndex);
}



NG_TetheredActiveNeighborSiteBinding::NG_TetheredActiveNeighborSiteBinding( char * name, 
							int n_reactants, 
							TemplateMolecule ** reactantTemplates,
							int DORreactantIndex,
							char * DORgroupName,
							int DORgroupValueIndex,
							double baseRate,
							const char *neighborSite,
							const char *recTetherSite,
							const char *recActiveSite,
							int cheActiveSiteIndex): DORrxn(name,n_reactants, reactantTemplates, DORreactantIndex, DORgroupName, DORgroupValueIndex, baseRate)
{
	
	this->neighborSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(neighborSite);
	this->recTetherSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(recTetherSite);
	this->recActiveSiteIndex = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(recActiveSite);
	this->cheActiveSiteIndex = cheActiveSiteIndex;
}
								
void NG_TetheredActiveNeighborSiteBinding::transformReactants(Molecule ** reactants, int nReactants)
{
	Molecule *che = reactants[0]->getBondedMolecule(neighborSiteIndex)->getBondedMolecule(recTetherSiteIndex);
	
	if(reactants[0]->isBindingSiteOpen(neighborSiteIndex)) {
		cerr<<"tethered binding to active neighbor site failed!  the neighbor that should be bound to cheR/B  wasn't there!"<<endl;
		this->printDetails();
		exit(1);
	}
	
	if(reactants[0]->getBondedMolecule(neighborSiteIndex)->isBindingSiteOpen(recTetherSiteIndex)) {
		cerr<<"tethered binding to active neighbor site failed!  the neighbor that should be bound to cheR/B  wasn't there!"<<endl;
		this->printDetails();
		exit(1);
	}
	
	
	
	Molecule::bind(reactants[0],recActiveSiteIndex, che, cheActiveSiteIndex);

}
		




NG_BR_MethDemeth::NG_BR_MethDemeth( char * name, 
							int n_reactants, 
							TemplateMolecule ** reactantTemplates,
							int DORreactantIndex,
							char * DORgroupName,
							int DORgroupValueIndex,
							double baseRate,
							const char *receptorMethStateName,
							const char *receptorActiveSiteName,
							int deltaMethState):DORrxn(name,n_reactants, reactantTemplates, DORreactantIndex, DORgroupName, DORgroupValueIndex, baseRate)
{
	this->receptorMethStateIndex= reactantTemplates[0]->getMoleculeType()->getStateIndex(receptorMethStateName);
	this->receptorActiveSiteName = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(receptorActiveSiteName);
	this->deltaMethState = deltaMethState;
}
								
void NG_BR_MethDemeth::transformReactants(Molecule ** reactants, int nReactants)
{
	//1st, update meth state
	reactants[0]->setState(receptorMethStateIndex, reactants[0]->getState(receptorMethStateIndex)+deltaMethState);
	Molecule::unbind(reactants[0],receptorActiveSiteName);
}




//////////////////////////  Methylation that does not depend on activity
NG_BR_MethDemeth2::NG_BR_MethDemeth2( char * name, 
							int n_reactants, 
							TemplateMolecule ** reactantTemplates,
							double rate,
							const char *receptorMethStateName,
							const char *receptorActiveSiteName,
							int deltaMethState):ReactionClass(name, n_reactants, reactantTemplates, rate)
{
	this->receptorMethStateIndex= reactantTemplates[0]->getMoleculeType()->getStateIndex(receptorMethStateName);
	this->receptorActiveSiteName = reactantTemplates[0]->getMoleculeType()->getBindingSiteIndex(receptorActiveSiteName);
	this->deltaMethState = deltaMethState;
}
								
void NG_BR_MethDemeth2::transformReactants(Molecule ** reactants, int nReactants)
{
	//1st, update meth state
	reactants[0]->setState(receptorMethStateIndex, reactants[0]->getState(receptorMethStateIndex)+deltaMethState);
	Molecule::unbind(reactants[0],receptorActiveSiteName);
}














void NG::rxn_retether(System *s, MoleculeType * recDimer, MoleculeType * cheR, MoleculeType * cheB)
{
	//R unbind the tether
	int n_reactants = 1;
	TemplateMolecule **reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheR);
	reactantTemplates[0]->addEmptyBindingSite("te");
	
	TemplateMolecule *receptor = new TemplateMolecule(recDimer);
	receptor->addEmptyBindingSite("tether");
	TemplateMolecule::bind(reactantTemplates[0], "av", receptor, "asite");
	
	ReactionClass *r = new NG_retether("R_retether_RECEPTOR",
		n_reactants,reactantTemplates,Param::ACTIVE_BOUND_CHER_bind_TETHER,
		recDimer->getBindingSiteIndex("tether"),
		cheR->getBindingSiteIndex("av"),
		cheR->getBindingSiteIndex("te"));
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	//R unbind the tether
	reactantTemplates = new TemplateMolecule *[n_reactants];
	reactantTemplates[0] = new TemplateMolecule(cheB);
	reactantTemplates[0]->addEmptyBindingSite("te");
	
	receptor = new TemplateMolecule(recDimer);
	receptor->addEmptyBindingSite("tether");
	TemplateMolecule::bind(reactantTemplates[0], "av", receptor, "asite");
	
	r = new NG_retether("ACTIVE_BOUND_CHEB_bind_TETHER",
		n_reactants,reactantTemplates,Param::ACTIVE_BOUND_CHEB_bind_TETHER,
		recDimer->getBindingSiteIndex("tether"),
		cheB->getBindingSiteIndex("av"),
		cheB->getBindingSiteIndex("te"));
	r->setTraversalLimit(2);
	s->addReaction(r);
	
	
	
}











NG_retether::NG_retether(	char * name, 
								int n_reactants, 
								TemplateMolecule ** reactantTemplates, 
								double rate,
								int recTetherSiteIndex,
								int cheActiveSiteIndex,
								int cheTetherSiteIndex): ReactionClass(name, n_reactants, reactantTemplates, rate)
{
	this->recTetherSiteIndex = recTetherSiteIndex;
	this->cheActiveSiteIndex = cheActiveSiteIndex;
	this->cheTetherSiteIndex = cheTetherSiteIndex;
}
void NG_retether::transformReactants(Molecule ** reactants, int nReactants)
{
	Molecule *receptor = reactants[0]->getBondedMolecule(cheActiveSiteIndex);
	Molecule::bind(reactants[0],cheTetherSiteIndex, receptor, recTetherSiteIndex);
}





//////////////////  CheA autophosphorylation
NG_AutoPhosA::NG_AutoPhosA( char * name, 
							int n_reactants, 
							TemplateMolecule ** reactantTemplates,
							int DORreactantIndex,
							char * DORgroupName,
							int DORgroupValueIndex,
							double baseRate,
							const char * stateName,
							int newStateValue): DORrxn(name,1, reactantTemplates, DORreactantIndex, DORgroupName, DORgroupValueIndex, baseRate)
{
	this->stateIndex = reactantTemplates[0]->getMoleculeType()->getStateIndex(stateName);
	this->newStateValue	= newStateValue;
}

void NG_AutoPhosA::transformReactants(Molecule ** reactants, int nReactants)
{
	reactants[0]->setState(stateIndex, newStateValue);
}




