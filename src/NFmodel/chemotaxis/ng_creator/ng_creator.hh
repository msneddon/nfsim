#ifndef NG_CREATOR_HH_
#define NG_CREATOR_HH_

#include "../receptor_cluster/receptor_cluster.hh"
#include "../ng/ng.hh"
#include "../chemotaxis.hh"
#include "../ng_param/ng_param.hh"


System * create(
		const char *systemName, 
		NGparam &p);



MoleculeType * createReceptorDimer(System *s, NGparam &p);
MoleculeType * createCheA(System *s, NGparam &p);
MoleculeType * createCheR(System *s, NGparam &p);
MoleculeType * createCheB(System *s, NGparam &p);
MoleculeType * createCheY(System *s, NGparam &p);
MoleculeType * createMotor(System *s, NGparam &p);


void createRxn_freeR_bind_unbind_active(System *s, 
		MoleculeType * receptorDimer, 
		MoleculeType * cheR,
		NGparam &p);

void createRxn_R_meth(System *s, 
		MoleculeType * receptorDimer, 
		MoleculeType * cheR,
		NGparam &p);

void createRxn_freeB_bind_unbind_active(System *s, 
		MoleculeType * receptorDimer, 
		MoleculeType * cheB,
		NGparam &p);

void createRxn_B_demeth(System *s, 
		MoleculeType * receptorDimer, 
		MoleculeType * cheB,
		NGparam &p);

void createRxn_CheA_auto_phos(System * s, 
		MoleculeType * cheA, 
		NGparam &p);

void createRxn_CheA_phos_CheY(System * s, 
		MoleculeType * cheA, 
		MoleculeType *cheY,
		NGparam &p);



void createRxn_add_CheY_autodephos(System *s, 
		MoleculeType * cheY,
		NGparam &p);

void createRxn_CheA_phos_CheB(System * s, 
					MoleculeType * cheA, 
					MoleculeType *cheB,
					NGparam &p);

void createRxn_add_CheB_autodephos(System *s, 
					MoleculeType * cheB,
					NGparam &p);


void createObservable_ReceptorMethLevel(MoleculeType *recDimer, NGparam &p);
void createObservable_MotorState(MoleculeType * motor, NGparam &p);


#endif /*NG_CREATOR_HH_*/
