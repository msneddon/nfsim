#ifndef TLBR_HH_
#define TLBR_HH_





#include "../../NFinput/NFinput.hh"
#include "../../NFcore/NFcore.hh"
#include "../../NFutil/NFutil.hh"

//!  Tests for running the TLBR system.
/*!
    @author Michael Sneddon
 */
namespace NFtest_tlbr {

	void run(map<string,string> &argMap);

	void runSystem(int n_L, int n_R, double cTot, double beta, double koff,
			string outputFileName, double simTime, double dt,
			bool outputObservables, bool outputAvg, bool outputFinalDist);

	MoleculeType * createL(System * s, int count);
	MoleculeType * createR(System * s, int count);

	void createFreeBindingRxns(System * s, MoleculeType * L, MoleculeType * R, double rate);
	void createUnbindingRxns(System * s, MoleculeType * R, double rate);
	void createCrossLinkingRxns(System * s, MoleculeType * L, MoleculeType *R, double rate);
}













#endif /*TLBR_HH_*/
