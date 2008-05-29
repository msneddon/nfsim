#include <math.h>
#include "receptor_cluster.hh"





using namespace ReceptorCluster;




RecCluster::RecCluster(const char *name, System *s, int methStateIndex, NGparam &p) : Group(name, s, methStateIndex)
{
	//First, initialize the free energy offsets
	TAR_freeEnergyOffset = new double[9];
	TSR_freeEnergyOffset = new double[9];
	
	for(int m=0; m<=8; m++) {
		TAR_freeEnergyOffset[m] = p.get_clusterTarFreeEnergyOffset(m);
		TSR_freeEnergyOffset[m] = p.get_clusterTsrFreeEnergyOffset(m);
	}
	
		
	//Of course we start with nothing in the cluster
	this->TAR_COUNT = 0;
	this->TSR_COUNT = 0;
	this->CHE_A_COUNT = 0;
			
	//Remember the ligand concentration and K values
	this->AspConc = p.get_aspartateConcentration();
	this->Asp_Koff_TAR = p.get_clusterTarAspKoff();
	this->Asp_Kon_TAR = p.get_clusterTarAspKon();
	this->Asp_Koff_TSR = p.get_clusterTsrAspKoff();
	this->Asp_Kon_TSR = p.get_clusterTsrAspKon();
		
	
	//Keep track of several values throughout the simulation
	this->TAR_LOG_TERM = 0;
	this->TSR_LOG_TERM = 0;
	updateTarLogValue();
	updateTsrLogValue();
	this->FreeEnergyOffsetSum = 0;
	
	//Value vector contains our P_ON and P_OFF values, so we have to
	//add those two values to the value vector.  I clear it first just for fun.
	value.clear();
	value.push_back(1);  //pushing P_on
	value.push_back(1);  //pushing P_off
}



RecCluster::RecCluster(char *name, System *s, int methStateIndex,
			double Asp_Koff_TAR,
			double Asp_Kon_TAR,
			double Asp_Koff_TSR,
			double Asp_Kon_TSR,
			double Asp_Conc ) : Group(name, s, methStateIndex)
{
	//First, initialize the free energy offsets
	this->initFreeEnergyOffset();
	
	//Of course we start with nothing in the cluster
	this->TAR_COUNT = 0;
	this->TSR_COUNT = 0;
	this->CHE_A_COUNT = 0;
		
	//Remember the ligand concentration and K values
	this->AspConc = Asp_Conc;
	this->Asp_Koff_TAR = Asp_Koff_TAR;
	this->Asp_Kon_TAR = Asp_Kon_TAR;
	this->Asp_Koff_TSR = Asp_Koff_TSR;
	this->Asp_Kon_TSR = Asp_Kon_TSR;
	
	//Keep track of several values throughout the simulation
	this->TAR_LOG_TERM = 0;
	this->TSR_LOG_TERM = 0;
	updateTarLogValue();
	updateTsrLogValue();
	this->FreeEnergyOffsetSum = 0;
	
	//Value vector contains our P_ON and P_OFF values, so we have to
	//add those two values to the value vector.  I clear it first just for fun.
	value.clear();
	value.push_back(1);  //pushing P_on
	value.push_back(1);  //pushing P_off
}

void RecCluster::initFreeEnergyOffset()
{
	//Free Energy Values in units KbT taken from:
	//   Endres, Precise adaptation in bacterial chemotaxis through "assistance neighborhoods"
	//   PNAS, Aug 29, 2006, vol 103, no 35, 13040-13044
	double offsetFactor = 1;
	
	TAR_freeEnergyOffset = new double[9];
	TAR_freeEnergyOffset[0] =  1.0*offsetFactor;
	TAR_freeEnergyOffset[1] =  0.5*offsetFactor;
	TAR_freeEnergyOffset[2] =  0.0*offsetFactor;
	TAR_freeEnergyOffset[3] = -0.3*offsetFactor;
	TAR_freeEnergyOffset[4] = -0.6*offsetFactor;
	TAR_freeEnergyOffset[5] = -0.85*offsetFactor;
	TAR_freeEnergyOffset[6] = -1.1*offsetFactor;
	TAR_freeEnergyOffset[7] = -2.0*offsetFactor;
	TAR_freeEnergyOffset[8] = -3.0*offsetFactor;
	
	TSR_freeEnergyOffset = new double[9];
	TSR_freeEnergyOffset[0] =  1.0*offsetFactor;
	TSR_freeEnergyOffset[1] =  0.5*offsetFactor;
	TSR_freeEnergyOffset[2] =  0.0*offsetFactor;
	TSR_freeEnergyOffset[3] = -0.3*offsetFactor;
	TSR_freeEnergyOffset[4] = -0.6*offsetFactor;
	TSR_freeEnergyOffset[5] = -0.85*offsetFactor;
	TSR_freeEnergyOffset[6] = -1.1*offsetFactor;
	TSR_freeEnergyOffset[7] = -2.0*offsetFactor;
	TSR_freeEnergyOffset[8] = -3.0*offsetFactor;
	

	//Will's parameters
//	double e0 = 3.36/19;
//	double e1 = -0.063;
//	TAR_freeEnergyOffset = new double[9];
//	TSR_freeEnergyOffset = new double[9];
//	for(int m=0; m<=8; m++)
//	{
//		TAR_freeEnergyOffset[m] = e0+e1*m;
//		TSR_freeEnergyOffset[m] = e0+e1*m;
//	}
	
}

RecCluster::~RecCluster() 
{
	delete [] TAR_freeEnergyOffset;
	delete [] TSR_freeEnergyOffset;
}

//When we call this, we are always adding receptor dimers!!!
//to add CheA to this group, call the function addToGroupWithoutListener
//so that we don't listen for CheA's methylation state
void RecCluster::addToGroup(Molecule * m)
{
	//cout<<"Trying to add: "<<endl; m->printDetails();
	
	
	//Depending on the receptor type,recalculate the log value, free
	//energy offset, and update the counts
	int ReceptorType = m->getState(m->getMoleculeType()->getStateIndex("type"));
	//cout<<"ReceptorType: " << ReceptorType<<endl;
	if(ReceptorType==TAR)
	{
		TAR_COUNT++;
		FreeEnergyOffsetSum += TAR_freeEnergyOffset[m->getState(stateIndex)];
		updateTarLogValue();
	}
	else if(ReceptorType==TSR)
	{
		TSR_COUNT++;
		FreeEnergyOffsetSum += TSR_freeEnergyOffset[m->getState(stateIndex)];
		updateTsrLogValue();
	}
	else
	{
		cerr<<"Trying to add: "<<endl; m->printDetails();
		cerr<<"Which is not a receptor of type TAR or TSR.  I just can't do that yet."<<endl;
		exit(1);
	}
	
	//Add this molecule to the group if we have successfully updated everything
	this->groupMembers.push_back(m);
	
	//Create a StateChangeListener to listen to the molecule
	new StateChangeListener(m,this,this->stateIndex);
	
	//determine the new p_on, which is our value (rateFactor) for this group
	double p_on = 1.0 / (1.0 + (exp(FreeEnergyOffsetSum) * TAR_LOG_TERM * TSR_LOG_TERM));
	
	
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		cerr<<"FreeEnergyOffsetSum: " << FreeEnergyOffsetSum <<endl;
		cerr<<"TAR_LOG_TERM: " << TAR_LOG_TERM <<endl;
		cerr<<"TSR_LOG_TERM: " << TSR_LOG_TERM <<endl;
		exit(1);
	}
	
	//and update the value of this group
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = (1-p_on);
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false;
//	cout<<"Free Energy Offset: " << FreeEnergyOffsetSum;
//	cout<<"  Tar Log Term: " << TAR_LOG_TERM<<" Tsr Log Term: " << TSR_LOG_TERM <<endl;
//	cout<<"Adding receptor of type: " << ReceptorType <<", ON: "<<p_on<<"  OFF: "<<(1-p_on)<<endl;
}



void RecCluster::updateTarLogValue()
{
	TAR_LOG_TERM = pow((1.0+(AspConc/Asp_Koff_TAR)) / (1.0+(AspConc/Asp_Kon_TAR)) ,TAR_COUNT);
//	cout<<"Updating Tar Log Term: " << TAR_LOG_TERM << endl;
}

void RecCluster::updateTsrLogValue()
{
	//TSR_LOG_TERM = pow( (Asp_Kon_TSR),TSR_COUNT);
	TSR_LOG_TERM = pow(  (1.0+(AspConc/Asp_Koff_TSR)) / (1.0+(AspConc/Asp_Kon_TSR))  ,TSR_COUNT);
//	cout<<"Updating Tsr Log Term: " << TSR_LOG_TERM << endl;
}


/* this gets called by a stateChangeListener when a state changes.  Remember - it
 * doesn't update reaction rates!  You have to update ReactionRates with a separate
 * call to updateReactionRates! (This is because we need to finish firing the 
 * reaction before we update rates)  */
void RecCluster::notify(Molecule *changedMolecule, int oldStateValue, int newStateValue)
{
	//cout<<"Notifying.."<<endl;
	
	//Depending on the type of receptor that changed, we have to update
	//the free energy offsets accordingly.  The ligand concentration hasn't changed
	//so we don't actually have to update the log values again.
	int ReceptorType = changedMolecule->getMoleculeType()->getStateIndex("type");
	
	if(newStateValue>8) { cout<<"Error in group!::"<<oldStateValue<<"   "<<newStateValue<<endl; exit(1); }
	if(newStateValue<0) { cout<<"Error in group!::"<<oldStateValue<<"   "<<newStateValue<<endl; exit(1); }
	if(ReceptorType==TAR)
	{
		FreeEnergyOffsetSum -= TAR_freeEnergyOffset[oldStateValue];
		FreeEnergyOffsetSum += TAR_freeEnergyOffset[newStateValue];
		this->updateTsrLogValue();
	}
	else if(ReceptorType==TSR)
	{
		FreeEnergyOffsetSum -= TSR_freeEnergyOffset[oldStateValue];
		FreeEnergyOffsetSum += TSR_freeEnergyOffset[newStateValue];
		this->updateTsrLogValue();
	}
	else
	{
		cerr<<"How did you manage to get here?!?  You are notifying a receptor cluster with ";
		cerr<<"A molecule that is not a receptor!"<<endl;
		exit(1);
	}
	
	//determine the new p_on, which is our value (rateFactor) for this group
	double p_on = 1 / (1 + (exp(FreeEnergyOffsetSum) * TAR_LOG_TERM * TSR_LOG_TERM));
	
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		cerr<<"FreeEnergyOffsetSum: " << FreeEnergyOffsetSum <<endl;
		cerr<<"TAR_LOG_TERM: " << TAR_LOG_TERM <<endl;
		cerr<<"TSR_LOG_TERM: " << TSR_LOG_TERM <<endl;
		exit(1);
	}
	
	
	//and update the value of this group
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = (1-p_on);
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false; 
	
	
	//cout<<"done"<<endl;
//////////////OLD stuff	

  // System *s = changedMolecule->getMoleculeType()->getSystem();
 ////  if(s->getCurrentTime() >1000)
  // {
 //  		s->printAllGroups();
  // 		exit(0);
  // }

   //cout<<"    Group " << groupName<<"_"<<Group_ID<<" was alerted of change. ";
//	
//	
//	//determin p_on and remember it as our value (rateFactor)
//	//double p_on = 1 / (1 + exp(energySum + logValue));
//	
//	
//	areReactionsUpToDate = false;
//	cout<<" and is: " << energySum;
//	
//	cout<<" log value is: " << logValue;
//	cout<<" value (p_on) was "<<value.at(P_ON_INDEX);
//	value.at(P_ON_INDEX) = p_on;
//	cout<<" and is "<<value.at(P_ON_INDEX)<<endl;
//	
//	cout<<" value (p_off) was "<<value.at(P_OFF_INDEX);
//	value.at(P_OFF_INDEX) = (1-p_on);
//	cout<<" and is "<<value.at(P_OFF_INDEX)<<endl;
}

//This is the function that gets called if a ligand concentration
//is changed.  By convention, values will be the same length as the number
//of ligand molecules we are trying to sense.  For now, this means that
//it is only of length one, and that one value is Aspartate concentration
void RecCluster::updateGroupProperty(double * values, int n_values)
{
	//By convention, value will be of length one and will be the ligand value in molar
	this->AspConc = values[0];
	updateTarLogValue();
	updateTsrLogValue();
	
	//determine the new p_on, which is our value (rateFactor) for this group
	double p_on = 1 / (1 + (exp(FreeEnergyOffsetSum) * TAR_LOG_TERM * TSR_LOG_TERM));
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		cerr<<"FreeEnergyOffsetSum: " << FreeEnergyOffsetSum <<endl;
		cerr<<"TAR_LOG_TERM: " << TAR_LOG_TERM <<endl;
		cerr<<"TSR_LOG_TERM: " << TSR_LOG_TERM <<endl;
		exit(1);
	}
	
	
	//and update the value of this group
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = (1-p_on);
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false; 
	
	//Finally, update the rates accordingly
	updateReactionRates();
}



void RecCluster::printDetails()
{
	cout<<"AN_Cluster " << this->Group_ID <<" : " <<endl;
	cout<<"   -Tar Count: " << this->TAR_COUNT <<", TSR Count: " <<this->TSR_COUNT;
	cout<<", CheA Count: "<<this->CHE_A_COUNT<<endl;
	cout<<"   -AspConcentration: " << this->AspConc <<" M "<<endl;
	
	
	int methSumTar = 0;
	int methSumTsr = 0;
	for(molIter = groupMembers.begin(); molIter != groupMembers.end(); molIter++ )
	{
		if((*molIter)->getMoleculeTypeName()=="ReceptorDimer") 
		{
		int type = (*molIter)->getState((*molIter)->getMoleculeType()->getStateIndex("type"));
		if(type==TAR)
		{
			methSumTar+=(*molIter)->getState(stateIndex);
			//cout<<"              tar: "<<(*molIter)->getState(stateIndex)<< endl;
		}
		else if(type==TSR)
		{
			methSumTsr+=(*molIter)->getState(stateIndex);
			//cout<<"tsr: "<<(*molIter)->getState(stateIndex)<< endl;
		}
		}
	}
	cout<<"   -TarMethLevel = "<<methSumTar<<", TsrMethLevel = " << methSumTsr;
	cout<<",   Total: "<<methSumTar+methSumTsr<<endl;
	cout<<"   -P_On = " << value.at(P_ON_INDEX) << ", P_Off: " << value.at(P_OFF_INDEX)<<endl;
	cout<<"   -Free Energy Offset Sum: " << this->FreeEnergyOffsetSum;
	cout<<", LogTarTerm = "<<this->TAR_LOG_TERM << ", LogTsrTerm = "<<this->TSR_LOG_TERM<<endl;
	cout<<endl;
	
	
}




////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

RecClusterWill::RecClusterWill(string name, System *s, int methStateIndex, NGparam &p) : Group(name, s, methStateIndex)
{

	this->e0 = p.get_clusterEnergy_e0();
	this->e1 = p.get_clusterEnergy_e1();
	
	
	this->methylationLevel = 0;
	
		
	//Of course we start with nothing in the cluster
	this->TAR_COUNT = 0;
	this->TSR_COUNT = 0;
	this->CHE_A_COUNT = 0;
			
	//Remember the ligand concentration and K values
	this->AspConc = p.get_aspartateConcentration();
	this->Asp_Koff_TAR = p.get_clusterTarAspKoff();
	this->Asp_Kon_TAR = p.get_clusterTarAspKon();
	this->Asp_Koff_TSR = p.get_clusterTsrAspKoff();
	this->Asp_Kon_TSR = p.get_clusterTsrAspKon();

	//cout<<AspConc<<"\t"<<Asp_Koff_TAR<<"\t"<<Asp_Kon_TAR<<"\t"<<Asp_Koff_TSR<<"\t"<<Asp_Kon_TSR<<endl;
	
	
	//Value vector contains our P_ON and P_OFF values, so we have to
	//add those two values to the value vector.  I clear it first just for fun.
	value.clear();
	value.push_back(1);  //pushing P_on
	value.push_back(1);  //pushing P_off
}





RecClusterWill::~RecClusterWill() 
{
}

//When we call this, we are always adding receptor dimers!!!
//to add CheA to this group, call the function addToGroupWithoutListener
//so that we don't listen for CheA's methylation state
void RecClusterWill::addToGroup(Molecule * m)
{
	//cout<<"Trying to add: "<<endl; m->printDetails();
	
	
	//Depending on the receptor type,recalculate the log value, free
	//energy offset, and update the counts
	int ReceptorType = m->getState(m->getMoleculeType()->getStateIndex("type"));
	//cout<<"ReceptorType: " << ReceptorType<<endl;
	if(ReceptorType==TAR)
	{
		TAR_COUNT++;
		this->methylationLevel = methylationLevel+m->getState(stateIndex);
	}
	else if(ReceptorType==TSR)
	{
		TSR_COUNT++;
		this->methylationLevel = methylationLevel+m->getState(stateIndex);
	}
	else
	{
		cerr<<"Trying to add: "<<endl; m->printDetails();
		cerr<<"Which is not a receptor of type TAR or TSR.  I just can't do that yet."<<endl;
		exit(1);
	}
	
	//Add this molecule to the group if we have successfully updated everything
	this->groupMembers.push_back(m);
	
	//Create a StateChangeListener to listen to the molecule
	new StateChangeListener(m,this,this->stateIndex);
	
	//determine the new p_on, which is our value (rateFactor) for this group
	double p_on = this->getPon();
	
	
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		exit(1);
	}
	
	//and update the value of this group
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = (1-p_on);
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false;
//	cout<<"Free Energy Offset: " << FreeEnergyOffsetSum;
//	cout<<"  Tar Log Term: " << TAR_LOG_TERM<<" Tsr Log Term: " << TSR_LOG_TERM <<endl;
//	cout<<"Adding receptor of type: " << ReceptorType <<", ON: "<<p_on<<"  OFF: "<<(1-p_on)<<endl;
}






/* this gets called by a stateChangeListener when a state changes.  Remember - it
 * doesn't update reaction rates!  You have to update ReactionRates with a separate
 * call to updateReactionRates! (This is because we need to finish firing the 
 * reaction before we update rates)  */
void RecClusterWill::notify(Molecule *changedMolecule, int oldStateValue, int newStateValue)
{
	//cout<<"Notifying.."<<endl;
	
	//Depending on the type of receptor that changed, we have to update
	//the free energy offsets accordingly.  The ligand concentration hasn't changed
	//so we don't actually have to update the log values again.
	
	//if(newStateValue>8) { cout<<"Error in group!::"<<oldStateValue<<"   "<<newStateValue<<endl; exit(1); }
	if(newStateValue<0) { cout<<"Error in group!::"<<oldStateValue<<"   "<<newStateValue<<endl; exit(1); }

	methylationLevel = methylationLevel - oldStateValue;
	methylationLevel = methylationLevel + newStateValue;
	
	//determine the new p_on, which is our value (rateFactor) for this group
	double p_on = getPon();
	double p_off = 1-p_on;
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		exit(1);
	}
	
	
	//and update the value of this group
	//cout<<p_on<<"\t"<<methylationLevel<<endl;
	//cout<<AspConc<<endl;
	//cout<<getFofL()<<endl;
	//cout<<e1*methylationLevel<<endl;
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = p_off;
	
	
	
//	cout<<"e0:"<<this->e0<<endl;
//	cout<<"e1:"<<this->e1<<endl;
//	
//	for(int i=0; i<2000; i+=10)
//	{
//		methylationLevel = i;
//		p_on = getPon();
//		cout<<i<<",\t"<<p_on<<";\n";
//	}
//	
//	exit(1);
//	cout<<"Updating activity: "<<p_on<<"  methLevel: "<<this->methylationLevel<<endl;
//	cout<<this->TAR_COUNT<<this->TSR_COUNT<<endl;
	
	
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false; 
	
	//cout<<methylationLevel<<"\t"<<AspConc<<"\t"<<p_on<<endl;
}

//This is the function that gets called if a ligand concentration
//is changed.  By convention, values will be the same length as the number
//of ligand molecules we are trying to sense.  For now, this means that
//it is only of length one, and that one value is Aspartate concentration
void RecClusterWill::updateGroupProperty(double * values, int n_values)
{
	//By convention, value will be of length one and will be the ligand value in molar
	this->AspConc = values[0];
	
	
	double p_on = getPon();
	double p_off = 1-p_on;
	
	if(!p_on>0 && !p_on<1) {
		cerr<<"Somehow we got an NAN when calculating the probability of a cluster to be on."<<endl;
		exit(1);
	}
	
	
	//and update the value of this group
	value.at(P_ON_INDEX) = p_on;
	value.at(P_OFF_INDEX) = p_off;
	
	//Things are no longer up to date, so remember to update them
	areReactionsUpToDate = false; 
	
	//Finally, update the rates accordingly
	updateReactionRates();
	//cout<<AspConc<<"\t"<<p_on<<endl;
}



void RecClusterWill::printDetails()
{
	cout<<"AN_Cluster " << this->Group_ID <<" : " <<endl;
	cout<<"   -Tar Count: " << this->TAR_COUNT <<", TSR Count: " <<this->TSR_COUNT;
	cout<<", CheA Count: "<<this->CHE_A_COUNT<<endl;
	cout<<"   -AspConcentration: " << this->AspConc <<" M "<<endl;
	
	
	int methSumTar = 0;
	int methSumTsr = 0;
	for(molIter = groupMembers.begin(); molIter != groupMembers.end(); molIter++ )
	{
		if((*molIter)->getMoleculeTypeName()=="ReceptorDimer") 
		{
		int type = (*molIter)->getState((*molIter)->getMoleculeType()->getStateIndex("type"));
		if(type==TAR)
		{
			methSumTar+=(*molIter)->getState(stateIndex);
			//cout<<"              tar: "<<(*molIter)->getState(stateIndex)<< endl;
		}
		else if(type==TSR)
		{
			methSumTsr+=(*molIter)->getState(stateIndex);
			//cout<<"tsr: "<<(*molIter)->getState(stateIndex)<< endl;
		}
		}
	}
	cout<<"   -TarMethLevel = "<<methSumTar<<", TsrMethLevel = " << methSumTsr;
	cout<<",   Total: "<<methSumTar+methSumTsr<<endl;
	cout<<"   -P_On = " << value.at(P_ON_INDEX) << ", P_Off: " << value.at(P_OFF_INDEX)<<endl;
	cout<<endl;
}







double RecClusterWill::getPon()
{
	return 1.0/( 1.0 + (exp(e1*methylationLevel) * getFofL()));
}


double RecClusterWill::getFofL()
{
	double tarTerm = pow(  (1.0+(AspConc/Asp_Koff_TAR)) / (1.0+(AspConc/Asp_Kon_TAR))  ,TAR_COUNT);
	double tsrTerm = pow(  (1.0+(AspConc/Asp_Koff_TSR)) / (1.0+(AspConc/Asp_Kon_TSR))  ,TSR_COUNT);
	return exp(e0)*tarTerm*tsrTerm;
}














