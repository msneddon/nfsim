/**
 * AN_clusters.cpp
 * This file contains functions to create clusters
 * 
 * 
 * 
 */

#include <string>
#include "receptor_cluster.hh"


using namespace ReceptorCluster;


/*
 * Creates as many clusters as possible of the same size, choosing the members of each cluster
 * randomly from the two vectors such that the number of tar or tsr dimers per cluster is
 * maintained.
 * 
 * Through the 'withAN' parameter, you can either create the clusters with the hex grid bond
 * pattern, or, if set to false, will only add them to clusters.
 * 
 * This returns the number of clusters created, or -1 if
 */
int ReceptorCluster::createRandomFixedSizeClusters(vector <Molecule *> &tarDimers,
	vector <Molecule *> &tsrDimers, int tarPerCluster, int tsrPerCluster,bool withAN)
{
	
	return 0;
}







int ReceptorCluster::createFixedSizeClustersWithConstantRatio(System *s, vector <Molecule *> &receptorDimers, NGparam &p)
{
	//First, calculate the size of the clusters that the dimers should be
	int size = p.get_clusterTarCountPerCluster() + p.get_clusterTsrCountPerCluster();
	
	//first create two new vectors to contain either Tar or Tsr
	vector <Molecule *>::iterator molIter;
	vector <Molecule *> tars;
	vector <Molecule *> tsrs;
	for(molIter=receptorDimers.begin(); molIter!=receptorDimers.end(); molIter++)
	{
		int ReceptorType = (*molIter)->getState((*molIter)->getMoleculeType()->getStateIndex("type"));
		if(ReceptorType==TAR)
			tars.push_back((*molIter));
		else if(ReceptorType==TSR)
			tsrs.push_back((*molIter));
		else
		{
			cerr<<"You are trying to create a cluster (in file AN_clusters.cpp) that has a ";
			cerr<<"receptor that is not a TAR or TSR!!!"<<endl;
			cerr<<"What were you thinking!?!?"<<endl;
			exit(1);	
		}	
	}
	
//	cout<<endl<<"Creating clusters of fixed size and fixed composition of size: " << size<<endl;
//	cout<<"Number of Tar Dimers: "<<tars.size()<<endl;
//	cout<<"Number of Tsr Dimers: "<<tsrs.size()<<endl;
//	cout<<"Looks like you can create "<<receptorDimers.size()/size<<"  clusters";
//	cout<<" with " << receptorDimers.size()%size << " leftover"<<endl;
//	cout<<"There will be " << tars.size()/(receptorDimers.size()/size) <<" Tar Dimers per cluster and"<<endl;
//	cout<<"There will be " << tsrs.size()/(receptorDimers.size()/size) <<" Tsr Dimers per cluster"<<endl;
	
	//now, check if we can create even clusters, if not, don't continue
	if(receptorDimers.size()%size!=0)
	{
		cerr<<"!!  You are trying to create receptor clusters of size "<<size<<", but you have ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		cerr<<"!!  That gives you "<< receptorDimers.size()%size <<" leftover! "<<endl;
		exit(1);
	}
	if(tars.size()%(receptorDimers.size()/size)!=0)
	{
		cerr<<"!!  You are trying to create "<<receptorDimers.size()/size<<" clusters, but you have ";
		cerr<<tars.size()<<" tar receptors to go around evenly!! " <<endl;
		cerr<<"!!  That gives you "<< tars.size()%(receptorDimers.size()/size) <<" tar leftover! "<<endl;
		exit(1);
	}
	if(tsrs.size()%(receptorDimers.size()/size)!=0)
	{
		cerr<<"!!  You are trying to create "<<receptorDimers.size()/size<<" clusters, but you have ";
		cerr<<tsrs.size()<<" tsr receptors to go around evenly!! " <<endl;
		cerr<<"!!  That gives you "<< tsrs.size()%(receptorDimers.size()/size) <<" tsr leftover! "<<endl;
		exit(1);
	}
	if(receptorDimers.size()<=0)
	{
		cerr<<"!!  Can't create a receptor cluster with no receptors! You gave me ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		exit(1);
	}
	
	//Some initial variables to get us started
	int clusterCount = receptorDimers.size()/size;
	int lastIndexTar = tars.size()-1;
	int lastIndexTsr = tsrs.size()-1;
	int randDimer = -1;
	int methylStateIndex = receptorDimers.at(0)->getMoleculeType()->getStateIndex(p.get_nameReceptorMethState());
	int tarPerCluster = tars.size()/(receptorDimers.size()/size);
	int tsrPerCluster = tsrs.size()/(receptorDimers.size()/size);
	
	vector <Molecule *> currentClusterMembers;
	//vector <Molecule *>::iterator molIter;
	vector <int> bondSeq;
	
	//We have to create a bond sequence if we want to create the hexgrid
	if(p.get_useNeighborRxns())
		ReceptorCluster::createStandardHexBondTraversalSeq(bondSeq,0, size);
	
	//Populate all the clusters
	cout<<"Making "<<clusterCount<<" Cluster(s)..."<<endl;
	
	for(int c=0; c<clusterCount; c++)
	{
		//cout<<endl<<"populating cluster " << c+1 << endl;
		
		//First do tars
		for(int k=0; k<tarPerCluster; k++)
		{
			//1) randomly select a dimer that has not been selected	
			randDimer = NFutil::RANDOM_INT(0,lastIndexTar+1);
			//cout<<"   -selecting: " << randDimer << " / " << lastIndexTar;
			//cout<<"  which is molecule: " << tars.at(randDimer)->getUniqueID() << endl;
			
			//2) put that dimer into a list
			currentClusterMembers.push_back(tars.at(randDimer));
			
			// swap, so that the random element moves to the end and we don't select it again
			Molecule * temp = tars.at(lastIndexTar);
			tars.at(lastIndexTar) = tars.at(randDimer);
			tars.at(randDimer) = temp;
			
			lastIndexTar--;
		}
		
		//Next do tsrs
		for(int k=0; k<tsrPerCluster; k++)
		{
			//1) randomly select a dimer that has not been selected	
			randDimer = NFutil::RANDOM_INT(0,lastIndexTsr+1);
			//cout<<"   -selecting: " << randDimer << " / " << lastIndexTsr;
			//cout<<"  which is molecule: " << tsrs.at(randDimer)->getUniqueID() << endl;
			
			//2) put that dimer into a list
			currentClusterMembers.push_back(tsrs.at(randDimer));
			
			// swap, so that the random element moves to the end and we don't select it again
			Molecule * temp = tsrs.at(lastIndexTsr);
			tsrs.at(lastIndexTsr) = tsrs.at(randDimer);
			tsrs.at(randDimer) = temp;
			
			lastIndexTsr--;
		}
		
		//Now, create the hex grid if we have to
		if(p.get_useNeighborRxns())
			createRandomHexNeighborhood(currentClusterMembers, bondSeq, 0);
			
		//Now, go through that list of receptors we selected and put them into a new group
		
		//Create a new group that will hold our AN_Cluster
		//Create a new group that will hold our AN_Cluster
		new RecCluster(CLUSTER_NAME,s,methylStateIndex,
				p.get_clusterTarAspKoff(),p.get_clusterTarAspKon(),
				p.get_clusterTsrAspKoff(),p.get_clusterTsrAspKon(),p.get_aspartateConcentration());
		
		/*new RecClusterWill(CLUSTER_NAME,s,methylStateIndex,p);*/
		
		//Put each molecule into that cluster
		for(molIter=currentClusterMembers.begin(); molIter!=currentClusterMembers.end(); molIter++)
			s->getGroup(s->getGroupCount()-1)->addToGroup((*molIter));
	
		//Clear the currentClusterMembers vector
		currentClusterMembers.clear();
	}
	return clusterCount;
	
	
	
}




/*
 * Creates as many clusters as possible of the same size.  The members of each cluster are
 * chosen in the proper ratio of Tar to Tsr.  Thus, if the overall ratio of Tar to Tsr is 1:2, then
 * in each cluster this ratio will be exactly preserved.
 * 
 * Through the 'withAN' parameter, you can either create the clusters with the hex grid bond
 * pattern, or, if set to false, will only add them to clusters.
 * 
 * This function returns the number of clusters created, or -1 if there was some
 * kind of error (most likely you cannot evenly divide the receptor dimers into clusters
 * with the given size)
 */
int ReceptorCluster::createFixedSizeClustersWithConstantRatio(vector <Molecule *> &receptorDimers, 
	int size, 
	bool withAN,
	System *s,
	double Asp_Koff_TAR,
	double Asp_Kon_TAR,
	double Asp_Koff_TSR,
	double Asp_Kon_TSR,
	double Asp_Conc )
{
	//first create two new vectors to contain either Tar or Tsr
	vector <Molecule *>::iterator molIter;
	vector <Molecule *> tars;
	vector <Molecule *> tsrs;
	for(molIter=receptorDimers.begin(); molIter!=receptorDimers.end(); molIter++)
	{
		int ReceptorType = (*molIter)->getState((*molIter)->getMoleculeType()->getStateIndex("type"));
		if(ReceptorType==TAR)
			tars.push_back((*molIter));
		else if(ReceptorType==TSR)
			tsrs.push_back((*molIter));
		else
		{
			cerr<<"You are trying to create a cluster (in file AN_clusters.cpp) that has a ";
			cerr<<"receptor that is not a TAR or TSR!!!"<<endl;
			cerr<<"What were you thinking!?!?"<<endl;
			exit(1);	
		}	
	}
	
	cout<<endl<<"Creating clusters of fixed size and fixed composition of size: " << size<<endl;
	cout<<"Number of Tar Dimers: "<<tars.size()<<endl;
	cout<<"Number of Tsr Dimers: "<<tsrs.size()<<endl;
	cout<<"Looks like you can create "<<receptorDimers.size()/size<<"  clusters";
	cout<<" with " << receptorDimers.size()%size << " leftover"<<endl;
	cout<<"There will be " << tars.size()/(receptorDimers.size()/size) <<" Tar Dimers per cluster and"<<endl;
	cout<<"There will be " << tsrs.size()/(receptorDimers.size()/size) <<" Tsr Dimers per cluster"<<endl;
	
	//now, check if we can create even clusters, if not, don't continue
	if(receptorDimers.size()%size!=0)
	{
		cerr<<"!!  You are trying to create receptor clusters of size "<<size<<", but you have ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		cerr<<"!!  That gives you "<< receptorDimers.size()%size <<" leftover! "<<endl;
		exit(1);
	}
	if(tars.size()%(receptorDimers.size()/size)!=0)
	{
		cerr<<"!!  You are trying to create "<<receptorDimers.size()/size<<" clusters, but you have ";
		cerr<<tars.size()<<" tar receptors to go around evenly!! " <<endl;
		cerr<<"!!  That gives you "<< tars.size()%(receptorDimers.size()/size) <<" tar leftover! "<<endl;
		exit(1);
	}
	if(tsrs.size()%(receptorDimers.size()/size)!=0)
	{
		cerr<<"!!  You are trying to create "<<receptorDimers.size()/size<<" clusters, but you have ";
		cerr<<tsrs.size()<<" tsr receptors to go around evenly!! " <<endl;
		cerr<<"!!  That gives you "<< tsrs.size()%(receptorDimers.size()/size) <<" tsr leftover! "<<endl;
		exit(1);
	}
	if(receptorDimers.size()<=0)
	{
		cerr<<"!!  Can't create a receptor cluster with no receptors! You gave me ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		exit(1);
	}
	
	//Some initial variables to get us started
	int clusterCount = receptorDimers.size()/size;
	int lastIndexTar = tars.size()-1;
	int lastIndexTsr = tsrs.size()-1;
	int randDimer = -1;
	int methylStateIndex = receptorDimers.at(0)->getMoleculeType()->getStateIndex("m");
	int tarPerCluster = tars.size()/(receptorDimers.size()/size);
	int tsrPerCluster = tsrs.size()/(receptorDimers.size()/size);
	
	vector <Molecule *> currentClusterMembers;
	//vector <Molecule *>::iterator molIter;
	vector <int> bondSeq;
	
	//We have to create a bond sequence if we want to create the hexgrid
	if(withAN)
		ReceptorCluster::createStandardHexBondTraversalSeq(bondSeq,0, size);
	
	//Populate all the clusters
	for(int c=0; c<clusterCount; c++)
	{
		//cout<<endl<<"populating cluster " << c+1 << endl;
		
		//First do tars
		for(int k=0; k<tarPerCluster; k++)
		{
			//1) randomly select a dimer that has not been selected	
			randDimer = NFutil::RANDOM_INT(0,lastIndexTar+1);
			//cout<<"   -selecting: " << randDimer << " / " << lastIndexTar;
			//cout<<"  which is molecule: " << tars.at(randDimer)->getUniqueID() << endl;
			
			//2) put that dimer into a list
			currentClusterMembers.push_back(tars.at(randDimer));
			
			// swap, so that the random element moves to the end and we don't select it again
			Molecule * temp = tars.at(lastIndexTar);
			tars.at(lastIndexTar) = tars.at(randDimer);
			tars.at(randDimer) = temp;
			
			lastIndexTar--;
		}
		
		//Next do tsrs
		for(int k=0; k<tsrPerCluster; k++)
		{
			//1) randomly select a dimer that has not been selected	
			randDimer = NFutil::RANDOM_INT(0,lastIndexTsr+1);
			//cout<<"   -selecting: " << randDimer << " / " << lastIndexTsr;
			//cout<<"  which is molecule: " << tsrs.at(randDimer)->getUniqueID() << endl;
			
			//2) put that dimer into a list
			currentClusterMembers.push_back(tsrs.at(randDimer));
			
			// swap, so that the random element moves to the end and we don't select it again
			Molecule * temp = tsrs.at(lastIndexTsr);
			tsrs.at(lastIndexTsr) = tsrs.at(randDimer);
			tsrs.at(randDimer) = temp;
			
			lastIndexTsr--;
		}
		
		//Now, create the hex grid if we have to
		if(withAN)
			createRandomHexNeighborhood(currentClusterMembers, bondSeq, 0);
			
		//Now, go through that list of receptors we selected and put them into a new group
		
		//Create a new group that will hold our AN_Cluster
		new RecCluster(CLUSTER_NAME,s,methylStateIndex,
				Asp_Koff_TAR,Asp_Kon_TAR,Asp_Koff_TSR,Asp_Kon_TSR,Asp_Conc);
		
		//Put each molecule into that cluster
		for(molIter=currentClusterMembers.begin(); molIter!=currentClusterMembers.end(); molIter++)
			s->getGroup(s->getGroupCount()-1)->addToGroup((*molIter));
	
		//Clear the currentClusterMembers vector
		currentClusterMembers.clear();
	}
	return clusterCount;
}


/*
 * Creates as many clusters as possible of the same size.  The members of each cluster are
 * randomly chosen from the receptorDimers vector, so 
 * 
 * Through the 'withAN' parameter, you can either create the clusters with the hex grid bond
 * pattern, or, if set to false, will only add them to clusters.
 * 
 * This function returns the number of clusters created, or -1 if there was some
 * kind of error (most likely you cannot evenly divide the receptor dimers into clusters
 * with the given size)
 */
int ReceptorCluster::createRandomFixedSizeClusters(vector <Molecule *> &receptorDimers, 
	int size, 
	bool withAN,
	System *s,
	double Asp_Koff_TAR,
	double Asp_Kon_TAR,
	double Asp_Koff_TSR,
	double Asp_Kon_TSR,
	double Asp_Conc )
{
	cout<<endl<<"Creating clusters of fixed size " << size<<endl;
	cout<<"Number of Receptor Dimers: "<<receptorDimers.size()<<endl;
	cout<<"Looks like you can create "<<receptorDimers.size()/size<<"  clusters";
	cout<<" with " << receptorDimers.size()%size << " leftover"<<endl;
	
	//First, check if we can create even clusters, if not, don't continue
	if(receptorDimers.size()%size!=0)
	{
		cerr<<"!!  You are trying to create receptor clusters of size "<<size<<", but you have ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		cerr<<"!!  That gives you "<< receptorDimers.size()%size <<" leftover! "<<endl;
		exit(1);
	}
	if(receptorDimers.size()<=0)
	{
		cerr<<"!!  Can't create a receptor cluster with no receptors! You gave me ";
		cerr<<receptorDimers.size()<<" receptors !! " <<endl;
		exit(1);
	}
	
	//Some initial variables to get us started
	int clusterCount = receptorDimers.size()/size;
	int lastIndex = receptorDimers.size()-1;
	int randDimer = -1;
	int methylStateIndex = receptorDimers.at(0)->getMoleculeType()->getStateIndex("m");
	
	vector <Molecule *> currentClusterMembers;
	vector <Molecule *>::iterator molIter;
	vector <int> bondSeq;
	
	//We have to create a bond sequence if we want to create the hexgrid
	if(withAN)
		createStandardHexBondTraversalSeq(bondSeq,0, size);
	
	//Populate all the clusters
	for(int c=0; c<clusterCount; c++)
	{
		//cout<<endl<<"populating cluster " << c+1 << endl;
		for(int k=0; k<size; k++)
		{
			//1) randomly select a dimer that has not been selected	
			randDimer = NFutil::RANDOM_INT(0,lastIndex+1);
			//cout<<"   -selecting: " << randDimer << " / " << lastIndex;
			//cout<<"  which is molecule: " << receptorDimers.at(randDimer)->getUniqueID() << endl;
			
			//2) put that dimer into a list
			currentClusterMembers.push_back(receptorDimers.at(randDimer));
			
			// swap, so that the random element moves to the end and we don't select it again
			Molecule * temp = receptorDimers.at(lastIndex);
			receptorDimers.at(lastIndex) = receptorDimers.at(randDimer);
			receptorDimers.at(randDimer) = temp;
			
			lastIndex--;
		}
		
		//Now, create the hex grid if we have to
		if(withAN)
			createRandomHexNeighborhood(currentClusterMembers, bondSeq, 0);
			
		//Now, go through that list of receptors we selected and put them into a new group
		
		//Create a new group that will hold our AN_Cluster
		new RecCluster(CLUSTER_NAME,s,methylStateIndex,
				Asp_Koff_TAR,Asp_Kon_TAR,Asp_Koff_TSR,Asp_Kon_TSR,Asp_Conc);
		
		//Put each molecule into that cluster
		for(molIter=currentClusterMembers.begin(); molIter!=currentClusterMembers.end(); molIter++)
			s->getGroup(s->getGroupCount()-1)->addToGroup((*molIter));
	
		//Clear the currentClusterMembers vector
		currentClusterMembers.clear();
	}
	return clusterCount;
}




void ReceptorCluster::distributeRandomCheAtoClusters_Evenly(System * s, vector <Molecule *> &cheA, int numOfClusters)
{
	int n_groups = s->getGroupCount();
	
	//Check if we can evenly distribute CheA amongst all the clusters
	if(cheA.size()%numOfClusters!=0)
	{
		cerr<<"!!  I cannot distribute "<<cheA.size()<<" molecules evenly over  ";
		cerr<<numOfClusters<<" clusters !! " <<endl;
		exit(1);
	}
	
	int cheAPerCluster = cheA.size()/numOfClusters;
	//cout<<"Putting "<<cheAPerCluster<<" CheA molecules in each cluster."<<endl;
	int gCounter = 0;
	int lastIndex = cheA.size()-1;
	int randCheA = -1;
	
	for(int g=0; g<n_groups; g++)
	{
		Group *group = s->getGroup(g);
		//if(strlen(group->getName())==strlen(CLUSTER_NAME))
		if(group->getName()==CLUSTER_NAME)
		{
			for(int k=0; k<cheAPerCluster; k++)
			{
				//1) randomly select a dimer that has not been selected	
				randCheA = NFutil::RANDOM_INT(0,lastIndex+1);
				//cout<<"   -selecting cheA: " << randCheA << " / " << lastIndex;
				//cout<<"  which is molecule: " << cheA.at(randCheA)->getUniqueID() << endl;
		
				//2) Add that cheA to the cluster
				group->addToGroupWithoutListener(cheA.at(randCheA));
		
				// swap, so that the random element moves to the end and we don't select it again
				Molecule * temp = cheA.at(lastIndex);
				cheA.at(lastIndex) = cheA.at(randCheA);
				cheA.at(randCheA) = temp;
		
				lastIndex--;
			}
			gCounter++;
		}
			
	}
	
	if(gCounter!=numOfClusters)
	{
		cerr<<"There were actually "<<gCounter<<" clusters, but you told me there were "<< numOfClusters;
		cerr<<"!!!!"<<endl<<"You lied to me!!  I can't add CheA to clusters in conditions like this!"<<endl;
		exit(1);
	}
}








/*
 * Creates a given number of clusters, each with different sizes.  Each receptor has an equal
 * probability of being added to any cluster, thus there is a simple poisson distribution of
 * cluster sizes.
 * 
 * Through the 'withAN' parameter, you can either create the clusters with the hex grid bond
 * pattern, or, if set to false, will only add them to clusters.
 */
int ReceptorCluster::createRandomVariableSizeClusters(vector <Molecule *> &receptorDimers, int numOfClusters, bool withAN)
{
	
	return 0;
}





