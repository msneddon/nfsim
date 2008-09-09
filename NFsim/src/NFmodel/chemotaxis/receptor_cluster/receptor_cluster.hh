//#ifndef RECEPTOR_CLUSTER_HH_
//#define RECEPTOR_CLUSTER_HH_
//
//
//#include "../chemotaxis.hh"
//#include "../ng_param/ng_param.hh"
//
//
//using namespace NFcore;
//
//namespace ReceptorCluster
//{
//
//
//	void createRandomHexNeighborhood(vector <Molecule *> & molecules, vector <int> bondSeq, int firstBsiteIndex);
//	void createStandardHexBondTraversalSeq(vector <int> &bondSeq, int firstBsiteIndex, unsigned int length);
//
//	int createRandomFixedSizeClusters(vector <Molecule *> &tarDimers,
//			vector <Molecule *> &tsrDimers, int tarPerCluster, int tsrPerCluster,bool withAN);
//
//	void distributeRandomCheAtoClusters_Evenly(System * s, vector <Molecule *> &cheA, int numOfClusters);
//	int createFixedSizeClustersWithConstantRatio(vector <Molecule *> &receptorDimers, 
//		int size, 
//		bool withAN,
//		System *s,
//		double Asp_Koff_TAR,
//		double Asp_Kon_TAR,
//		double Asp_Koff_TSR,
//		double Asp_Kon_TSR,
//		double Asp_Conc );
//	
//	/** New Version!!! Use this from now on! */
//	int createFixedSizeClustersWithConstantRatio(System *s, vector <Molecule *> &receptorDimers, NGparam &p);
//	/*		int size, 
//			bool withAN,
//			System *s,
//			double Asp_Koff_TAR,
//			double Asp_Kon_TAR,
//			double Asp_Koff_TSR,
//			double Asp_Kon_TSR,
//			double Asp_Conc ); */
//	
//	
//	int createRandomFixedSizeClusters(vector <Molecule *> &receptorDimers, 
//		int size, 
//		bool withAN,
//		System *s,
//		double Asp_Koff_TAR,
//		double Asp_Kon_TAR,
//		double Asp_Koff_TSR,
//		double Asp_Kon_TSR,
//		double Asp_Conc );
//	int createRandomVariableSizeClusters(vector <Molecule *> &receptorDimers, int numOfClusters, bool withAN);
//
//	
//	
//	///////////////////////////////////////////////////////////////////////
//		// Receptor Cluster group definition
//		class RecCluster:public Group
//		{
//			public:
//			
//				RecCluster(char *name, System *s, int methStateIndex,
//					double Asp_Koff_TAR,
//					double Asp_Kon_TAR,
//					double Asp_Koff_TSR,
//					double Asp_Kon_TSR,
//					double Asp_Conc
//				);
//				
//				/** New version!! use me from now on!! */
//				RecCluster(const char *name, System *s, int methStateIndex, NGparam &p);
//				
//				
//				virtual ~RecCluster();
//				
//				virtual void addToGroup(Molecule * m);
//				virtual void notify(Molecule *changedMolecule, int oldStateValue, int newStateValue);
//				
//				/* function that lets us update ligand concentration for this group */
//				virtual void updateGroupProperty(double * values, int n_values);
//				virtual void printDetails();
//				
//			protected:
//				//Free energy offset values
//				void initFreeEnergyOffset();
//				double * TAR_freeEnergyOffset;
//				double * TSR_freeEnergyOffset;
//				
//				//Counts of what's in this cluster
//				int TAR_COUNT;
//				int TSR_COUNT;
//				int CHE_A_COUNT;
//				
//				//Kon and Koff values
//				double Asp_Koff_TAR;
//				double Asp_Kon_TAR;
//				double Asp_Koff_TSR;
//				double Asp_Kon_TSR;
//				
//				//Ligand Concentration
//				double AspConc;
//				
//				//Internal state values to remember and functions to calculate them
//				double TAR_LOG_TERM;
//				double TSR_LOG_TERM;
//				void updateTarLogValue();
//				void updateTsrLogValue();
//				
//				double FreeEnergyOffsetSum;
//		};
//
//
//		
//		class RecClusterWill:public Group
//		{
//			public:
//				
//				/** New version!! use me from now on!! */
//				RecClusterWill(string name, System *s, int methStateIndex, NGparam &p);
//				virtual ~RecClusterWill();
//				
//				
//				virtual void addToGroup(Molecule * m);
//				virtual void notify(Molecule *changedMolecule, int oldStateValue, int newStateValue);
//				
//				/* function that lets us update ligand concentration for this group */
//				virtual void updateGroupProperty(double * values, int n_values);
//				virtual void printDetails();
//				
//			protected:
//				
//				double getPon();
//				double getFofL();
//				
//				double e0;
//				double e1;
//				int methylationLevel;
//				
//				
//				//Counts of what's in this cluster
//				int TAR_COUNT;
//				int TSR_COUNT;
//				int CHE_A_COUNT;
//				
//				//Kon and Koff values
//				double Asp_Koff_TAR;
//				double Asp_Kon_TAR;
//				double Asp_Koff_TSR;
//				double Asp_Kon_TSR;
//				
//				//Ligand Concentration
//				double AspConc;
//		};
//		
//}
//
//
//
//
//
//#endif /*RECEPTOR_CLUSTER_HH_*/
