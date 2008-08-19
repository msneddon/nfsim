#ifndef NFMODULE_HH_
#define NFMODULE_HH_












	class Module {
		public:
			Module();
			~Module();
			
			
		
		protected:
		
			//members
			MoleculeType *moleculeTypes;
			Molecule **molecules;
			Module *modules;
			
			
			//observables
			
			//functions
			unsigned int n_functions;
			string *functionNames;
			mu::Parser *p; 
			
			
			
			
			//identification of this module
			int uniqueModuleID;
			int moduleType;
			string moduleTypeName;
			
			//get at the parent module
			int parentModuleID;
			Module *parentModule;
			
			
			static unsigned int MODULE_COUNT;
	};



Module::MODULE_COUNT = 0;




















#endif /*NFMODULE_HH_*/
