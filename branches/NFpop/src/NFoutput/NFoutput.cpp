#include "NFoutput.hh"

#include <queue>
#include <list>

using namespace NFcore;

Outputter::Outputter(string filename, System *s) {
	this->filename = filename;
	this->s = s;
}
Outputter::~Outputter() {
	if(outputFileStream.is_open()) {
		outputFileStream.flush();
		outputFileStream.close();
	}
}





//! Careful!  potentially dangerous function.  use with care.
/*!
 */
void clearMoleculeComplexIds(System *s) {
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		for(int j=0; j<mt->getMoleculeCount(); j++) {
			mt->getMolecule(j)->setComplexID(-1);
		}
	}
}






DumpSystem::DumpSystem(System *s, vector <double> dumpTimes, string pathToFolder, bool verbose)
{
	this->s=s;
	this->currentDumpTimeIndex = 0;
	for(unsigned int i=0;i<dumpTimes.size(); i++) {
		this->dumpTimes.push_back(dumpTimes.at(i));
	}
	this->verbose = verbose;
	this->pathToFolder=pathToFolder;
}
DumpSystem::~DumpSystem() {};

void DumpSystem::tryToDump(double simTime)
{
	//If we've dumped all of our load
	if(currentDumpTimeIndex>=(int)dumpTimes.size()) {
		return;
	}

	while(simTime>=dumpTimes.at(currentDumpTimeIndex)) {

		if(verbose) cout<<"dumping at: "<<dumpTimes.at(currentDumpTimeIndex)<<endl;
		dumpHeaderFile(dumpTimes.at(currentDumpTimeIndex));
		dumpMoleculeTypeFiles(dumpTimes.at(currentDumpTimeIndex));

		if(!s->isUsingComplex()) clearMoleculeComplexIds(s);

		currentDumpTimeIndex++;
		if(currentDumpTimeIndex>=(int)dumpTimes.size()) break;
	}
}








void DumpSystem::dumpHeaderFile(double dumpTime) {
	string dumpFileName = pathToFolder + s->getName()+"_nf."+NFutil::toString(dumpTime)+".dump.head";

	//cout<<"writing file: "<<dumpFileName<<endl;;

	ofstream ofs;
	//ios_base::out -- Set for output only, instead of for input/output
	//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there
	ofs.open((dumpFileName).c_str(), ios_base::out | ios_base::trunc);

	//Make sure we could open the stream (this is really the only thing that can go around)
	if(!ofs.is_open()) {
		cout<<"--------"<<endl;
		cout<<"-Error dumping header file.";
		cout<<"-Could not open stream to file: '"<<dumpFileName<<"'"<<endl;
		cout<<"-The path is probably incorrect or does not exist."<<endl;
		return;
	}

	ofs<<"## Automatically generated file from NFsim.  Do not edit manually.\n";
	ofs<<"## To read this file and the associated information, use the\n";
	ofs<<"## included Matlab script readNFdump.m.\n\n";

	ofs<<"\n>> Time #######################################\n";
	ofs<<dumpTime<<"\n";

	ofs<<"\n>> MoleculeTypeDef ############################\n";
	ofs<<"##\tindex\tname\tnumOfInstances\tnumOfComponents\n";
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		ofs<<"\t"<<i<<"\t"<<mt->getName()<<"\t"<<mt->getMoleculeCount()<<"\t";
		ofs<<mt->getNumOfComponents()<<"\t"<<mt->getNumOfTypeIFunctions()<<"\n";
	}

	ofs<<"\n>> MoleculeTypeComponents #####################\n";
	ofs<<"##\ttype\tindex\tname\n";
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		for(int k=0; k<mt->getNumOfComponents(); k++) {
			ofs<<"\t"<<i<<"\t"<<k<<"\t"<<mt->getComponentName(k)<<"\n";
		}
	}

	ofs<<"\n>> MoleculeTypeFunctions ######################\n";
	ofs<<"##\ttype\tindex\tname\n";
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		for(int k=0; k<mt->getNumOfTypeIFunctions(); k++) {
			ofs<<"\t"<<i<<"\t"<<k<<"\t"<<mt->getTypeILocalFunction(k)->getName()<<"\n";
		}
	}

	ofs<<"\n>> EOF ##########################################\n";
	ofs.flush();
	ofs.close();
}
void DumpSystem::dumpMoleculeTypeFiles(double dumpTime) {

	double complexCount = 0;
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		string dumpFileName = pathToFolder+s->getName() + "_nf."+NFutil::toString(dumpTime)+".dump."+NFutil::toString(i);

		ofstream ofs;

		//ios_base::out -- Set for output only, instead of for input/output
		//ios_base::binary --  Set output to binary
		//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there
		ofs.open((dumpFileName).c_str(), ios_base::out | ios_base::binary | ios_base::trunc);

		//Make sure we could open the stream (this is really the only thing that can go around)
		if(!ofs.is_open()) {
			cout<<"--------"<<endl;
			cout<<"-Error dumping header file.";
			cout<<"-Could not open stream to file: '"<<dumpFileName<<"'"<<endl;
			cout<<"-The path is probably incorrect or does not exist."<<endl;
			return;
		}


		double NOBOND = -1;
		list <Molecule *> molList; list <Molecule *>::iterator molIter;
		for(int j=0; j<mt->getMoleculeCount(); j++) {

			//First, write out this molecules unique ID
			Molecule *m=mt->getMolecule(j);
			double uId = (double)m->getUniqueID();
			ofs.write((char *)&uId,sizeof(double));

			//Second, output the complex number so we can quickly recover
			//the complexes
			if(!s->isUsingComplex()) {
				if(m->getComplexID()==-1) { //if we haven't visited this guy yet...
					m->traverseBondedNeighborhood(molList,ReactionClass::NO_LIMIT);
					for(molIter=molList.begin(); molIter!=molList.end(); molIter++) {
						(*molIter)->setComplexID((int)complexCount);
					}
					ofs.write((char *)&complexCount,sizeof(double));
					complexCount++;
					molList.clear();
				} else {
					double thisUid = m->getComplexID();
					ofs.write((char *)&thisUid,sizeof(double));
				}
			} else {
				complexCount = m->getComplexID();
				ofs.write((char *)&complexCount,sizeof(double));
			}



			//Now export all the components, marking its state and
			//binding partner if any
			for(int k=0; k<mt->getNumOfComponents(); k++) {
				double stateValue = (double)m->getComponentState(k);
				ofs.write((char *)&stateValue,sizeof(double));
				if(!m->isBindingSiteBonded(k)) {
					ofs.write((char *)&NOBOND,sizeof(double));
				} else {
					Molecule *m2=m->getBondedMolecule(k);
					double uId2 = (double)m2->getUniqueID();
					ofs.write((char *)&uId2,sizeof(double));
				}
			}

			for(int k=0; k<mt->getNumOfTypeIFunctions(); k++) {
				double val = m->getLocalFunctionValue(k);
				LocalFunction *lf = m->getLocalFunction(k);
				double localValue = lf->evaluateOn(m,LocalFunction::MOLECULE);
				ofs.write((char *)&val,sizeof(double));
				ofs.write((char *)&localValue,sizeof(double));
			//		ofs<<"\t"<<i<<"\t"<<k<<"\t"<<mt->getTypeIILocalFunction(k)->getName()<<"\n";
			}



		}


		ofs.flush();
		ofs.close();

	}
}

