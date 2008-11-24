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




DumpMoleculeType::DumpMoleculeType(string filename, System *s, MoleculeType *mt) : Outputter(filename,s) {
	this->mt=mt;

	//open up the fileStream
	outputFileStream.open(filename.c_str());
}
DumpMoleculeType::~DumpMoleculeType() {
	if(outputFileStream.is_open()) {
		outputFileStream.flush();
		outputFileStream.close();
	}
}





void DumpMoleculeType::outputHeader() {
	outputFileStream<<"#\tTime";
	outputFileStream<<"\tMolecule("<<mt->getName()<<")";
	outputFileStream<<"\tComplex";
	for(int c=0; c<mt->getNumOfComponents(); c++ ) {
		outputFileStream<<"\t"<<mt->getComponentName(c)<<"(state)";
		outputFileStream<<"\t"<<mt->getComponentName(c)<<"(bond)";
	}
	outputFileStream<<endl;
}









/*
 * A version of the  breadth first search that looks for only a single
 * moleculeType, and only continues the search if the moleculeType is
 * connected only through other moleculeTypes.
 */
void breadthFirstSearchSeqStyle(list <Molecule *> &members, Molecule *m)
{
	MoleculeType *mt = m->getMoleculeType();
	//Create the queues (for efficiency, now queues are a static attribute of Molecule...)
	queue <Molecule *> q;
	queue <int> d;
	list <Molecule *>::iterator molIter;
	int currentDepth = 0;

	//First add this molecule
	q.push(m);
	members.push_back(m);
	d.push(currentDepth+1);
	m->hasVisitedMolecule=true;

	//Look at children until the queue is empty
	while(!q.empty())
	{
		//Get the next parent to look at (currentMolecule)
		Molecule *cM = q.front();
		currentDepth = d.front();
		q.pop();
		d.pop();


		//Loop through the bonds
		int cMax = cM->getMoleculeType()->getNumOfComponents();
		for(int c=0; c<cMax; c++)
		{
			if(cM->isBindingSiteBonded(c))
			{
				Molecule *neighbor = cM->getBondedMolecule(c);
				if(mt==neighbor->getMoleculeType()) {
					if(!neighbor->hasVisitedMolecule)
					{
						neighbor->hasVisitedMolecule=true;
						members.push_back(neighbor);
						q.push(neighbor);
						d.push(currentDepth+1);
					}
				}
			}
		}
	}


	//clear the has visitedMolecule values
	for( molIter = members.begin(); molIter != members.end(); molIter++ )
  		(*molIter)->hasVisitedMolecule=false;
}





void DumpMoleculeType::output() {

	if(s->isUsingComplex()) {
		cout<<"output function for complexes not working yet."<<endl;


	} else {

		int currentComplex = 0;
		list <Molecule *> molecules;
		list <Molecule *>::iterator molIter;
		for(int i=0; i<mt->getMoleculeCount(); i++) {
			Molecule *m = mt->getMolecule(i);

			//If the complexID is less than the current complex, than this Molecule
			//has already been dumped and was already identified.
			//cout<<"currentComplex: " <<currentComplex<<endl;
			//cout<<"looking at molecule: "<<m->getUniqueID()<<" cId: "<<m->getComplexID()<<endl;
			if(m->getComplexID()<currentComplex && m->getComplexID()>=0) { continue; }
			//cout<<"assembling complex"<<endl;

			//If we have not seen this molecule yet, traverse the neighborhood
			//and grab all the molecules that are connected to this molecule
			breadthFirstSearchSeqStyle(molecules, m);



			//Make sure it is of the correct type, because we are only outputting
			//molecules of a certain type

			//cout<<"starting the list: "<<molecules.size()<<"\n";

			for(molIter=molecules.begin(); molIter!=molecules.end(); molIter++) {


				//Make sure we mark it in this complex so we don't check it again
				(*molIter)->setComplexID(currentComplex);

				//Only output if it is of the correct type
				if((*molIter)->getMoleculeType()==mt) {
					//cout<<"\t"<<(*molIter)->getUniqueID()<<endl;
					outputFileStream<<"\t"<<s->getCurrentTime();
					//outputFileStream<<"\t"<<(*molIter)->getMoleculeTypeName()<<"_"<<(*molIter)->getUniqueID();
					outputFileStream<<"\t"<<(*molIter)->getUniqueID();
					outputFileStream<<"\t"<<currentComplex;
					for(int c=0; c<mt->getNumOfComponents(); c++ ) {
						//outputFileStream<<"\t"<<mt->getComponentStateName(c,(*molIter)->getComponentState(c));
						outputFileStream<<"\t"<<(*molIter)->getComponentState(c);
						if((*molIter)->isBindingSiteBonded(c)) {
							Molecule *m2 = (*molIter)->getBondedMolecule(c);
							//outputFileStream<<"\t"<<m2->getMoleculeTypeName()<<"_"<<m2->getUniqueID();
							outputFileStream<<"\t"<<m2->getUniqueID();
						} else {
							outputFileStream<<"\t-1";
						}
					}
					outputFileStream<<"\n";
				}


			} //end loop over bonded neighborhood
			molecules.clear();
		//	outputFileStream<<"finished the list\n";


			currentComplex++;
		} //end loop over all molecules

		//Reset complex ids... may be able to eliminate this loop, but for
		//now it is safer
		for(int i=0; i<mt->getMoleculeCount(); i++) {
			mt->getMolecule(i)->setComplexID(-1);
		}
	}
	outputFileStream.flush();

}





//This is called after we output all molecules if we are not
//using complex bookkeeping
void clearMoleculeComplexIds(System *s) {
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		for(int j=0; j<mt->getMoleculeCount(); j++) {
			mt->getMolecule(j)->setComplexID(-1);
		}
	}
}






DumpSystem::DumpSystem(System *s, vector <double> dumpTimes)
{
	this->s=s;
	this->currentDumpTimeIndex = 0;
	for(unsigned int i=0;i<dumpTimes.size(); i++) {
		this->dumpTimes.push_back(dumpTimes.at(i));
	}
}
//DumpSystem::~DumpSystem() {};

void DumpSystem::tryToDump(double simTime)
{

	//cout<<"Trying to dump at time: "<<simTime<<endl;

	if(currentDumpTimeIndex>=dumpTimes.size()) {
		//cout<<"but I already dumped all of my load."<<endl;exit(1);
		return;
	}

	bool dumped = false;
	while(simTime>=dumpTimes.at(currentDumpTimeIndex)) {

		dumped = true;
		cout<<"dumping at: "<<dumpTimes.at(currentDumpTimeIndex)<<endl;
		dumpHeaderFile(dumpTimes.at(currentDumpTimeIndex));
		dumpMoleculeTypeFiles(dumpTimes.at(currentDumpTimeIndex));

		if(!s->isUsingComplex()) clearMoleculeComplexIds(s);

		currentDumpTimeIndex++;
		if(currentDumpTimeIndex>=dumpTimes.size()) break;
	}
	//if(!dumped) {
	//	cout<<"no need to dump yet."<<endl;
	//}
}




string toString(double x)
{
	std::ostringstream o;
	if (!(o << x)) {
		cout<<endl; cerr<<"Error converting double to string."<<endl;
		exit(1);
	}
	return o.str();
}
string toString(int x)
{
	std::ostringstream o;
	if (!(o << x)) {
		cout<<endl; cerr<<"Error converting double to string."<<endl;
		exit(1);
	}
	return o.str();
}



void DumpSystem::dumpHeaderFile(double dumpTime) {
	string dumpFileName = s->getName()+"_nf."+toString(dumpTime)+".dump.head";
	cout<<"writing file: "<<dumpFileName<<endl;;


	ofstream ofs;
	//try {

	//ios_base::out -- Set for output only, instead of for input/output
	//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there
	ofs.open((dumpFileName).c_str(), ios_base::out | ios_base::trunc);

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
			ofs<<"\t"<<i<<"\t"<<k<<"\t"<<mt->getTypeIILocalFunction(k)->getName()<<"\n";
		}
	}

	ofs<<"\n>> EOF ##########################################\n";
	ofs.flush();
	ofs.close();
//	} catch ()





}
void DumpSystem::dumpMoleculeTypeFiles(double dumpTime) {

	double complexCount = 0;
	for(int i=0; i<s->getNumOfMoleculeTypes(); i++) {
		MoleculeType *mt = s->getMoleculeType(i);
		string dumpFileName = s->getName() + "_nf."+toString(dumpTime)+".dump."+toString(i);
		cout<<"writing file: "<<dumpFileName<<endl;

		ofstream ofs;
			//try {

		//ios_base::out -- Set for output only, instead of for input/output
		//ios_base::binary --  Set output to binary
		//ios_base::trunc --  Truncate the file - that is overwrite anything that was already there
		ofs.open((dumpFileName).c_str(), ios_base::out | ios_base::binary | ios_base::trunc);

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

		}


		ofs.flush();
		ofs.close();

	}


}
































