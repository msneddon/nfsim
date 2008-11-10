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





