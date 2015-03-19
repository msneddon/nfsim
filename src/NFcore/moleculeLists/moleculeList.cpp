#include "moleculeList.hh"


using namespace NFcore;

MoleculeList::MoleculeList(MoleculeType *mt, int init_capacity, int finalCapacity)
{
	this->n_molecules = 0;
	this->lastAllocated = 0;
	this->capacity = init_capacity;
	this->finalCapacity=finalCapacity;
	this->mt = mt;

	this->molPos = new int [init_capacity];
	this->mArray = new Molecule * [init_capacity];

	for(int i=0; i<this->capacity; i++)
	{
		mArray[i]=new Molecule (mt,i);
		molPos[i]=i;
	}
}

MoleculeList::~MoleculeList()
{
	for(int i=0; i<capacity; i++)
	{
		delete mArray[i];
		molPos[i]=0;
	}
	delete [] mArray;
	delete [] molPos;
	this->n_molecules = 0;
	this->capacity = 0;
}


Molecule *MoleculeList::at(int index) const
{
	return mArray[index];
}

int MoleculeList::size() const
{
	return n_molecules;
}

int MoleculeList::create(Molecule *&m)
{
	//Check if we are going to exceed capacity
	if(n_molecules>=capacity)
	{
		int newCapacity;
		if(capacity>400000) {
			newCapacity = capacity+50000;
		} else {
			newCapacity=capacity*2;
		}

		if(capacity>finalCapacity && finalCapacity!=MoleculeList::NO_LIMIT) {
			cout.flush();
			cerr<<"\n\nError in Simulation!  Creating space for "<<capacity;
			cerr<<" copies of a MoleculeType: '"<<mt->getName()<<"'.\n\n";
			cerr<<"There is currently an imposed limit of: "<<finalCapacity<< " molecules \nper MoleculeType. ";
			cerr<<"This is done to keep your operating system \nfrom crashing, due to excessive system size.";
			cerr<<"  If you need \nto have more molecules, rerun with the -gml [int] flag \nto increase the limit.";
			cerr<<"  For instance, to increase the limit \nto 1 million, write: -gml 1000000.\n\n";
			cerr<<"Better luck next time!"<<endl;
			exit(1);
		}

		//Copy everything over to new arrays that are double the size
		Molecule ** new_mArray = new Molecule *[newCapacity];
		int * new_molPos = new int [newCapacity];
		for(int i=0; i<capacity; i++)  {
			new_mArray[i] = mArray[i];
			new_molPos[i] = molPos[i];
		}
		for(int i=capacity; i<newCapacity; i++)  {
			new_mArray[i] = new Molecule (mt,i);
			new_molPos[i] = i;
		}

		//Swap the copied data with the real data and double the capacity
		delete [] mArray;
		delete [] molPos;
		mArray = new_mArray;
		molPos = new_molPos;
		capacity=newCapacity;
	}

	//Increase the number of reactants, and return the activated mappingSet
	n_molecules++;
	m = mArray[n_molecules-1];

	//cout<<"ADDING!!!"<<endl;
	//printDetails();
	return m->getMolListId();
}


void MoleculeList::remove(int listId, Molecule *m)
{
	// I think this is redundant (see below).  --Justin
	//Make sure this mappingSet is not empty
	//if(n_molecules==0) {
	//	cerr<<"Trying to remove from an empty MoleculeList!!"<<endl;
	//	exit(1);
	//}

	//First, get the position of the mappingSet we need to remove
	int pos = molPos[listId];

	//Make sure the position is valid (not out of bounds of the List)
	if(pos+1>(n_molecules)) {
		// Handle this graciously, but don't abort!
		cout << "!! Warning in MoleculeList: attempt to remove a dead molecule!\n"
		     << "   This may occur when complex bookkeeping is disabled and two reactant patterns\n"
		     << "   with delete transforms match the same complex. Enable complex bookkeeping (-cb)\n"
		     << "   and see if this message disappears.\n"
		     << "   (trying to remove: " << listId << " in pos " << pos <<" but size is: " << size() << endl;
		//cout<<"Error in MoleculeList:  you can't remove a molecule that is not in the simulation! (trying to remove: "<< listId << " in pos " << pos <<" but size is: "<<size()<<endl;
		//m->printDetails();
		//exit(1);
		return;
	}

	//If the array has only one element, or we just happened to select the last element,
	//then just remove the last element without a swap
	if( pos+1 == (n_molecules) ) {
		removeLast();
		return;
	}

	//Otherwise, we have to swap with the last element in the list
	Molecule *tempMol = mArray[pos];
	mArray[pos] = mArray[n_molecules-1];
	mArray[n_molecules-1] = tempMol;

	//Careful here!
	molPos[listId] = n_molecules-1;
	molPos[mArray[pos]->getMolListId()] = pos;


	//Remember to remove
	n_molecules--;


	//cout<<"REMOVING!!! "<< listId<<" at pos: "<<pos<<endl;
	//printDetails();
}

void MoleculeList::removeLast()
{
	if(n_molecules<=0) {
		cerr<<"Trying to remove the last element from an empty MoleculeList!!"<<endl;
		exit(1);
	}

	n_molecules--;
}



void MoleculeList::printDetails()
{
	//Used for debuggin'...
	cout<<"ReactantList that contains: "<<size()<<" MappingSets and has a capacity for "<<capacity<<" total sets."<<endl;

	for(int i=0; i<capacity; i++)
	{
		if(i<10) cout<<" ";
		cout<<"["<<i<<"]: "<<molPos[i];
		if(i<(n_molecules)) {
			cout<<"\t\tpos="<<i<<"(mol="<<mArray[i]->getMolListId()<<") ";
		}

		if(i==n_molecules-1) cout<<"  _"<<endl;
		else cout<<endl;
	}
	cout<<endl;
}
