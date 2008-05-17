#include "reactantList.hh"


using namespace NFcore;

ReactantList::ReactantList(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity=50)
{

	this->n_mappingSets = 0;
	this->capacity = init_capacity;
	this->reactantIndex = reactantIndex;
	this->ts=ts;
	this->mappingSets = new MappingSet *[init_capacity];
	this->msPositionMap = new unsigned int [init_capacity];
	for(int i=0; i<this->capacity; i++)
	{
		mappingSets[i]=ts->generateBlankMappingSet(reactantIndex,i);
		msPositionMap[i]=i;
	}
}

ReactantList::~ReactantList()
{
	for(int i=0; i<capacity; i++)
	{
		delete mappingSets[i];
		msPositionMap[i]=0;
	}
	delete [] mappingSets;
	delete [] msPositionMap;
	this->n_mappingSets = 0;
	this->capacity = 0;
}


void ReactantList::pickRandom(MappingSet *&ms)
{
	unsigned int rand = NFutil::RANDOM_INT(0,n_mappingSets);
	ms = mappingSets[rand];
}



int ReactantList::size() const
{
	return n_mappingSets;
}


MappingSet * ReactantList::pushNextAvailableMappingSet()
{
	//Check if we are going to exceed capacity
	if(n_mappingSets>=capacity)
	{
		//Copy everything over to new arrays that are double the size
		MappingSet ** new_mappingSets = new MappingSet *[capacity*2];
		unsigned int * new_msPositionMap = new unsigned int [capacity*2];
		for(int i=0; i<capacity; i++)  {
			new_mappingSets[i] = mappingSets[i];
			new_msPositionMap[i] = msPositionMap[i];
		}
		for(int i=capacity; i<capacity*2; i++)  {
			new_mappingSets[i] = ts->generateBlankMappingSet(reactantIndex, i);
			new_msPositionMap[i] = i;
		}
		
		//Swap the copied data with the real data and double the capacity
		delete [] mappingSets;
		delete [] msPositionMap;
		mappingSets = new_mappingSets;
		msPositionMap = new_msPositionMap;
		capacity*=2;
	}
	
	//Increase the number of reactants, and return the activated mappingSet
	n_mappingSets++;
	return mappingSets[n_mappingSets-1];
}


void ReactantList::popLastMappingSet()
{
	if(n_mappingSets<=0) {
		cerr<<"Trying to pop an empty ReactantList!!"<<endl;
		exit(1);
	}
	
	//Clear out the mappingSet (just in case) and decrease the count
	mappingSets[n_mappingSets-1]->clear();
	n_mappingSets--;
}

void ReactantList::removeMappingSet(unsigned int mappingSetId)
{
	//Make sure this mappingSet is not empty
	if(n_mappingSets==0) {
		cerr<<"Trying to remove from an empty ReactantList!!"<<endl;
		if(n_mappingSets==0) {
				cerr<<"Trying to pop an empty ReactantList!!"<<endl;
				exit(1);
			}exit(1);
	}
	
	//First, get the position of the mappingSet we need to remove
	int pos = msPositionMap[mappingSetId];
	
	//Make sure the position is valid (not out of bounds of the List)
	if(pos+1>(n_mappingSets)) {
		cout<<"Error in ReactantList:  you can't remove a mappingSet that has been cleared! (trying to remove: "<< mappingSetId << " in pos " << pos <<" but size is: "<<size()<<endl;
		printDetails();
		
		return;
	}
	
	//If the array has only one element, or we just happened to select the last element,
	//then just remove the last element without a swap
	if( pos+1 == (n_mappingSets) ) {
		popLastMappingSet();
		return;
	}
	
	//Otherwise, we have to swap with the last element in the list
	MappingSet *tempMappingSet = mappingSets[pos];
	mappingSets[pos] = mappingSets[n_mappingSets-1];
	mappingSets[n_mappingSets-1] = tempMappingSet;
	
	//Careful here!  We have to swap values in the msPositionMap
	//so that msPositionMap[mappingId] points correctly to where
	//that MappingSet now lives.  Notice that this is a different
	//swap than above in the MappingSets array!
	msPositionMap[mappingSetId] = n_mappingSets-1;
	msPositionMap[mappingSets[pos]->getId()] = pos;
	
	//Make sure we clear what we don't need
	tempMappingSet->clear();

	//Remember to remove 
	n_mappingSets--;
}


void ReactantList::printDetails()
{
	//Used for debuggin'...
	cout<<"ReactantList that contains: "<<size()<<" MappingSets and has a capacity for "<<capacity<<" total sets."<<endl;
	
	for(int i=0; i<capacity; i++)
	{
		if(i<10) cout<<" ";
		cout<<"["<<i<<"]: "<<msPositionMap[i];
		if(i<(n_mappingSets)) {
			cout<<"\t\tpos="<<i<<"(mol="<<mappingSets[i]->get(0)->getMolecule()->getUniqueID()<<") ";
			}
			
		if(i==n_mappingSets-1) cout<<"  _"<<endl;
		else cout<<endl;
	}
	cout<<endl;
}
