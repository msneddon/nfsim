


#include <math.h>
#include "reactantTree.hh"




using namespace NFcore;

ReactantTree::ReactantTree(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity)
{
	//cout<<"Creating Reactant Tree..."<<endl;
	//cout<<"  Needs to be able to store "<< maxElementCount <<" molecules, ";

	//first, figure out the initial capacity we would like to have
	//and start out with zero mapping sets
	this->n_mappingSets = 0;
	this->capacity = init_capacity;
	if(capacity<4) capacity=4;
	
	//For a binary tree, the depth is the rounded up value of the log base 10 of
	//the capacity divided by the log base 10 of 2.
	this->treeDepth = (unsigned int)ceil(log(capacity)/log(2));
	
	//We have a complete binary tree, so we must recalculate the capacity
	//from the depth of the tree (this is the number of nodes we can store
	this->capacity = (unsigned int) pow(2,treeDepth);
	
	//Calculate the number of nodes we need (and thus the length of the tree arrays)
	this->numOfNodes = ((unsigned int) pow(2,treeDepth+1)) - 1;
	
	//Calculate the first index (into the tree! not the mappingset array!) that can store a molecule element
	//This is, conveniently, the capacity of the tree (gotta love binary!)
	this->firstMoleculeTreeIndex = this->capacity;
	
	//Set up and init our tree arrays
	this->leftRateFactorSum = new double [numOfNodes+1];
	this->leftElementCount = new unsigned int [numOfNodes+1];
	this->rightElementCount = new unsigned int [numOfNodes+1];
	for(unsigned int i=0; i<=numOfNodes; i++) {
		this->leftRateFactorSum[i] = 0;
		this->leftElementCount[i] = 0;
		this->rightElementCount[i] = 0;
	}
	
	//Remember to keep track of the TransformationSet information
	this->reactantIndex = reactantIndex;
	this->ts=ts;
	
	//Set up the mappingSets array to store our mappingSets by using
	//the transformation set to generate blank mapping set objects
	this->mappingSets = new MappingSet * [this->capacity];
	for(unsigned int i=0; i<this->capacity; i++) {
		this->mappingSets[i] = ts->generateBlankMappingSet(reactantIndex,i);
	}
	
	
	//cout<<"so setting a limit of " << this->capacity <<" molecules. "<<endl;
	//cout<<"  The depth of the tree will be "<< treeDepth << " and contain ";
	//cout<<numOfNodes<<" nodes."<<endl;
	//cout<<"  The first molecule can be found at tree index "<< firstMoleculeTreeIndex <<"."<<endl;
}

ReactantTree::~ReactantTree() 
{
	delete [] leftRateFactorSum;
	delete [] leftElementCount;
	delete [] rightElementCount;
	for(unsigned int i=0; i<this->capacity; i++) {
		delete mappingSets[i];
		mappingSets[i]=0;
	}
	delete [] mappingSets;
}
		




MappingSet * ReactantTree::pushNextAvailableMappingSet(double rateFactor)
{
	//Check that we didn't go over the max
	//if(numOfMolecules >= maxElementCount)
	//	cerr<<"Error in NFReactantTree!!!  Adding more than I can take! "<<endl;
		
	unsigned int cn = 1; // index of current node
	
	//Keep going down the tree until we reach the bottom
	while(cn < firstMoleculeTreeIndex)
	{
		//Pick the side of the tree that has the least number
		//of elements, or the left side if they are equal
		if( leftElementCount[cn] <= rightElementCount[cn])
		{
			//Inserting left, so we have to remember the rateFactor...	
			leftElementCount[cn]++;
			leftRateFactorSum[cn] += rateFactor;
			cn = 2*cn;
		}
		else
		{
			//Inserting right, so just remember that...
			rightElementCount[cn]++;
			cn = 2*cn+1;
		}
	}
	
	leftRateFactorSum[cn] = rateFactor;
	leftRateFactorSum[0] += rateFactor;
	unsigned int mappingSetId = cn - firstMoleculeTreeIndex;
	
	
	
	//count the fact that we inserted, then return the result
	n_mappingSets++;
	return mappingSets[mappingSetId];
}
			
			
void ReactantTree::removeMappingSet(unsigned int mappingSetId)
{
	//Go to that position in the tree, and work up and out
	unsigned int cn = mappingSetId + firstMoleculeTreeIndex;

	
	//Get the rate factor from the bottom of the tree and set it to zero
	double rateFactor = leftRateFactorSum[cn];
	leftRateFactorSum[cn] = 0;
	
	//Make sure the tree isn't already empty (otherwise, we've got problems!)
	//and if we happen to empty out the tree, then set the rateFactorSum to exactly
	//zero so we clear away any round off errors
	if(n_mappingSets==0) {
		cerr<<"Trying to remove from an empty ReactantTree!!  Quitting!"<<endl;
		exit(1);
	} else if(n_mappingSets==1) {
		leftRateFactorSum[0] = 0;
	} else {
		leftRateFactorSum[0] -= rateFactor;
	}
	
	
	//Work our way back up to the root
	while(cn>1)
	{
		unsigned int parent = 	cn/2;
		if(cn%2==0)  //Then I was the left child, and we have to make adjustments
		{
			leftElementCount[parent]--;
			leftRateFactorSum[parent] -= rateFactor;
		}
		else  //I was the right child, and we don't have to make adjustments
		{
			rightElementCount[parent]--;
		}				
		
		cn = parent;
	}
		
	//Get rid of the molecule
	n_mappingSets--;
	
	//Notice that we can clear the mapping set if you're paranoid about errors, but
	//really, this doesn't have to happen
	mappingSets[mappingSetId]->clear();
}



void ReactantTree::updateValue(unsigned int mappingSetId, double newRateFactor)
{
	//Here we start from the bottom, removing the old value and adding the new value
	//Go to that position in the tree, and work up and out
	unsigned int cn = mappingSetId + firstMoleculeTreeIndex;
		
	//Get the rate factor from the bottom of the tree and set it to the new rate factor
	double oldRateFactor = leftRateFactorSum[cn];
	leftRateFactorSum[cn] = newRateFactor;
	leftRateFactorSum[0] -= oldRateFactor;
	leftRateFactorSum[0] += newRateFactor;
	
	//Work our way back up to the root
	while(cn>1)
	{
		unsigned int parent = 	cn/2;
		if(cn%2==0)  //Then I was the left child, and we have to make adjustments
		{
			leftRateFactorSum[parent] -= oldRateFactor;
			leftRateFactorSum[parent] += newRateFactor;
		}
		//In this case, the right child doesn't have to do anything
		
		cn = parent;
	}
	
	
}


void ReactantTree::pickReactantFromValue(MappingSet *&ms, double value, double baseRate) const
{
	//First a quick check to make sure we are in bounds (commented out unless we need it)
	//if(value > (leftRateFactorSum[0]*baseRate) )
	//{
	//	cerr<<"Something went wrong::: in NFReactantTree, trying to select a molecule";
	//	cerr<<" with a value greater than the size the total sum"<<endl;
	//	cerr<<" value: " << value;
	//	cerr<<" rateFactorSum: " << leftRateFactorSum[0] << " and total " << (leftRateFactorSum[0]*baseRate) << endl;
	//}
	
	//Start from the top of the tree, and based on the given value, determine what to choose
	unsigned int cn = 1; // index of current node
	
	//Keep going down the tree until we reach the bottom
	while(cn < firstMoleculeTreeIndex)
	{
		//Pick the side of the tree to go down based on the value
		if( value <= (leftRateFactorSum[cn] * baseRate) )
		{
			// Go down the left path
			cn = 2*cn;
		}
		else
		{
			// Go down the right path, but first we have to subtract out all 
			//the left path sums so that we deal with only the remainder
			value -= (leftRateFactorSum[cn] * baseRate);
			cn = 2*cn+1;
		}
	}
	
	//Now we should have the value of cn that gives our molecule, so return it
	unsigned int mappingSetId = cn - firstMoleculeTreeIndex;
	
	//Error check again to ensure things are still ok
	//if(reactants[molArrayIndex]==0)
	//{
	//	cerr<<"!!! in NFReactantTree !!!  could not find a molecule from a given value!"<<endl;
	//	cerr<<"I have: " << numOfMolecules << " molecules! "<<endl;
	//}
		
	ms = mappingSets[mappingSetId];
}
		

		
void ReactantTree::printDetails() const 
{
	cout<<"Reactant Tree Details:"<<endl;
	cout<<"  I have stored " << n_mappingSets << " mappingSets with a total rate sum of ";
	cout<<leftRateFactorSum[0] <<endl;
	for(unsigned int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<i;
	cout<<endl;
	for(unsigned int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<leftRateFactorSum[i];
	cout<<endl;
	for(unsigned int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<leftElementCount[i];
	cout<<endl;
	
	for(unsigned int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<rightElementCount[i];
	cout<<endl;

}
	
