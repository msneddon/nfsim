
#include "NFrateChangeReactions.hh"
#include <math.h>




NFReactantTree::NFReactantTree(unsigned long int maxElementCount)
{
	//cout<<"Creating Reactant Tree..."<<endl;
	//cout<<"  Needs to be able to store "<< maxElementCount <<" molecules, ";
	if(maxElementCount<4) maxElementCount=4;
	
	//Get the depth of the tree, (can cast here because depth will always be a small integer)
	this->treeDepth = (unsigned int)ceil(log(maxElementCount)/log(2)) ;
	
	//Calculate the max number of elements we can store (number of leaves)
	this->maxElementCount = (unsigned int) pow(2,treeDepth);
	
	//Calculate the number of nodes we need (and thus the length of the tree arrays)
	this->numOfNodes = ((unsigned int) pow(2,treeDepth+1)) - 1;
	
	//Calculate the first index that can store a molecule element
	this->firstMoleculeTreeIndex = this->maxElementCount;
	
	//Set up and init our tree arrays
	this->leftRateFactorSum = new double [numOfNodes+1];
	this->leftElementCount = new unsigned int [numOfNodes+1];
	this->rightElementCount = new unsigned int [numOfNodes+1];
	
	for(unsigned int i=0; i<=numOfNodes; i++)
	{
		this->leftRateFactorSum[i] = 0;
		this->leftElementCount[i] = 0;
		this->rightElementCount[i] = 0;
	}
	
	//Set up our molecule array.  The index into this array is the index into the
	//tree minus the firstMoleculeTreeIndex.
	this->reactants = new Molecule * [this->maxElementCount];
	//for(unsigned int i=0; i<maxElementCount; i++) this->reactants[i] = 0;
	this->numOfMolecules = 0;
	
	//cout<<"so setting a limit of " << this->maxElementCount <<" molecules. "<<endl;
	//cout<<"  The depth of the tree will be "<< treeDepth << " and contain ";
	//cout<<numOfNodes<<" nodes."<<endl;
	//cout<<"  The first molecule can be found at tree index "<< firstMoleculeTreeIndex <<"."<<endl;
}

NFReactantTree::~NFReactantTree() 
{
	delete [] leftRateFactorSum;
	delete [] leftElementCount;
	delete [] rightElementCount;
	delete [] reactants;
}
		
int NFReactantTree::insert(Molecule * m, double rateFactor) 
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
	unsigned int molArrayIndex = cn - firstMoleculeTreeIndex;
	
	reactants[molArrayIndex] = m;
	
	//count the fact that we inserted, then return the result
	numOfMolecules++;
	
	return molArrayIndex;
}

void NFReactantTree::remove(Molecule * m, unsigned int rxnListIndex) 
{
	//Go to that position in the tree, and work up and out
	unsigned int cn = rxnListIndex + firstMoleculeTreeIndex;

	//Quick error check
	//if(reactants[rxnListIndex]!=m)
	//	cout<<"we've got problems in remove"<<endl;
		
	//Get the rate factor from the bottom of the tree and set it to zero
	double rateFactor = leftRateFactorSum[cn];
	leftRateFactorSum[cn] = 0;
	
	if(numOfMolecules<=1)
		leftRateFactorSum[0] = 0;
	else
		leftRateFactorSum[0] -= rateFactor;
	
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
	reactants[rxnListIndex] = 0;
	numOfMolecules--;
}

void NFReactantTree::updateValue(Molecule * m, unsigned int rxnListIndex, double newRateFactor)
{
	//Here we start from the bottom, removing the old value and adding the new value
	//Go to that position in the tree, and work up and out
	unsigned int cn = rxnListIndex + firstMoleculeTreeIndex;
	
	//Do an error check here
	//if(reactants[rxnListIndex]!=m) 
	//	cout<<"we've got problems in update"<<endl;
		
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
		
Molecule * NFReactantTree::getReactantFromValue(double value, double baseRate) const 
{
	//First a quick check to make sure we are in bounds (commented out unless we need it)
	//if(value > (leftRateFactorSum[0]*baseRate) )
	//{
	//	cerr<<"Something went wrong::: in NFReactantTree, trying to select a molecule";
	//	cerr<<" with a value greater than the size the total sum"<<endl;
	//	cerr<<" value: " << value;
	//	cerr<<" rateFactorSum: " << leftRateFactorSum[0] << " and total " << (leftRateFactorSum[0]*baseRate) << endl;
	//}
	
	//Start from the top of the tree, and based on the given value, determine
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
	unsigned int molArrayIndex = cn - firstMoleculeTreeIndex;
	
	//Error check again to ensure things are still ok
	//if(reactants[molArrayIndex]==0)
	//{
	//	cerr<<"!!! in NFReactantTree !!!  could not find a molecule from a given value!"<<endl;
	//	cerr<<"I have: " << numOfMolecules << " molecules! "<<endl;
	//}
		
	return reactants[molArrayIndex];
}
		
		
void NFReactantTree::printDetails() const 
{
	cout<<"Reactant Tree Details:"<<endl;
	cout<<"  I have stored " << numOfMolecules << " molecules with a total sum of ";
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
	
