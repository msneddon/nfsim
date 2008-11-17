/*

 */


#include "reactantTree.hh"
#include <math.h>


using namespace NFcore;


ReactantTree::ReactantTree(unsigned int reactantIndex, TransformationSet *ts, unsigned int init_capacity) {
	cout<<"Creating reactant tree... "<<endl;
	cout<<"  Initial Capacity: "<< init_capacity <<" molecules, ";

	//First set basic properties of the tree
	this->reactantIndex=reactantIndex;
	this->ts=ts;

	//set the initial size of the tree
	if(init_capacity<4) maxElementCount=4;
	else maxElementCount = init_capacity;

	//Get the depth of the tree, (can cast here because depth will always be a small integer)
	this->treeDepth = (unsigned int)ceil((double)log((double)maxElementCount)/(double)log((double)2)) ;

	//Calculate the max number of elements we can store (number of leaves)
	this->maxElementCount = (unsigned int) ((double)pow((double)2,(double)treeDepth));

	//Calculate the number of nodes we need (and thus the length of the tree arrays)
	this->numOfNodes = ((unsigned int) ((double)pow((double)2,(double)(treeDepth+1))) ) - 1;

	//Calculate the first index that can store a molecule element - this is useful
	//to have later
	this->firstMappingTreeIndex = this->maxElementCount;

	//Set up and initiate our tree arrays
	this->leftRateFactorSum = new double [numOfNodes+1];
	this->leftElementCount = new int [numOfNodes+1];
	this->rightElementCount = new int [numOfNodes+1];

	for(int i=0; i<=numOfNodes; i++)
	{
		this->leftRateFactorSum[i] = 0;
		this->leftElementCount[i] = 0;
		this->rightElementCount[i] = 0;
	}

	//Set up our mappingSet array.  This array acts as a simple list to manage the
	//mappingSets that exist in the tree.
	this->mappingSets= new MappingSet * [this->maxElementCount];
	for(int i=0; i<maxElementCount; i++)
		mappingSets[i] = ts->generateBlankMappingSet(this->reactantIndex,i);



	msPositionMap = new int [this->maxElementCount];
	for(int i=0; i<maxElementCount; i++)
		msPositionMap[i]=i;  //we start with each element in its rightful position

	msTreePositionMap = new int [this->maxElementCount];
	for(int i=0; i<maxElementCount; i++)
		msTreePositionMap[i]=-1;  //-1 signifies it is not in the tree

	reverseMsTreePositionMap = new int [this->maxElementCount];
	for(int i=0; i<maxElementCount; i++)
		reverseMsTreePositionMap[i]=-1;  //-1 signifies that this position in the tree is empty



	this->n_mappingSets = 0;

	cout<<"so setting a limit of " << this->maxElementCount <<" molecules. "<<endl;
	cout<<"  The depth of the tree will be "<< treeDepth << " and contain ";
	cout<<numOfNodes<<" nodes."<<endl;
	cout<<"  The first molecule can be found at tree index "<< firstMappingTreeIndex <<"."<<endl;

}

void ReactantTree::expandTree(int newCapacity) {
	//Step 1: reallocate new arrays to store the tree, which is the exact same procedure
	//as creating the tree to begin with.  I name everything with the xx_ prefix to make
	//sure I'm working with new variables, not existing member variables.
	int xx_maxElementCount = newCapacity;
	int xx_treeDepth = (unsigned int)ceil((double)log((double)xx_maxElementCount)/(double)log((double)2));
	xx_maxElementCount = (unsigned int) ((double)pow((double)2,(double)xx_treeDepth));
	int xx_numOfNodes = ((unsigned int) ((double)pow((double)2,(double)(xx_treeDepth+1))) ) - 1;

	int xx_firstMappingTreeIndex = xx_maxElementCount;

	double *xx_leftRateFactorSum = new double [xx_numOfNodes+1];
	int *xx_leftElementCount = new int [xx_numOfNodes+1];
	int *xx_rightElementCount = new int [xx_numOfNodes+1];

	for(int i=0; i<=xx_numOfNodes; i++) {
		xx_leftRateFactorSum[i] = 0;
		xx_leftElementCount[i] = 0;
		xx_rightElementCount[i] = 0;
	}

	MappingSet **xx_mappingSets= new MappingSet * [xx_maxElementCount];
	////Take special precaution here!!  we don't want to actually reallocate the mappingSets!
	/// because then we would have to recompare each molecule to this template again!
	/// Instead, we will initialze only the end of this array, and fill in the rest
	/// of the array with the original elements, putting them in the proper position
	/// based on their id
	for(int i=0; i<this->maxElementCount; i++){
		xx_mappingSets[this->mappingSets[i]->getId()] = this->mappingSets[i];
	}
	for(int i=this->maxElementCount; i<xx_maxElementCount; i++) {
		xx_mappingSets[i] = ts->generateBlankMappingSet(this->reactantIndex,i);
	}

	/* original allocation procedure, for reference:
	 * for(int i=0; i<xx_maxElementCount; i++)
	 * 		xx_mappingSets[i] = ts->generateBlankMappingSet(this->reactantIndex,i);*/

	int *xx_msPositionMap = new int [xx_maxElementCount];
	for(int i=0; i<xx_maxElementCount; i++) xx_msPositionMap[i]=i;
	int *xx_msTreePositionMap = new int [xx_maxElementCount];
	for(int i=0; i<xx_maxElementCount; i++) xx_msTreePositionMap[i]=-1;
	int *xx_reverseMsTreePositionMap = new int [xx_maxElementCount];
	for(int i=0; i<xx_maxElementCount; i++) xx_reverseMsTreePositionMap[i]=-1;
	int xx_n_mappingSets = 0;

	//Step 2: Take each mapping set from the original tree, and reinsert them into the new
	//tree.  We cannot just copy over the elements for two reasons: 1) most importantly, the position
	//in the tree will be different because





	//Step 3: Delete all the arrays that we are no longer using to free up the memory
	delete [] this->leftRateFactorSum;
	delete [] this->leftElementCount;
	delete [] this->rightElementCount;
	delete [] this->mappingSets; //remember, just delete the array!  not the actual mappingSets here!
	delete [] this->msPositionMap;
	delete [] this->msTreePositionMap;
	delete [] this->reverseMsTreePositionMap;


	//Step 4: copy the newly created arrays over the original arrays
	this->maxElementCount = xx_maxElementCount;
	this->treeDepth = xx_treeDepth;
	this->numOfNodes = xx_numOfNodes;

	this->leftRateFactorSum = xx_leftRateFactorSum;
	this->leftElementCount = xx_leftElementCount;
	this->rightElementCount = xx_rightElementCount;

	this->mappingSets = xx_mappingSets;

	this->msPositionMap = xx_msPositionMap;

	this->msTreePositionMap = xx_msTreePositionMap;

	this->reverseMsTreePositionMap = xx_reverseMsTreePositionMap;

	this->n_mappingSets=xx_n_mappingSets;

	this->firstMappingTreeIndex = xx_firstMappingTreeIndex;
}


ReactantTree::~ReactantTree()
{
	for(int i=0; i<maxElementCount; i++)
	{
		delete mappingSets[i];
		msTreePositionMap[i]=0;
	}

	delete [] leftRateFactorSum;
	delete [] leftElementCount;
	delete [] rightElementCount;
}






MappingSet * ReactantTree::pushNextAvailableMappingSet() {
	//Check that we didn't go over the max
	if(n_mappingSets >= maxElementCount) {
		cerr<<"Error in ReactantTree!!!  Adding more than I can take! "<<endl;

		//we have to make ourselves bigger...


	}

	n_mappingSets++;
	return mappingSets[n_mappingSets-1];
}


void ReactantTree::confirmPush(int mappingSetId, double rateFactor) {

	unsigned int cn = 1; // index of current node

	//This is where we actually add the pushed mappingSet onto the tree.
	//Keep going down the tree until we reach the bottom, we know we
	//are at the bottom because the current node index will be greater
	// than the firstMoleculeTreeIndex
	while(cn < firstMappingTreeIndex)
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

	//Find the position that we actually want to insert the mapping
	unsigned int msTreeArrayPosition= cn - firstMappingTreeIndex;

	//update our arrays to remember this position in the tree
	//msPositionMap[mappingSetId];  //this does not change
	msTreePositionMap[mappingSetId]=msTreeArrayPosition;  //remember what position in the tree this mappingSet is at
	reverseMsTreePositionMap[msTreeArrayPosition]=mappingSetId;  //remember what mappingSet is at this tree position
}


void ReactantTree::popLastMappingSet() {
	if(n_mappingSets<=0) {
		cerr<<"Trying to pop an empty ReactantTree!!"<<endl;
		exit(1);
	}

	//We check here if the mappingSet that we tried to push was in fact confirmed (by seeing
	//if it had a place in the tree).  If it did have a place in the tree, the MappingSet can
	//only be removed by calling remove, not by popping.
	//cout<<msTreePositionMap[mappingSets[n_mappingSets-1]->getId()]<<endl;
	if(msTreePositionMap[mappingSets[n_mappingSets-1]->getId()]>=0) {
		this->printDetails();
		cout<<"Can't pop the last mappingSet if it was already confirmed to be in the tree!"<<endl;

		exit(1);
	}

	//Clear out the mappingSet (just in case) and decrease the count
	mappingSets[n_mappingSets-1]->clear();
	n_mappingSets--;
}


void ReactantTree::removeMappingSet(unsigned int mappingSetId) {

	//first get the position of this mappingSet in the tree
	int msTreeArrayPosition = msTreePositionMap[mappingSetId];

	//First some error checking...
	if(msTreeArrayPosition<0) {
		cerr<<"Trying to remove a MappingSet from a reactantTree when the MappingSet is not in the tree!"<<endl;
		exit(1);
	}
	if(n_mappingSets==0) {
		cerr<<"Trying to remove from an empty ReactantTree!!"<<endl;
		exit(1);
	}

	//Go to that position in the tree, and work up and out
	unsigned int cn = msTreeArrayPosition + firstMappingTreeIndex;
	cout<<"Removing tree index: "<<cn<<endl;

	//Get the rate factor from the bottom of the tree and set it to zero
	double rateFactor = leftRateFactorSum[cn];
	leftRateFactorSum[cn] = 0;

	if(n_mappingSets<=1)
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

	//Now, remove this guy from the tree array by telling the arrays
	//that this leaf in the tree is empty and this mappingSet is not
	//in the tree
	msTreePositionMap[mappingSetId] = -1;
	reverseMsTreePositionMap[msTreeArrayPosition] = -1;


	//At this point, the tree is up to date, but we still have to get rid of the empty mappingSet
	//by swapping it with the end of the list...  This will allow us to reuse the mappingSet without
	//destroying it and creating it again later...


	//So first, get the position of the mappingSet we need to remove
	int pos = msPositionMap[mappingSetId];

	//Make sure the position is valid (not out of bounds of the List)
	if(pos+1>(n_mappingSets)) {
		cout<<"Error in ReactantTree:  you can't remove a mappingSet that has been cleared! (trying to remove: ";
		cout<< mappingSetId << " in pos " << pos <<" but size is: "<<size()<<endl;
		exit(1);
	}

	//If the array has only one element, or we just happened to select the last element,
	//then just remove the last element without a swap - we can do this by popping the
	//last mapping set.  This will work, because we already cleared the tree at this location
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
	//that MappingSet now lives.
	msPositionMap[mappingSetId] = n_mappingSets-1;
	msPositionMap[mappingSets[pos]->getId()] = pos;

	//Make sure we clear what we don't need
	tempMappingSet->clear();

	//Remember to mark the removal on our counter...
	n_mappingSets--;
}


void ReactantTree::getReactantFromValue(MappingSet *&ms, double value, double baseRate) {

	//First a quick check to make sure we are in bounds (commented out unless we need to debug)
	//if(value > (leftRateFactorSum[0]*baseRate) )
	//{
	//	cerr<<"Something went wrong::: in NFReactantTree, trying to select a molecule";
	//	cerr<<" with a value greater than the size the total sum"<<endl;
	//	cerr<<" value: " << value;
	//	cerr<<" rateFactorSum: " << leftRateFactorSum[0] << " and total " << (leftRateFactorSum[0]*baseRate) << endl;
	//}

	//Start from the top of the tree, and based on the given value, determine
	//where we should end up...
	unsigned int cn = 1; // index of current node

	//Keep going down the tree until we reach the bottom
	while(cn < firstMappingTreeIndex)
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
	unsigned int msTreeArrayPosition = cn - firstMappingTreeIndex;

	//Given the position in the tree, retrieve the mappingSetId
	unsigned int mappingSetId = reverseMsTreePositionMap[msTreeArrayPosition];

	//Given the mappingSetId, we can get the position of the mappingSet in the list
	//of mappingSets
	unsigned int pos = msPositionMap[mappingSetId];

	//Finally, given this postion, we can retrieve the mappingSet
	ms = mappingSets[pos];


	//Error check again to ensure things are still ok...
	//if(reactants[molArrayIndex]==0)
	//{
	//	cerr<<"!!! in ReactantTree !!!  could not find a molecule from a given value!"<<endl;
	//	cerr<<"I have: " << numOfMolecules << " molecules! "<<endl;
	//}
}



void ReactantTree::updateValue(unsigned int mappingSetId, double newRateFactor)
{
	//Here we start from the bottom, removing the old value and adding the new value
	//Go to that position in the tree, and work up and out

	//So first, get the index of this map in the tree which we will set as the
	//current node.  Then we will work back up.
	unsigned int treeIndex = msTreePositionMap[mappingSetId];
	if(treeIndex<0 || treeIndex>maxElementCount) {
		cout<<"Error in ReacantTree! Trying to update a node that is not in the tree!"<<endl;
		exit(1);
	}
	unsigned int cn = treeIndex + this->maxElementCount;

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

	//Ok, we are up to date.  Nothing else changes here...
}


void ReactantTree::printDetails() {

	cout<<endl<<endl<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<"Printing ReactantTree: size="<<size()<<endl<<endl;
	cout<<"MS Array: [ ";
	for(int i=0; i<maxElementCount; i++)
		cout<<mappingSets[i]->getId()<<" ";
	cout<<"]"<<endl;

	cout<<endl<<"To map mappingSetId to position in the MS Array:"<<endl;
	cout<<"MS Pos Map: [ ";
	for(int i=0; i<maxElementCount; i++)
		cout<<msPositionMap[i]<<" ";
	cout<<"]"<<endl;

	cout<<endl<<"To map mappingSetId to position in the tree:"<<endl;
	cout<<"MS TreePos Map: [ ";
	for(int i=0; i<maxElementCount; i++)
		cout<<msTreePositionMap[i]<<" ";
	cout<<"]"<<endl;


	cout<<endl<<"To map position in the tree to mappingSetId:"<<endl;
	cout<<"Reverse MS TreePos Map: [ ";
	for(int i=0; i<maxElementCount; i++)
		cout<<reverseMsTreePositionMap[i]<<" ";
	cout<<"]"<<endl;

	cout<<leftRateFactorSum[0] <<endl;
	for(int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<i;
	cout<<endl;
	for(int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<leftRateFactorSum[i];
	cout<<endl;
	for(int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<leftElementCount[i];
	cout<<endl;

	for(int i=0; i<=numOfNodes; i++)
		cout<<"\t"<<rightElementCount[i];
	cout<<endl;


	//int * msPositionMap;

	//Given a mappingSet Id, this tells us what position it is in
	//in the actual tree
	//int * msTreePositionMap;

	//Given a tree position index, this tells us the mappingSet Id
	//which we can use to get our mappingSet out of the mappingSets
	//array
	//int * reverseMsTreePositionMap;


}


//int NFReactantTree::insert(Molecule * m, double rateFactor)
//{
//	//Check that we didn't go over the max
//	if(numOfMolecules >= maxElementCount)
//		cerr<<"Error in ReactantTree!!!  Adding more than I can take! "<<endl;
//
//	unsigned int cn = 1; // index of current node
//
//	//Keep going down the tree until we reach the bottom
//	while(cn < firstMoleculeTreeIndex)
//	{
//		//Pick the side of the tree that has the least number
//		//of elements, or the left side if they are equal
//		if( leftElementCount[cn] <= rightElementCount[cn])
//		{
//			//Inserting left, so we have to remember the rateFactor...
//			leftElementCount[cn]++;
//			leftRateFactorSum[cn] += rateFactor;
//			cn = 2*cn;
//		}
//		else
//		{
//			//Inserting right, so just remember that...
//			rightElementCount[cn]++;
//			cn = 2*cn+1;
//		}
//	}
//
//	leftRateFactorSum[cn] = rateFactor;
//	leftRateFactorSum[0] += rateFactor;
//	unsigned int molArrayIndex = cn - firstMoleculeTreeIndex;
//
//	reactants[molArrayIndex] = m;
//
//	//count the fact that we inserted, then return the result
//	numOfMolecules++;
//
//	return molArrayIndex;
//}
//
//void NFReactantTree::remove(Molecule * m, unsigned int rxnListIndex)
//{
//	//Go to that position in the tree, and work up and out
//	unsigned int cn = rxnListIndex + firstMoleculeTreeIndex;
//
//	//Quick error check
//	//if(reactants[rxnListIndex]!=m)
//	//	cout<<"we've got problems in remove"<<endl;
//
//	//Get the rate factor from the bottom of the tree and set it to zero
//	double rateFactor = leftRateFactorSum[cn];
//	leftRateFactorSum[cn] = 0;
//
//	if(numOfMolecules<=1)
//		leftRateFactorSum[0] = 0;
//	else
//		leftRateFactorSum[0] -= rateFactor;
//
//	//Work our way back up to the root
//	while(cn>1)
//	{
//		unsigned int parent = 	cn/2;
//		if(cn%2==0)  //Then I was the left child, and we have to make adjustments
//		{
//			leftElementCount[parent]--;
//			leftRateFactorSum[parent] -= rateFactor;
//		}
//		else  //I was the right child, and we don't have to make adjustments
//		{
//			rightElementCount[parent]--;
//		}
//
//		cn = parent;
//	}
//
//	//Get rid of the molecule
//	reactants[rxnListIndex] = 0;
//	numOfMolecules--;
//}
//
//void NFReactantTree::updateValue(Molecule * m, unsigned int rxnListIndex, double newRateFactor)
//{
//	//Here we start from the bottom, removing the old value and adding the new value
//	//Go to that position in the tree, and work up and out
//	unsigned int cn = rxnListIndex + firstMoleculeTreeIndex;
//
//	//Do an error check here
//	//if(reactants[rxnListIndex]!=m)
//	//	cout<<"we've got problems in update"<<endl;
//
//	//Get the rate factor from the bottom of the tree and set it to the new rate factor
//	double oldRateFactor = leftRateFactorSum[cn];
//	leftRateFactorSum[cn] = newRateFactor;
//	leftRateFactorSum[0] -= oldRateFactor;
//	leftRateFactorSum[0] += newRateFactor;
//
//	//Work our way back up to the root
//	while(cn>1)
//	{
//		unsigned int parent = 	cn/2;
//		if(cn%2==0)  //Then I was the left child, and we have to make adjustments
//		{
//			leftRateFactorSum[parent] -= oldRateFactor;
//			leftRateFactorSum[parent] += newRateFactor;
//		}
//		//In this case, the right child doesn't have to do anything
//
//		cn = parent;
//	}
//}
//
//Molecule * NFReactantTree::getReactantFromValue(double value, double baseRate) const
//{
//	//First a quick check to make sure we are in bounds (commented out unless we need it)
//	//if(value > (leftRateFactorSum[0]*baseRate) )
//	//{
//	//	cerr<<"Something went wrong::: in NFReactantTree, trying to select a molecule";
//	//	cerr<<" with a value greater than the size the total sum"<<endl;
//	//	cerr<<" value: " << value;
//	//	cerr<<" rateFactorSum: " << leftRateFactorSum[0] << " and total " << (leftRateFactorSum[0]*baseRate) << endl;
//	//}
//
//	//Start from the top of the tree, and based on the given value, determine
//	unsigned int cn = 1; // index of current node
//
//	//Keep going down the tree until we reach the bottom
//	while(cn < firstMoleculeTreeIndex)
//	{
//		//Pick the side of the tree to go down based on the value
//		if( value <= (leftRateFactorSum[cn] * baseRate) )
//		{
//			// Go down the left path
//			cn = 2*cn;
//		}
//		else
//		{
//			// Go down the right path, but first we have to subtract out all
//			//the left path sums so that we deal with only the remainder
//			value -= (leftRateFactorSum[cn] * baseRate);
//			cn = 2*cn+1;
//		}
//	}
//
//	//Now we should have the value of cn that gives our molecule, so return it
//	unsigned int molArrayIndex = cn - firstMoleculeTreeIndex;
//
//	//Error check again to ensure things are still ok
//	//if(reactants[molArrayIndex]==0)
//	//{
//	//	cerr<<"!!! in NFReactantTree !!!  could not find a molecule from a given value!"<<endl;
//	//	cerr<<"I have: " << numOfMolecules << " molecules! "<<endl;
//	//}
//
//	return reactants[molArrayIndex];
//}
//
//
//void NFReactantTree::printDetails() const
//{
//	cout<<"Reactant Tree Details:"<<endl;
//	cout<<"  I have stored " << numOfMolecules << " molecules with a total sum of ";
//	cout<<leftRateFactorSum[0] <<endl;
//	for(unsigned int i=0; i<=numOfNodes; i++)
//		cout<<"\t"<<i;
//	cout<<endl;
//	for(unsigned int i=0; i<=numOfNodes; i++)
//		cout<<"\t"<<leftRateFactorSum[i];
//	cout<<endl;
//	for(unsigned int i=0; i<=numOfNodes; i++)
//		cout<<"\t"<<leftElementCount[i];
//	cout<<endl;
//
//	for(unsigned int i=0; i<=numOfNodes; i++)
//		cout<<"\t"<<rightElementCount[i];
//	cout<<endl;
//
//}




