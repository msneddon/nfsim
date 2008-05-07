/**
 * AN_hexgrid.cpp
 * This file contains functions to generate a hexgrid of bonds that will
 * allow us to explicitly treat assistance neighborhoods.
 * 
 * 
 * 
 */

#include "receptor_cluster.hh"



/**
 * This function creates a random bonded neighborhood in a hex grid way, meaning it
 * creates a pattern where each molecule is bonded to six others in a regular repeating
 * way, where the base unit looks something like:
 *
 *     M-M-M                                             0   1   2
 *     |\|/|                                              \  |  /           
 *     | M |    and we assume bonds are labeled as:        \ | /
 *     |/|\|                                                 M
 *     M-M-M                                               / | \
 *                                                        /  |  \
 *                                                       5   4   3
 *
 * Then, this method takes the bondSeq vector to specify how to navigate the complex
 * and  attaches all the molecules together in a regular pattern until all molecules
 * have been incorporated.  The standard way to navigate the complex is by starting at
 * the center and working in concentric circles around the grid so that the grid is
 * as packed as possible.
 * 
 * Note that this method is random in that the given vector of molecules is randomized
 * such that for any particular run, this method will create a new arrangement of the
 * molecules (of course the molecules will always be in the same grid pattern as specified
 * by the bondSeq vector.
 *
 * Also note, this method does not alter the bondSeq vector, so it is safe to use that
 * vector multiple times.
 */
void ReceptorCluster::createRandomHexNeighborhood(vector <Molecule *> & molecules, vector <int> bondSeq, int firstBsiteIndex)
{
	int degree = 6;  // in a hex mesh, each molecule is connected to 6 other molecules
	int offset = degree/2; //an offset value to get the index of the bond on the adjacent molecule
	
	//First, shuffle the elements of molecules randomly
	//Here, I use the linear time Fisher-Yates shuffle
	for(int n=(molecules.size()-1); n>0; n--)
	{
		// first pick a random number between [0 and n+1), which lets us pick n
		int r = NFutil::RANDOM_INT(0,n+1);
		// then swap the random position with the last element of the vector
		Molecule * temp = molecules.at(n);
		molecules.at(n) = molecules.at(r);
		molecules.at(r) = temp;
	}
	
	//Finally, build the hex grid
	int bondSeqIndex = 0;
	unsigned int n_molecules = molecules.size();
	unsigned int nextMolecule = 0;
	Molecule * focus = molecules.at(nextMolecule);
	nextMolecule++;
	//int totalBondCount = 0;
	
	while(nextMolecule<n_molecules)
	{
		//cout<<endl<<"----"<<endl;
		//cout<<"focus is now " << focus->getUniqueID()<<endl;
		for(int b=0; b<degree; b++)
		{
			//If the site is open, we have to add another molecule
			if(focus->isBindingSiteOpen(b+firstBsiteIndex))
			{
				//Bind the site to a new random molecule
				int b2 = ((b+offset)%degree)+firstBsiteIndex;
				Molecule::bind(focus,b+firstBsiteIndex,molecules.at(nextMolecule++), b2);
				//cout<<"binding F " << focus->getUniqueID()<<"_"<<b+firstBsiteIndex;
				//cout<< " to " << molecules.at(nextMolecule-1)->getUniqueID() << "_"<<b2<<endl;
				//totalBondCount++;
			}
			
			//Now do the back binding as long as we are past the first b site
			if(b>0)
			{
				Molecule * m1 = focus->getBondedMolecule(b+firstBsiteIndex);
				Molecule * m2 = focus->getBondedMolecule((b-1)+firstBsiteIndex);
				int b1 = (((b+offset)%degree)+1)%degree;
				if(m1->isBindingSiteOpen(b1))
				{
					int b2 = (((b-1)+offset)%degree)-1;
					if(b2<0) b2=degree-1;
					Molecule::bind(m1,b1+firstBsiteIndex,m2,b2+firstBsiteIndex);
					//cout<<"backBinding " << m1->getUniqueID()<<"_"<<b1+firstBsiteIndex;
					//cout<< " to " << m2->getUniqueID() << "_"<<b2+firstBsiteIndex <<endl;
					//totalBondCount++;
				} 
			}
			
			//Finally, complete the wheel with one more bond if we are at the last element
			if(b==(degree-1))
			{
				Molecule * m1 = focus->getBondedMolecule(b+firstBsiteIndex);
				Molecule * m2 = focus->getBondedMolecule(firstBsiteIndex);
				int b1 = ((b+offset)%degree)-1;
				if(b1<0) b1=degree-1;
				if(m1->isBindingSiteOpen(b1))
				{
					int b2 = (((offset)%degree)+1)%degree;
					Molecule::bind(m1,b1+firstBsiteIndex,m2,b2+firstBsiteIndex);
					//cout<<"lastBinding " << m1->getUniqueID()<<"_"<<b1+firstBsiteIndex;
					//cout<< " to " << m2->getUniqueID() << "_"<<b2+firstBsiteIndex <<endl;
					//totalBondCount++;
				}
			}
			if(nextMolecule>=n_molecules) break;
		}
		
		//Now, switch focus to the next focus
		focus = focus->getBondedMolecule(bondSeq.at(bondSeqIndex++));
		//cout<<"focus switching through bond " << bondSeq.at(bondSeqIndex-1) << endl;
		
	}
	//cout<<"bonds created: "<<totalBondCount<<endl;
}



/**
 * This function creates a sequence of bond numbers saved in the given vector that allow
 * you to traverse a hexgrid, starting at the very center of a hexgrid, in concentric circles
 * thereby creating a complete hex neighborhood that is as compact as possible.  These are the
 * order of bonds that you would traverse as you connect the grid.
 * 
 * It assumes the bonds are numbered as follows, but can be 
 * offset by the firstBsiteIndex parameter:
 * 
 *                0   1   2
 *                 \  |  /
 *                  \ | /
 *                    M
 *                  / | \
 *                 /  |  \
 *                5   4   3
 * 
 * Thus, bond 0 will always bond to 3 and so on.
 * This bond sequence will look something like
 * 0,  2,3,4,5,0,1,0,  2,2,3,3,4,4,5,5,0,0,1,1,  0 , ...
 * and will be at least as long as the length parameter given.  
 */
void ReceptorCluster::createStandardHexBondTraversalSeq(vector <int> &bondSeq, int firstBsiteIndex, unsigned int length)
{
	bondSeq.clear();
	int repeatCount = 1;
	int currentBond = 2;
	int degree = 6;  //it must equal 6 for a hex grid
	while(bondSeq.size()<length)
	{
		bondSeq.push_back(firstBsiteIndex);
		//cout<<"putting " << firstBsiteIndex << endl;
		while(currentBond<degree)
		{
			for(int k=0; k<repeatCount; k++)
			{
				bondSeq.push_back(currentBond + firstBsiteIndex);
				//cout<<"putting " << currentBond + firstBsiteIndex << endl;
				if(bondSeq.size()>length) break;
			}
			currentBond++;
			if(bondSeq.size()>length) break;
		}
		currentBond = 0;
		while(currentBond<2)
		{
			for(int k=0; k<repeatCount; k++)
			{
				bondSeq.push_back(currentBond + firstBsiteIndex);
				//cout<<"putting " << currentBond + firstBsiteIndex << endl;
				if(bondSeq.size()>length) break;
			}
			currentBond++;
			if(bondSeq.size()>length) break;
		}
		currentBond=2;
		repeatCount++;
	}
}
