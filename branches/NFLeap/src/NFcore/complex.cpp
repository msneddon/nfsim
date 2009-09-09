#include <iostream>
#include "NFcore.hh"
#include <string>

using namespace std;
using namespace NFcore;


Complex::Complex(System * s, int ID_complex, Molecule * m)
{
	this->system = s;
	this->ID_complex = ID_complex;
	this->complexMembers.push_back(m);
}

Complex::~Complex()
{
}

bool Complex::isAlive() {
	if(complexMembers.size()==0) return false;
	return (*complexMembers.begin())->isAlive();
}


int Complex::getMoleculeCountOfType(MoleculeType *m)
{
	int count = 0;
	for( molIter = complexMembers.begin(); molIter != complexMembers.end(); molIter++ )
	{
		if((*molIter)->getMoleculeTypeName()==m->getName())
		{
			//cout<<count<<" Match! : "<<m->getName()<<" with "<< (*molIter)->getMoleculeTypeName()<<endl;
			count++;
		}
	}
  	return count;
}


void Complex::printDetails()
{
	cout<<"   -Complex "<<ID_complex<<": ("<<complexMembers.size()<<") -";
	for( molIter = complexMembers.begin(); molIter != complexMembers.end(); molIter++ )
	{
  		cout<<" "<<(*molIter)->getMoleculeTypeName()<<"_";
		cout<<"_u"<<(*molIter)->getUniqueID();
	}
  	cout<<endl;
}


void Complex::printDetailsLong()
{
	cout<<"   -Complex "<<ID_complex<<": ("<<complexMembers.size()<<") --------------------------\n";
	for( molIter = complexMembers.begin(); molIter != complexMembers.end(); molIter++ )
	{
  		(*molIter)->printDetails();
  		cout<<"degree check: " << (*molIter)->getDegree()<<endl;
	}
  	cout<<endl;
}

void Complex::getDegreeDistribution(vector <int> &degreeDist)
{
	for( molIter = complexMembers.begin(); molIter != complexMembers.end(); molIter++ )
	{
  		int d = (*molIter)->getDegree();
  		while(d>=(int)degreeDist.size())
  			degreeDist.push_back(0);
  		degreeDist.at(d)++;
	}
}

void Complex::printDegreeDistribution()
{
	vector <int> degreeDist;
	vector <int>::iterator degIter;
	getDegreeDistribution(degreeDist);
	cout<<"Degree Distribution for complex "<< ID_complex<<", size: "<<complexMembers.size()<<endl;
	cout<<"  Degree:";
	for(int d=0; d<(int)degreeDist.size(); d++)
		cout<<"\t"<<d;
	cout<<endl<<"  Count:";
	for( degIter = degreeDist.begin(); degIter != degreeDist.end(); degIter++ )
  		cout<<"\t"<<(*degIter);
	cout<<endl;
}


void Complex::refactorToNewComplex(int new_ID_complex)
{
	for( molIter = complexMembers.begin(); molIter != complexMembers.end(); molIter++ )
  		(*molIter)->moveToNewComplex(new_ID_complex);
}

/* for binding, we want to merge a new complex, c, with our complex, this */
void Complex::mergeWithList(Complex * c)
{
	c->refactorToNewComplex(this->ID_complex);
	this->complexMembers.splice(complexMembers.end(),c->complexMembers);
	system->notifyThatComplexIsAvailable(c->getComplexID());
}



/* class to decide when a molecule is in the wrong complex (will tell us to delete this molecule
 * from this complex */
class IsInWrongComplex
{
	public:
		IsInWrongComplex(int currentComplexID) : ID(currentComplexID) {};
		bool operator() (Molecule * m) const { return m->getComplexID()!=ID; };

	private:
		int ID;
};

//int counter=0;
//int totalSizeSum=0;
//int avgTraversalSize = 0;
/* for unbinding, we have to figure out the elements of the new complex,
 * put those elements in the new complex, renumber the complex_id for those
 * molecules, and delete those molecules from this complex.  wheh! */
void Complex::updateComplexMembership(Molecule * m)
{
	//Check if this molecule is indeed in this complex first, can be removed later for
	//optimization
	if(m->getComplexID()!=this->ID_complex) { cerr<< "ERROR IN COMPLEX!!! "<<endl; return; }



	//Get list of things this molecule is still connected to
	list <Molecule *> members;
	m->traverseBondedNeighborhood(members, ReactionClass::NO_LIMIT);

	//counter++;
	//cout<<"traversing neighborhood: "<<counter<<endl;
	//totalSizeSum+=members.size();
	//avgTraversalSize = totalSizeSum/counter;
	//cout<<members.size()<<endl;
	//cout<<"average: "<<avgTraversalSize<<endl;


	//Check if we even need to create a new complex (if not, return)
	if(members.size()==(unsigned)this->getComplexSize())
	{
		//cout<<"still in same complex, so no new complex"<<endl;
		return;
	}

	//Get the next available complex
	Complex *newComplex = system->getNextAvailableComplex();
	//cout<<" forming new complex:  next available: " <<newComplex->getComplexID()<<endl;

	//renumber our complex elements
	list <Molecule *>::iterator molIter;
	for( molIter = members.begin(); molIter != members.end(); molIter++ ) {
		(*molIter)->moveToNewComplex(newComplex->getComplexID());
	}

	//put our new complex elements into that complex
	newComplex->complexMembers.splice(newComplex->complexMembers.end(),members);
	//cout<<"size of list now: " << members.size() <<endl;

	//remove all molecules from this that don't have the correct complex id
	complexMembers.remove_if(IsInWrongComplex(this->ID_complex));



	//update new complex in reactions?

	//

	//done!
}





