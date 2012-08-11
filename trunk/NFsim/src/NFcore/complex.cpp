#include <iostream>
#include "NFcore.hh"
#include <string>

// for labeling
#include <sstream>
#include <algorithm>
#include "../nauty24/nausparse.h"

using namespace std;
using namespace NFcore;

const int Node::IS_MOLECULE = -1;

Complex::Complex(System * s, int ID_complex, Molecule * m)
	: is_canonical( false ), canonical_label("")
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
	cout << "label: " << getCanonicalLabel() << endl;
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
	// turn off canonical flag
	this->unsetCanonical();
	c->unsetCanonical();

	// move molecules in c to this complex
	c->refactorToNewComplex(this->ID_complex);
	this->complexMembers.splice(complexMembers.end(),c->complexMembers);
	(system->getAllComplexes()).notifyThatComplexIsAvailable(c->getComplexID());
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

	unsetCanonical();

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
	// NETGEN -- redirected call to ComplexList object at system->allComplexes
	Complex *newComplex = (system->getAllComplexes()).getNextAvailableComplex();
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


// get the canonical label for this complex
string Complex::getCanonicalLabel ( )
{
	if (!is_canonical)
		generateCanonicalLabel();

	return canonical_label;
}


// generate a canonical label
/*  1) construct a Nauty sparse graph representation of the complex.
    2) call Nauty to get canonical node order.
    3) construct label based on canonical order.
    Steps 1 and 3 have complexity around E log(V), where E is the number of edges
    and V is the number of vertices (nodes).
 */
void Complex::generateCanonicalLabel ( )
{
    #if DEBUG_NAUTY==1
    std::cout << "find_canonical_order" << std::endl;
    #endif

    // handle special case: complexes with 0 members
    if ( complexMembers.size() == 0 )
    {   // set canonical label
        canonical_label = string("");
        is_canonical = true;
        return;
    }

    // stringstream object where label is constructed
    stringstream  labelstream;
    // declare containers
    vector < Node * >  nodes;
    map < node_t, Node * >  node_index;
    // declare iterators
    vector < Node * >::iterator       node_iter;
    map < node_t, Node * >::iterator  node_index_iter;
    // other variables
    bool      nauty_required;
    Node      *curr_node, *prev_node;
    Molecule  *mol;
    MoleculeType *moltype;
    int       icomp, icomp_nbh;
    // nauty related variables
    int      nv, m, nde;
    int      v_index, e_index;
    int      *lab, *ptn, *orbits;
    setword  *workspace;

	// construct Node objects
    nde = 0;  /* keep track of edges */
	for ( molIter = complexMembers.begin(); molIter != complexMembers.end(); ++molIter )
	{
		mol = *molIter;
		moltype = mol->getMoleculeType();
		nodes.push_back( new Node(mol, Node::IS_MOLECULE) );
        node_index.insert( node_index_t( node_t(mol, Node::IS_MOLECULE), nodes.back() ) );

    	for ( icomp=0; icomp < moltype->getNumOfComponents(); ++icomp )
		{
			nodes.push_back( new Node(mol, icomp) );
            node_index.insert( node_index_t( node_t(mol, icomp), nodes.back() ) );

            nde += 2;  /* add two edges (mol->comp) and (comp->mol) */
			if ( mol->isBindingSiteBonded(icomp) )  ++nde;
        }
	}

    // get number of vertices
    nv = nodes.size();
    // sort nodes (N logN)
    std::sort ( nodes.begin(), nodes.end(), Node::less_by_label );
    // index nodes (N)
    for ( v_index=0, node_iter = nodes.begin();  node_iter != nodes.end();  ++node_iter, ++v_index )
        (*node_iter)->setIndex( v_index );


    // declare various data elements for Nauty
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    statsblk stats;
    // Select option for canonical labelling
    options.getcanon   = TRUE;
    options.defaultptn = FALSE;
    // SIMPLE GRAPH:
    options.digraph = FALSE;
    // define m (see Nauty documentation)
    m = (nv + WORDSIZE - 1) / WORDSIZE;
    nauty_check( WORDSIZE, m, nv, NAUTYVERSIONID );

    // allocate memory for sparse graph
    SG_DECL(sg);
    SG_DECL(cg);
    SG_ALLOC( sg, nv, nde, "malloc" );
    sg.nv  = nv;   /* Number of vertices */
    sg.nde = nde;  /* Number of directed edges */
    // allocate memory for other Nauty data structures
    orbits     = new int [nv];
    workspace  = new setword [10*m];
    lab        = new int [nv];
    ptn        = new int [nv];

    // build sparse graph representation for Nauty ( N + sum_N[ deg(n) logN ] )
    nauty_required = false;
    for ( e_index=0, v_index=0, node_iter = nodes.begin();  node_iter != nodes.end();  ++node_iter, ++v_index )
    {
        curr_node = *node_iter;
        mol   = curr_node->getMolecule();
        moltype = mol->getMoleculeType();
        icomp = curr_node->getComponent();

        // this vertex is labeled v_index
        lab[v_index]  = v_index;
        // edges for this vertex begin at e_index
        sg.v[v_index] = e_index;

        // compare current node to the previous node
        if ( v_index != 0 )
        {
            if ( Node::less_by_label(prev_node, curr_node) )
            {
                ptn[v_index-1] = 0;
            }
            else
            {
                nauty_required = true;
                ptn[v_index-1] = 1;
            }
        }

        // build edges
        if ( curr_node->isMolecule() )
        {   // this is a molecule
    	    for ( int icomp_nbh=0; icomp_nbh < moltype->getNumOfComponents(); ++icomp_nbh )
		    {
                node_index_iter = node_index.find( node_t(mol, icomp_nbh) );
                sg.e[ e_index ] = node_index_iter->second->getIndex();
                ++e_index;
            }
        }
        else
        {   // this is a component
            node_index_iter = node_index.find( node_t( mol, Node::IS_MOLECULE ) );
            sg.e[ e_index ] = node_index_iter->second->getIndex();
            ++e_index;

			if ( mol->isBindingSiteBonded( icomp ) )
			{
                node_index_iter = node_index.find( node_t( mol->getBondedMolecule(icomp),
                                                       mol->getBondedMoleculeBindingSiteIndex(icomp) ) );
                sg.e[ e_index ] = node_index_iter->second->getIndex();
                ++e_index;
			}
        }

        // calculate degree of this node
        sg.d[v_index] = ( v_index == 0 ? e_index : e_index-(sg.v[v_index-1]+sg.d[v_index-1]) );
        // remember this node for comparison
        prev_node = curr_node;
    }

    // assignment the last value of the partition to 0.
    ptn[nv-1] = 0;


    #if DEBUG_NAUTY==1
    cout << "lab: ";
    for ( int i = 0; i != nv; ++i )  cout << lab[i] << ",";
    cout << "\nptn: ";
    for ( int i = 0; i != nv; ++i )  cout << ptn[i] << ",";
    cout << "\nd: ";
    for ( int i = 0; i != nv; ++i )  cout << sg.d[i] << ",";
    cout << "\nv: ";
    for ( int i = 0; i != nv; ++i )  cout << sg.v[i] << ",";
    cout << "\ne: ";
    for ( int i = 0; i != nde; ++i )  cout << sg.e[i] << ",";
    cout << "\nnv:  " << sg.nv;
    cout << "\nnde: " << sg.nde << endl;
    #endif


    if ( nauty_required )
    {  /*  Label sg, result in cg and labelling in lab.
        *  It is not necessary to pre-allocate space in cg1 and cg2, but
        *  they have to be initialised as we did above.
        */
        nauty( (graph*)&sg, lab, ptn, NULL, orbits, &options, &stats,
                      workspace, 10*m, m, nv, (graph*)&cg );
    }

    #if DEBUG_NAUTY==1
    cout << "lab, post-nauty: ";
    for ( int i = 0; i != nv; ++i )  cout << lab[i] << ",";
    cout << std::endl;
    #endif

    // Reindex nodes (N)
    for ( v_index = 0;  v_index < nv;  ++v_index )
        nodes.at( lab[v_index] )->setIndex( v_index );

    // Free Nauty data
    SG_FREE( sg );
    SG_FREE( cg );
    delete [] orbits;
    delete [] workspace;
    delete [] lab;
    delete [] ptn;

    // sort nodes
    std::sort ( nodes.begin(), nodes.end(), Node::less_by_index );

    // Build label  ( N + sum_N[ Edges(n)*logN ] )
    for ( node_iter = nodes.begin();  node_iter != nodes.end();  ++node_iter )
    {
        curr_node = *node_iter;
        mol = curr_node->getMolecule();
        moltype = mol->getMoleculeType();
        icomp = curr_node->getComponent();

        labelstream << curr_node->getLabel();

        if ( curr_node->isMolecule() )
        {   // this is a molecule
    	    for ( icomp_nbh=0; icomp_nbh < moltype->getNumOfComponents(); ++icomp_nbh )
		    {
                node_index_iter = node_index.find( node_t(mol, icomp_nbh) );
                labelstream << "!" << node_index_iter->second->getIndex();
            }
        }
        else
        {   // this is a component
            node_index_iter = node_index.find( node_t( mol, Node::IS_MOLECULE ) );
            labelstream << "!" << node_index_iter->second->getIndex();

			if ( mol->isBindingSiteBonded( icomp ) )
			{
                node_index_iter = node_index.find( node_t( mol->getBondedMolecule(icomp),
                                                       mol->getBondedMoleculeBindingSiteIndex(icomp) ) );
                labelstream << "!" << node_index_iter->second->getIndex();
			}
        }
        labelstream << ",";
    }

    // Free nodes
    for ( node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter )
        delete *node_iter;

    // set canonical label
    canonical_label = labelstream.str();
    is_canonical = true;
}

