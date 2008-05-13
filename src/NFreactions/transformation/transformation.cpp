#include "transformation.hh"


using namespace NFcore;

unsigned int Transformation::currentIndex=0;









Transformation::Transformation(ReactionClass *r, unsigned int type, unsigned int n_mappings)
{
	this->type = type;
	this->n_mappings = n_mappings;
	this->t_index=Transformation::currentIndex++;
}


Transformation::~Transformation()
{
	delete [] templateMappings;
	
}

void Transformation::printDetails()
{
	cout<<t_index<<":: Transformation of type: "<<getTypeName()<<endl;
	cout<<"\tn_mappings: "<<n_mappings<<endl;
	if(templateMappings==NULL)
		cout<<"\tNo templateMappings are linked to this Transformation yet."<<endl;
	else
		for(unsigned int m=0; m<n_mappings; m++)
		{
			if(templateMappings[m]!=NULL)
			{
				cout<<"\tMapping "<<m+1<<" :: Linked to moleculeType ";
				cout<<templateMappings[m]->getMoleculeType()->getName();
				if(templateMappings[m]->getMappingType()==Mapping::STATE)
				{
					cout<<" - ( state = ";
					cout<<templateMappings[m]->getMoleculeType()->getStateName(templateMappings[m]->getMappingIndex());
					cout<<" )"<<endl;
				}
				else
				{
					cout<<" - ( bond = ";
					cout<<templateMappings[m]->getMoleculeType()->getBindingSiteName(templateMappings[m]->getMappingIndex());
					cout<<" )"<<endl;
				}
			}
		}
}

const char * Transformation::getTypeName() const
{
	if(type==Transformation::STATE_CHANGE)
		return "State Change";
	else if(type==Transformation::BINDING)
		return "Binding";
	else if(type==Transformation::UNBINDING)
		return "Unbinding";
	else if(type==Transformation::DECREMENT_STATE)
		return "Decrement State";
	else if(type==Transformation::INCREMENT_STATE)
		return "Increment State";
	else
		return "???";
	
}



Transformation * Transformation::genStateChangeTransform(TemplateMolecule *tm, const char *stateName, int stateValue, ReactionClass *r)
{
		TemplateMapping * tMapping = new TemplateMapping(tm,Mapping::STATE,stateName);
		Transformation *t = new StateChangeTransform(tMapping, stateValue,r);
		r->registerTransformation(t);
		return t;
}

Transformation * Transformation::genBindingTransform(TemplateMolecule *tm1, TemplateMolecule *tm2, const char *bSiteName1, const char *bSiteName2, ReactionClass *r)
{
	
		TemplateMapping * tMapping1 = new TemplateMapping(tm1,Mapping::BOND,bSiteName1);
		TemplateMapping * tMapping2 = new TemplateMapping(tm2,Mapping::BOND,bSiteName2);
		Transformation *t = new BindingTransform(tMapping1,tMapping2,r);
		r->registerTransformation(t);
		return t;
	
}



Transformation * Transformation::genUnbindingTransform(TemplateMolecule *tm, const char *bSiteName, ReactionClass *r)
{
	//	ReactionClass *r = tm->getReactionClass();
		TemplateMapping * tMapping = new TemplateMapping(tm,Mapping::BOND,bSiteName);
		Transformation *t = new UnbindingTransform(tMapping,r);
		r->registerTransformation(t);
		return t;
}















StateChangeTransform::StateChangeTransform(TemplateMapping *tm, int newStateValue, ReactionClass *r) :
	Transformation(r,Transformation::STATE_CHANGE,1)
{
	this->newStateValue = newStateValue;
	this->templateMappings = new TemplateMapping * [1] ;
	this->templateMappings[0] = tm;
	tm->setTransformation(this);
	
}

void StateChangeTransform::transform(Mapping **m)
{
	if(m[0]==NULL)
	{
		cout<<"StateChangeTransform was given a null mapping to transform!!"<<endl;
		this->printDetails();
		exit(1);
	}
	if(m[0]->getMolecule()==NULL)
	{
		cout<<"StateChangeTransform was given a mapping without a molecule!!"<<endl;
		this->printDetails();
		exit(1);
	}
	m[0]->getMolecule()->setState(m[0]->getIndex(),newStateValue);
}



BindingTransform::BindingTransform(TemplateMapping *tm1, TemplateMapping *tm2, ReactionClass *r) :
	Transformation(r,Transformation::BINDING,2)
{
	this->templateMappings = new TemplateMapping * [2] ;
	this->templateMappings[0] = tm1;
	this->templateMappings[1] = tm2;
	tm1->setTransformation(this);
	tm2->setTransformation(this);
}

void BindingTransform::transform(Mapping **m)
{
	if(m[0]==NULL || m[1]==NULL)
	{
		cout<<"StateChangeTransform was given a null mapping to transform!!"<<endl;
		this->printDetails();
		exit(1);
	}
	if(m[0]->getMolecule()==NULL || m[1]->getMolecule()==NULL)
	{
		cout<<"StateChangeTransform was given a mapping without a molecule!!"<<endl;
		this->printDetails();
		exit(1);
	}
	
	int bSite0 = m[0]->getIndex();
	int bSite1 = m[1]->getIndex();
	Molecule::bind(m[0]->getMolecule(), bSite0, m[1]->getMolecule(), bSite1);
}



UnbindingTransform::UnbindingTransform(TemplateMapping *tm, ReactionClass *r) :
	Transformation(r,Transformation::UNBINDING,1)
{
	this->templateMappings = new TemplateMapping * [1] ;
	this->templateMappings[0] = tm;
	tm->setTransformation(this);
}

void UnbindingTransform::transform(Mapping **m)
{
	if(m[0]==NULL)
	{
		cout<<"UnbindingTransform was given a null mapping to transform!!"<<endl;
		this->printDetails();
		exit(1);
	}
	if(m[0]->getMolecule()==NULL)
	{
		cout<<"UnbindingTransform was given a mapping without a molecule!!"<<endl;
		this->printDetails();
		exit(1);
	}
	
	int bSite = m[0]->getIndex();
	Molecule::unbind(m[0]->getMolecule(), bSite);
}
