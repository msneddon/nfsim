/*
 * observable.hh
 *
 *  Created on: May 26, 2009
 *      Author: msneddon
 */


#ifndef OBSERVABLE_HH_
#define OBSERVABLE_HH_


#include "NFcore.hh"


namespace NFcore
{

	class ReactionClass;
	class Molecule;
	class TemplateMolecule;
	class Complex;


	//!  Tracks the counts of predefined observables in the simulation.
	/*!
	    Observables keep track of the counts of specific Molecule configurations
	    in the system.  Observables use TemplateMolecules to determine if a Molecule
	    configuration should be counted as an Observable.  Users have the choice of
	    computing Observables on the fly so that they are incrementally updated after
	    each event or recomputing all Observables before each output step.
        @author Michael Sneddon
	 */
	class Observable
	{
		public:
			Observable(string name);
			virtual ~Observable();

			virtual Observable * clone()=0;



			void add();
			void straightAdd();
			void clear() { count=0; };

			/* add multiple new matches, all at once. useful for counter updates --justin */
			void add( int n_matches );
			void straightAdd( int n_matches );
			/* remove multiple new matches, all at once. useful for counter updates --justin */
			void subtract( int n_matches );
			void straightSubtract( int n_matches );

			int getCount() const { return (int)count; };
			string getName() const { return obsName; };
			int getType() const { return type; };

			void getTemplateMoleculeList(int &n_templates, TemplateMolecule **&tmList);

			void addReferenceToMyself(mu::Parser *p);
			void addReferenceToMyself(string referenceName, mu::Parser *p);
			void addDependentRxn(ReactionClass *r);

			virtual int isObservable(Molecule *m) const = 0;
			virtual int isObservable(Complex *c) const = 0;

			//Indentifiers
			static const int NO_RELATION = -1;
			static const int EQUALS = 0;
			static const int NOT_EQUALS = 1;
			static const int GREATER_THAN = 2;
			static const int LESS_THAN = 3;
			static const int GREATOR_OR_EQUAL_TO = 4;
			static const int LESS_THAN_OR_EQUAL_TO = 5;

			static const int NO_TYPE = 0;
			static const int MOLECULES = 1;
			static const int SPECIES = 2;

		protected:
			string obsName;
			int type;
			double count;


			int n_templates;
			TemplateMolecule ** templateMolecules;

			int n_dependentRxns;
			ReactionClass ** dependentRxns;
	};

	class MoleculesObservable : public Observable
	{
		public:
			MoleculesObservable(string name, TemplateMolecule *tm);
			MoleculesObservable(string name, vector <TemplateMolecule *> &tmList);

			virtual ~MoleculesObservable();

			virtual Observable * clone();

			virtual int isObservable(Molecule *m) const;
			virtual int isObservable(Complex *c) const;
	};

	class SpeciesObservable : public Observable
	{
		public:
			SpeciesObservable(string name, vector <TemplateMolecule *> &tmList, vector <string> &stochRelation, vector <int> &stochQuantity);

			virtual ~SpeciesObservable();

			virtual Observable * clone();

			virtual int isObservable(Molecule *m) const;
			virtual int isObservable(Complex *c) const;

		protected:

			// information for processing stochiometric observables
			int *relation;
			int *quantity;
	};
}


#endif /* OBSERVABLE_HH_ */
