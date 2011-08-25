#ifndef TRANSFORMATION_HH_
#define TRANSFORMATION_HH_


#include "../NFreactions.hh"

namespace NFcore
{

	class MoleculeCreator;
	class SpeciesCreator;
	class Mapping;
	class Transformation;
	class AddMoleculeTransform;
	class AddSpeciesTransform;


	/*!
	 	Static object that generates transformations with the given sites.
	    @author Michael Sneddon
	 */
	class TransformationFactory {
		public:
			/*!
			 	Generates a state change Transformation for transforming the given state at the given
			 	index into the new state value.
			    @author Michael Sneddon
			 */
			static Transformation * genStateChangeTransform(unsigned int cIndex, int newValue);

			/*!
			 	Generates a binding transformation for one of the two binding sites in the binding Transform.  You will
			 	have to tell it where the other Transformation object lives (TransformationSet has this information).
			    @author Michael Sneddon
			 */
			static Transformation * genBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex);

			/*!
			 	Generates a binding transformation for one of the two binding sites in the binding Transform, where at least
			 	one of the molecules is a newly created object, so that null condition checks are omitted. 	You will still have
			 	to tell it where the other Transformation object lives (TransformationSet has this information).
			    @author Michael Sneddon
			*/
			static Transformation * genNewMoleculeBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex);


			// deprecated
			//static Transformation * genBindingSeparateComplexTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex);



			/*!
			 	Generates the second half of a binding transform.  The other site already knows about this site, so all you
			 	need here is the index of the binding site that must be bonded.
			    @author Michael Sneddon
			 */
			static Transformation * genBindingTransform2(unsigned int bSiteIndex);

			/*!
			 	Generates an unbinding transformation for a particular binding site.  Only
			 	one of the binding sites must be specified as a Transformation - the other
			 	is automatically taken care of.
			    @author Michael Sneddon
			 */
			static Transformation * genUnbindingTransform(unsigned int bSiteIndex);

			/*!
			 	Generates an Add Molecule transformation.
			    @author Michael Sneddon

			    Changes by JustinHogg, 3Mar2011:
			    Unlike other generators, this now returns a specific child class.
			    AddMoleculeTransforms are handled seperately from other transforms, so
			    we don't need to follow the polymorphic interface. Also, we need to
			    call the method "create_and_map" which is specifc to AddMoleculeTransform.
			 */
			static AddSpeciesTransform  * genAddSpeciesTransform(  SpeciesCreator  * sc );
			static AddMoleculeTransform * genAddMoleculeTransform( MoleculeCreator * mc );


			/*!
			 	Generates a removal of a molecule from the system.  The removalType specifies
			 	how the molecule should be removed (either everything that is connected, or just
			 	the molecules, or just the molecules conditional on how it is connected)
			    @author Michael Sneddon
			 */
			static Transformation * genRemoveMoleculeTransform(int removalType);
			/*!
			 	Generates an empty transformation.  This is used in cases where there is
			 	a reactant that is not transformed in a reaction, but that still needs
			 	to be counted and marked so that the rate of the reaction is correct.
			    @author Michael Sneddon
			 */
			static Transformation * genEmptyTransform();

			/*!
			 	Generates an IncrementState transformation.
			    @author Michael Sneddon
			 */
			static Transformation * genIncrementStateTransform(unsigned int cIndex);

			/*!
			 	Generates an DecrementState transformation.
			    @author Michael Sneddon
			*/
			static Transformation * genDecrementStateTransform(unsigned int cIndex);

			/*!
			 	Generates an IncrementPopulation transformation.
			    @author Justin Hogg
			*/
			//static Transformation * genIncrementPopulationTransform();


			/*!
			 	Generates an DecrementPopulation transformation.
			    @author Justin Hogg
			*/
			static Transformation * genDecrementPopulationTransform();


			/*! Indicates that a delete transform deletes the entire connected species */
			static const int COMPLETE_SPECIES_REMOVAL = 0;

			/*! Indicates that a delete transform deletes only the pointed-to molecule */
			static const int DELETE_MOLECULES = 1;

			/*! Delete only pointed-to molecules, only if deleting it creates only one remaining species*/
			static const int DELETE_MOLECULES_NO_KEYWORD = 2;






			/*!	Indicates a state change transformation or mapping onto a state	*/
			static const unsigned int STATE_CHANGE = 0;

			/*!	Indicates a binding transformation or mapping onto a binding site	*/
			static const unsigned int BINDING = 1;

			/*!	Indicates an unbinding transformation or mapping onto a binding site	*/
			static const unsigned int UNBINDING = 2;

			/*!	Indicates a removal transform or mapping onto an entire Molecule	*/
			static const unsigned int REMOVE = 3;

			/*!	Indicates an addition transform	*/
			static const unsigned int ADD = 4;

			/*!	Indicates no transformation is needed (or is the second partner
			    in a binding transform and so should be skipped when applying transforms	*/
			static const unsigned int EMPTY = 5;

			/*!	Indicates an increment transformation */
			static const unsigned int INCREMENT_STATE = 6;


			/*!	Indicates a decrement transformation */
			static const unsigned int DECREMENT_STATE = 7;



			////////////////////////////////
			/*!	Return a pointer to the transformation that is needed by a local function */
			static Transformation * genLocalFunctionReference(string PointerName, int type, TemplateMolecule *tm);
			static const unsigned int LOCAL_FUNCTION_REFERENCE = 8;
			////////////////////////////////

			/*!	Indicates a population increment transform */
			static const unsigned int INCREMENT_POPULATION = 9;

			/*!	Indicates a population decrement transform */
			static const unsigned int DECREMENT_POPULATION = 10;

	};



	/*!
	 	Abstract transformation object that other types of transformations inherit from.
	    @author Michael Sneddon
	 */
	class Transformation {
		public:
			Transformation(int type) {this->type=type;};
			virtual ~Transformation() {};
			int getType() const { return type; }
			virtual void apply(Mapping *m, MappingSet **ms) = 0;
			virtual int getComponentIndex() const = 0;
			virtual int getRemovalType() { return -1; };
			// returns false if it does not meet a null condition, true if the reaction
			// should be rejected do to a null condition
			virtual bool checkForNullCondition(Mapping *m, MappingSet **ms) { return false; };
		protected:
			int type;
	};


	class LocalFunctionReference : public Transformation {
		public:
			LocalFunctionReference(string PointerName, int scope, TemplateMolecule *tm);
			virtual ~LocalFunctionReference() {};
			virtual void apply(Mapping *m, MappingSet **ms) {};
			virtual int getComponentIndex() const { return -1; };

			TemplateMolecule *getTemplateObject() const {return tm;};
			string getPointerName() const { return PointerName; };
			int getFunctionScope() const {return scope; };

			static const unsigned int SPECIES_FUNCTION=0;
			static const unsigned int SINGLE_MOLECULE_FUNCTION=1;
		protected:
			string PointerName;
			int scope;
			TemplateMolecule *tm;
	};



	class EmptyTransform : public Transformation {
		public:
			EmptyTransform() : Transformation(TransformationFactory::EMPTY){ this->cIndex=-1; };
			EmptyTransform(int cIndex) : Transformation(TransformationFactory::EMPTY){ this->cIndex=cIndex; };
			virtual ~EmptyTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms) {};
			virtual int getComponentIndex() const { return cIndex; };
		protected:
			int cIndex;
	};


	class StateChangeTransform : public Transformation {
		public:
			StateChangeTransform(int cIndex, int newValue);
			virtual ~StateChangeTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {return cIndex;};
		protected:
			int cIndex;
			int newValue;
	};

	class BindingTransform : public Transformation {
		public:
			BindingTransform(int cIndex, int otherReactantIndex, int otherMappingIndex);
			virtual ~BindingTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {return cIndex;};

			virtual bool checkForNullCondition(Mapping *m, MappingSet **ms);
		protected:
			int cIndex;
			int otherReactantIndex;
			int otherMappingIndex;
	};

	class NewMoleculeBindingTransform : public BindingTransform {
			public:
				NewMoleculeBindingTransform(int cIndex, int otherReactantIndex, int otherMappingIndex);
				virtual ~NewMoleculeBindingTransform() {};
				virtual bool checkForNullCondition(Mapping *m, MappingSet **ms);
	};

	/*! Deprecated!!!  molecularity is now checked at a more basic level with checkForNullCondition!! */
//	class BindingSeparateComplexTransform : public BindingTransform {
//			public:
//				BindingSeparateComplexTransform(int cIndex, int otherReactantIndex, int otherMappingIndex) :
//					BindingTransform(cIndex, otherReactantIndex, otherMappingIndex) {};
//				virtual ~BindingSeparateComplexTransform() {};
//				virtual void apply(Mapping *m, MappingSet **ms);
//				virtual int getComponentIndex() const {return cIndex;};
//	};

	class UnbindingTransform : public Transformation {
		public:
			UnbindingTransform(int cIndex);
			virtual ~UnbindingTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {return cIndex;};
		protected:
			int cIndex;
	};


	class AddSpeciesTransform : public Transformation {
		public:
			AddSpeciesTransform( SpeciesCreator * sc );
			virtual ~AddSpeciesTransform();
			virtual void apply( Mapping *m, MappingSet **ms );
			virtual int getComponentIndex() const {cerr<<"You should not get a component index from an AddMoleculeTransform!!"<<endl; return -1;};

		protected:
			SpeciesCreator * sc;
	};


	class AddMoleculeTransform : public Transformation
	{
		public:
			AddMoleculeTransform( MoleculeCreator * _mc );
			virtual ~AddMoleculeTransform();
			virtual void apply( Mapping * m, MappingSet ** ms ) { cerr<<"apply method should not be called from an AddMoleculeTranform!!"<<endl;};
			virtual int getComponentIndex() const { cerr<<"You should not get a component index from an AddMoleculeTransform!!"<<endl; return -1;};

			// adds molecule and points mapping set to that new molecule
			void apply_and_map( MappingSet * ms );
			// is this a population type?
			bool isPopulationType() const;
			// get pointer to population molecule
			Molecule * get_population_pointer() const;

		protected:
			MoleculeCreator * mc;
			Molecule *        new_molecule;
	};

	class RemoveMoleculeTransform : public Transformation {
		public:
			RemoveMoleculeTransform(int removalType);
			virtual ~RemoveMoleculeTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {cout<<"You should not get a component index from a RemoveMoleculeTransform!!"<<endl; exit(1); return -1;};
			virtual int getRemovalType() { return removalType; };

		protected:
			int removalType;

	};


	class DecrementPopulationTransform : public Transformation
	{
		public:
			DecrementPopulationTransform();
			virtual ~DecrementPopulationTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const { return cIndex; };
		protected:
			int cIndex;
	};


	class IncrementStateTransform : public Transformation {
		public:
			IncrementStateTransform(unsigned int stateIndex);
			virtual ~IncrementStateTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {return cIndex;};
		protected:
			int cIndex;
	};

	class DecrementStateTransform : public Transformation {
		public:
			DecrementStateTransform(unsigned int stateIndex);
			virtual ~DecrementStateTransform() {};
			virtual void apply(Mapping *m, MappingSet **ms);
			virtual int getComponentIndex() const {return cIndex;};
		protected:
			int cIndex;
	};



}







#endif /*TRANSFORMATION_HH_*/
