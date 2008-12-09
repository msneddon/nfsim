#ifndef TRANSFORMATIONSET_HH_
#define TRANSFORMATIONSET_HH_


#include "../NFreactions.hh"

using namespace std;

namespace NFcore
{
	//Forward declarations
	class Transformation;
	class TemplateMolecule;
	class Molecule;
	class SpeciesCreator;


	//!  Maintains a set of Transformation objects for a ReactionClass
	/*!
	 	This class is one of the core pieces of a ReactionClass.  It maintains the
	 	set of Transformation objects for each of the reactants in the ReactionClass.
	 	Therefore, it knows about the TemplateMolecules of the ReactionClass; it knows
	 	how to generate MappingSet objects for the ReactantLists; it knows how to
	 	create MappingGenerator objects for TemplateMolecules; it knows how to execute
	 	the set of Transformations; and finally, it knows how to generate the list of
	 	product molecules that are generated from firing a reaction rule.  To use this
	 	class, first create your TemplateMolecule objects for your ReactionClass.  Create
	 	a TransformationSet from the TemplateMolecule vector (the vector should be created
	 	such that each element of the vector is a reactant in the ReactionClass).  Add
	 	as many Transformations to the set as you like using the TransformationSet static
	 	functions.  Finally, call the finalize() function to finish setting upt the
	 	TransformationSet.  Then you can use this TransformationSet to create a ReactionClass
	 	with ease.  Really, this is all easier than it sounds! (see simple_system.hh and
	 	simple_system.cpp )
	    @author Michael Sneddon
	 */
	class TransformationSet
	{
		public:

			/*!
			 	Creates a new TransformationSet from the given vector of TemplateMolecules.
			 	This takes care of setting everything up in the TemplateMolecules (such as
			 	adding MappingGenerators).  Be sure to call finalize() before using this
			 	in a ReactionClass!
				@author Michael Sneddon
			*/
			TransformationSet(vector <TemplateMolecule *> reactantTemplates);

			/*!
				Destroys the TransformationSet and associated Transformation objects.
				@author Michael Sneddon
			*/
			~TransformationSet();

			/*!
				Adds a state change transformation on the given TemplateMolecule (that must have been included
				in the original vector of TemplateMolecules) along with the stateName and final value
				of the state to be transformed.
				@author Michael Sneddon
			*/
			bool addStateChangeTransform(TemplateMolecule *t, string cName, int finalStateValue);

			/*!
				Adds a state change transformation on the given TemplateMolecule (that must have been included
				in the original vector of TemplateMolecules) along with the stateName and final value
				of the state to be transformed.
				@author Michael Sneddon
			*/
			bool addStateChangeTransform(TemplateMolecule *t, string cName, string finalStateValue);


			bool addIncrementStateTransform(TemplateMolecule *t, string cName);
			bool addDecrementStateTransform(TemplateMolecule *t, string cName);


			/*!
				Adds a binding reaction between the two given TemplateMolecules at the specified
				binding sites.
				@author Michael Sneddon
			*/
			bool addBindingTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2);

			/*!
				Adds a binding reaction between the two given TemplateMolecules at the specified
				binding sites with the constraint that the two molecules are not connected.  Note: This only stops
				the binding transform!  It does not prevent the entire reaction!  That is not programmed in yet!
				@author Michael Sneddon
			*/
			bool addBindingSeparateComplexTransform(TemplateMolecule *t1, string bSiteName1, TemplateMolecule *t2, string bSiteName2);



			/*!
				Adds an unbinding reaction at the given site of the given TemplateMolecule
				@author Michael Sneddon
			*/
			bool addUnbindingTransform(TemplateMolecule *t, string bSiteName, TemplateMolecule *t2, string bSiteName2);


			/*!
				Adds a delete rule to the given TemplateMolecule.  This has not been
				implemented yet, so it may take other parameters.
				@author Michael Sneddon
			*/
			bool addDeleteMolecule(TemplateMolecule *t);

			/*!
				Adds a create molecule rule, but has not been implemented yet.
				@author Michael Sneddon
			*/
			bool addAddMolecule(SpeciesCreator *sc);


			/*!
				Adds a create molecule rule, but has not been implemented yet.
				@author Michael Sneddon
			*/
			bool addLocalFunctionReference(TemplateMolecule *t, string PointerName, int scope);


			/*!
				Call this (in a ReactionClass) to transform the array of MappingSets (one
				MappingSet per reactant in the correct position in the array, please!).
				@author Michael Sneddon
			*/
			bool transform(MappingSet **mappingSets);

			/*!
				Generates a blank MappingSet (blank in the sense that it is not mapped to
				any Molecules yet) from the list of Transformations.  This function is
				called by ReactantList and ReactantTree to populate the lists of MappingSets
				that are created at the beginning of a simulation.
				@author Michael Sneddon
			*/
			MappingSet *generateBlankMappingSet(unsigned int reactantIndex, unsigned int mappingSetId);


			/*!
				Returns the TemplateMolecule of a given reactant.  This function is primarily
				used for initial initializations of a ReactionClass.
				@author Michael Sneddon
			*/
			TemplateMolecule * getTemplateMolecule(unsigned int reactantIndex) const { return reactants[reactantIndex]; };

			/*!
				Get the number of reactants in the rule governed by this TransformationSet.  This
				function is really only used for initial initializations of a ReactionClass.
				@author Michael Sneddon
			*/
			unsigned int getNreactants() const {return n_reactants; };

			/*!
				This function sets up the actual array of Transformation objects for this TransformationSet
				and makes sure that no one adds additional Transformations once this TransformationSet is
				in use.  TransformationSets also do not allow you to do certain things (like add it
				to a ReactionClass) until the Set is finalized by a call to this function.
				@author Michael Sneddon
			*/
			void finalize();

			/*!
				Check whether or not the finalize() function has been called.
				@author Michael Sneddon
			*/
			bool isFinalized() const { return finalized; };

			/*!
				From an array of MappingSet objects, this function populates the list of Molecules
				with all the Molecules that can be affected by this Transformation.  The traversalLimit
				parameters sets the depth at which the Molecules should be searched (starting at the
				original TemplateMolecules given in the initial vector of TemplateMolecules).  So giving
				a value of 1 will only give you the immediate reactant Molecules.  Setting this to two
				will explore down one level of bonds, and so on.
				@author Michael Sneddon
			*/
			bool getListOfProducts(MappingSet **mappingSets, list<Molecule *> &products, int traversalLimit);

			bool hasSymUnbindingTransform() const { return hasSymUnbinding; };


			int getNumOfTransformations(int reactantIndex) const { return transformations[reactantIndex].size();};
			Transformation *getTransformation(int reactantIndex, int index) const { return transformations[reactantIndex].at(index); };




		protected:

			/*!
				Used for error checking when setting up a TransformationSet.  This finds a TemplateMolecule
				in the TemplateMolecule vector given in the constructor and finds the reactant index under
				which it was found.  This makes sure a particular TemplateMolecule that you are trying to
				add a transformation to exists and does not exist in multiple places.
				@author Michael Sneddon
			*/
			int find(TemplateMolecule *t);

			/*!	Remembers if the finalize function has been called	*/
			bool finalized;

			/*!	Keeps track of the number of reactants in the rule monitored by this TransformationSet	*/
			unsigned int n_reactants;

			/*!	The array of TemplateMolecules that represent the reactants		*/
			TemplateMolecule ** reactants;

			/*!	A vector that holds the actual Transformation objects	*/
			vector <Transformation *> *transformations;


			vector <Transformation *> addMoleculeTransformations;
			bool hasSymUnbinding;
	};

}












#endif /*TRANSFORMATIONSET_HH_*/
