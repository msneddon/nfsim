#ifndef TRANSFORMATION_HH_
#define TRANSFORMATION_HH_


#include "../NFreactions.hh"

namespace NFcore
{

	class SpeciesCreator;

	//!  Keeps information about how to perform one transform of a Molecule.
	/*!
	 	This is the base unit used to store information about transformations.  Based on
	 	its "type", it is able to store state changes, binding, unbinding, molecule creation,
	 	and molecule deletion.  It is managed and owned by a TransformationSet object which
	 	uses Transformation objects to actually carry out the operation.  Transforms are
	 	created by a TransformationSet which calls one of the static functions that generates
	 	new Transformations.  Transformation Objects exist for each transformation in the
	 	System and are not deleted until the simulation is over.
	    @author Michael Sneddon
	 */
	class Transformation
	{
		public:
			
			/*!
			 	Deletes this Transformation and clears the data associated with it.
			    @author Michael Sneddon
			 */
			~Transformation();
			
			/*!
			 	Returns the type of this transformation (for instance, state change, binding, etc).
			 	The list of possible types are static constant integers in this class, but in general
			 	you do not need to keep track of them as you can only create transformations
			 	via the static generator functions.
			    @author Michael Sneddon
			 */
			unsigned int getType() const { return type; };
			
			/*!
			 	Returns the State or Binding Site index of the component that is transformed.  This
			 	value may be irrelevant for certain types of transformations, and what it represents
			 	can be determined by the type of transformation. 
			    @author Michael Sneddon
			 */
			unsigned int getStateOrSiteIndex() const { return stateORsiteIndex; };
			
			/*!
			 	For a state change transformation, this will return the new value of the state
			 	that is transformed.  In other words, the state will change to this value.  For other
			 	transformations, this value is meaningless.
			    @author Michael Sneddon
			 */
			int getNewStateValue() const { return newStateValue; };
			
			/*!
			 	For Binding transformations, this remembers which reactant this Transformation object
			 	will bind to.  For any binding reaction, there are two Transformation objects created for
			 	each site.  This is the first bit of information needed to get the other site.  The other
			 	piece of information can be obtained from the getPartnerMappingIndes() function.
			    @author Michael Sneddon
			 */
			unsigned int getPartnerReactantIndex() const { return otherReactantIndex; };
			
			/*!
			 	For Binding transformations, this remembers which particular Mapping this Transformation
			 	object will bind to.  You will need this, along the the PartnerReactantIndex, (from the
			 	getPartnerReactantIndex() function) to completely specify the partner.
			    @author Michael Sneddon
			 */
			unsigned int getPartnerMappingIndex() const { return otherMappingIndex; };
			
			
			
			
			void createSpecies();
			
			
			
			
			
			/*!
			 	Generates a state change Transformation for transforming the given state at the given
			 	index into the new state value.
			    @author Michael Sneddon
			 */
			static Transformation * genStateChangeTransform(unsigned int stateIndex, int newStateValue);
			
			/*!
			 	Generates a binding transformation for one of the two binding sites in the binding Transform.  You will
			 	have to tell it where the other Transformation object lives (TransformationSet has this information).
			    @author Michael Sneddon
			 */
			static Transformation * genBindingTransform1(unsigned int bSiteIndex, unsigned int otherReactantIndex, unsigned int otherMappingIndex);
			
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
			 	Generates an Add Molecule transformation.  Currently, this is not yet
			 	implemented.
			    @author Michael Sneddon
			 */
			static Transformation * genAddMoleculeTransform(SpeciesCreator *sc);
			
			/*!
			 	Generates a removal of a molecule from the system.  Currently this is
			 	not yet implemented.
			    @author Michael Sneddon
			 */
			static Transformation * genRemoveMoleculeTransform();
			
			/*!
			 	Generates an empty transformation.  This is used in cases where there is
			 	a reactant that is not transformed in a reaction, but that still needs
			 	to be counted and marked so that the rate of the reaction is correct.
			    @author Michael Sneddon
			 */
			static Transformation * genEmptyTransform();
			
			
			
			
			
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
			static const unsigned int SKIP = 5; 
			
			
		protected:
			
			/*!
				Constructs a new Transformation object.  This is protected because Transformations
				should only be created by the static public generator functions of this class.  So
				don't use this unless you know exactly what you are doing!
				@author Michael Sneddon
			*/
			Transformation();
			
			/*!	Remembers the type of transformation this is	*/
			unsigned int type;
			
			/*!	Remember the index of the state or binding site that this transform affects	*/
			unsigned int stateORsiteIndex;
			
			//For a state change transformation
			/*!	Keeps the final state of a state change reaction	*/
			int newStateValue;
			
			//For a binding reaction
			/*!	Keeps the reactant index of the binding partner of this binding transform	*/
			unsigned int otherReactantIndex;
			/*!	Keeps the mapping index of the binding partner of this binding transform	*/
			unsigned int otherMappingIndex;
			
			//For a creation of a new molecule, not yet impelemented
			SpeciesCreator *sc;
	};
}







#endif /*TRANSFORMATION_HH_*/
