#ifndef REACTANTCONTAINER_HH_
#define REACTANTCONTAINER_HH_




#include "../NFreactions.hh"

namespace NFcore
{
	//Forward Declarations
	class TransformationSet;
	class MappingSet;
	class ReactantContainer;


	//!  Interface for all containers of Reactants (MappingSets) needed by ReactionClasses
	/*!
	 *
	 *
	 */
	class ReactantContainer {


		public:
			ReactantContainer() { hasClonedMappings=false; };
			virtual ~ReactantContainer() {};

			/*!
			 */
			virtual int size() const = 0;


			/*!
			 */
			virtual MappingSet * pushNextAvailableMappingSet() = 0;;

			/*!
			 */
			virtual void popLastMappingSet() = 0;;

			/*!
			 */
			virtual void removeMappingSet(unsigned int mappingSetId) = 0;;


			/*!
			 */
			virtual MappingSet * getMappingSet(unsigned int mappingSetId) const = 0;;

			/*!
			 */
			virtual void printDetails() const = 0;;

			/*!
			 */
			void notifyPresenceOfClonedMappings() { hasClonedMappings=true; };

			/*!
			 */
			bool getHasClonedMappings() const { return hasClonedMappings; };



		protected:
			bool hasClonedMappings;


	};
}













#endif /* REACTANTTREE_HH_ */
