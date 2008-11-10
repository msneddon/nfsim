/*
 * NFoutput.hh
 *
 *  Created on: Sep 28, 2008
 *      Author: msneddon
 */
#ifndef NFOUTPUT_HH_
#define NFOUTPUT_HH_



#include "../NFcore/NFcore.hh"

namespace NFcore
{
	class MoleculeType;
	class System;

	class Outputter {
		public:
			Outputter(string filename, System *s);
			virtual ~Outputter();

			virtual void output()=0;
			virtual void outputHeader()=0;
			virtual void tryToDump(double simTime)=0;
		protected:
			string filename;
			System *s;
			ofstream outputFileStream;
	};

	class DumpMoleculeType : public Outputter {
		public:
			DumpMoleculeType(string filename, System *s, MoleculeType *mt);
			virtual ~DumpMoleculeType();

			virtual void output();
			virtual void outputHeader();
		protected:
			MoleculeType *mt;
	};

	class DumpAllComplexes : public Outputter {
		public:
			DumpAllComplexes(string filename, System *s);
			virtual ~DumpAllComplexes();

			virtual void output();
			virtual void outputHeader();
			virtual void tryToDump(double simTime);
		protected:
			MoleculeType *mt;
	};

}








#endif /* NFOUTPUT_HH_ */
