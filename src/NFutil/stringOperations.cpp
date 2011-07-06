#include "NFutil.hh"



using namespace NFutil;










void NFutil::trim(string& str) {
		size_t startpos = str.find_first_not_of(" \t");
		size_t endpos = str.find_last_not_of(" \t");
		if(( string::npos == startpos ) || (string::npos == endpos)) {
			str="";
		} else {
			str = str.substr(startpos,endpos-startpos+1);
		}
	}
