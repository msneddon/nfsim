#include "NFutil.hh"



using namespace NFutil;










void NFutil::trim(string& str) {
		string::size_type startpos = str.find_first_not_of(" \t");
		string::size_type endpos = str.find_last_not_of(" \t");
		if(( string::npos == startpos ) || (string::npos == endpos)) {
			str="";
		} else {
			str = str.substr(startpos,endpos-startpos+1);
		}
	}
