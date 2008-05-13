#include "NFutil.hh"



using namespace NFutil;



	
	double NFutil::convertToDouble(const std::string& s)
	{
		bool failIfLeftoverChars = true;
		std::istringstream i(s);
		double x;
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
			throw std::runtime_error("error in NFutil::convertToDouble(\"" + s + "\")");
		return x;
	}
	
	
	
	
	int NFutil::convertToInt(const std::string& s)
	{
		bool failIfLeftoverChars = true;
		std::istringstream i(s);
		int x;
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
			throw std::runtime_error("error in NFutil::convertToInt(\"" + s + "\")");
		return x;
	}