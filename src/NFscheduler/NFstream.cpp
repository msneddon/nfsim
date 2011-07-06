
#include "NFstream.h"

NFstream::NFstream() 
{
    check_mpi();
}

NFstream::NFstream(const char* filename, ios_base::openmode mode) 
{
    check_mpi();
    open(filename, mode);
}

NFstream::NFstream(bool useFile):useFile_(useFile) {}

NFstream::NFstream(string filename):useFile_(false)
{
    strname_ = filename;
}

NFstream::~NFstream() {}

void NFstream::check_mpi()
{
    useFile_ = true;
#ifdef NF_MPI
    useFile_ = false;
#endif
}

void NFstream::setUseFile(bool useFile)
{
    useFile_ = useFile;
}

string NFstream::getStrName()
{
    return strname_;
}

void NFstream::open(const char* filename, ios_base::openmode mode) 
{
    if (useFile_)
	file_.open(filename, mode); 
    else {
	strname_ = filename;
    }
}

void NFstream::close() 
{
    file_.close();
}

ostream& NFstream::write(const char* s , streamsize n)
{
    if (useFile_)
	return file_.write(s, n);
    else 
	return str_.write(s, n);
}

ostream& NFstream::flush()
{
    if (useFile_)
	return file_.flush();
    else
	return str_.flush();
}

ios_base::fmtflags NFstream::setf(ios_base::fmtflags fmtfl) 
{
    if (useFile_)
	return file_.setf(fmtfl);
    else 
	return str_.setf(fmtfl);
}

ios_base::fmtflags NFstream::setf(ios_base::fmtflags fmtfl, ios_base::fmtflags mask) 
{
    if (useFile_)
	return file_.setf(fmtfl, mask);
    else 
	return str_.setf(fmtfl, mask);
}

streamsize NFstream::precision() const
{
    if (useFile_)
	return file_.precision();
    else
	return str_.precision();
}

streamsize NFstream::precision(streamsize prec)
{
    if (useFile_)
	return file_.precision(prec);
    else
	return str_.precision(prec);
}
    
string NFstream::str() const 
{
    return str_.str();
}

bool NFstream::is_open()
{
    if (useFile_)
	return file_.is_open();
    else 
	return true;
}

void NFstream::test() 
{
    NFstream f1;
    
#ifdef NF_MPI
    cout << "NF_MPI defined\n";
#endif

//     ofstream o1("regular-ofstream.txt");
//     o1 << "works" << endl;

    string sss = "<test string ok>\n";
    f1.open("test.out");
    f1 << endl;
    f1 << "Hello NFstream. " << endl << sss << "ofstream seems to be working.." << endl;

    NFstream s1(false);
    s1 << endl;
    s1 << "Hello NFstream. " << endl << "stringstream seems to be working.." << endl;

    string s = s1.str();
    cout << "content of stringstream s1:" << endl << s;
    
}

NFstream& NFstream::operator<<(NFstream& (*func)(NFstream &)) 
{
    return ((*func))(*this);
}

// friend functions
template<class T>
NFstream& operator<<(NFstream& nfstream, const T& value) 
{
    if (nfstream.useFile_) 
	nfstream.file_ << value;
    else
	nfstream.str_ << value;

    return nfstream;
}

NFstream& endl(NFstream& nfstream) 
{
    if (nfstream.useFile_)
	nfstream.file_ << endl;
    else
	nfstream.str_ << endl;
   
    return nfstream;
}

