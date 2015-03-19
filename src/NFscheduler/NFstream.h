
#ifndef _NFSTREAM_H_
#define _NFSTREAM_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class NFstream
{
public:
    NFstream();
    NFstream(const char* filename, ios_base::openmode mode = ios_base::out);
    NFstream(bool useFile);
    NFstream(string filename);

    ~NFstream();

    void setUseFile(bool useFile);
    string getStrName();

    void open(const char* filename, ios_base::openmode mode = ios_base::out);
    void close();
    ostream& write(const char* s , streamsize n);
    ostream& flush();
    bool is_open();
    string str() const;

    ios_base::fmtflags setf(ios_base::fmtflags fmtfl);
    ios_base::fmtflags setf(ios_base::fmtflags fmtfl, ios_base::fmtflags mask);
    streamsize precision() const;
    streamsize precision(streamsize prec);

    static void test();

    NFstream& operator<<(NFstream& (*func)(NFstream &));

    friend NFstream& endl (NFstream& nfstream);

    template<class T>
    friend NFstream& operator<< (NFstream& nfstream, const T& value);
    
private:
    ofstream file_;
    stringstream str_;

    bool useFile_;
    string strname_;

    void check_mpi();
};

template<class T>
NFstream& operator<< (NFstream& nfstream, const T& value);

NFstream& endl (NFstream& nfstream);

#endif /* _NFSTREAM_H_ */
