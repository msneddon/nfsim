#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#ifdef NF_MPI
#include "mpi.h"
#endif

#include "NFstream.h"
#include "../NFsim.hh"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define MASTER 0
#define TAG_MSG  99   
#define TAG_DATA 98
#define MSG_DATA_SIZE (1<<20)
#define MAX_MPI_SIZE  4096

enum {cmd_free=0, cmd_job, cmd_pre_data_ack, cmd_data_ack, rpt_ready, rpt_pre_data, rpt_data, rpt_done};

struct job {
	string filename;
	int processors;
	vector<string> argument;
	vector<string> argval;
	vector<string> parameters;
	vector<double> values;
};

struct scan {
	vector<string> parameter;
	vector<double> min;
	vector<double> max;
	vector<int> steps;
};

struct model {
	string filename;
	int processors;
	int replicates;
	vector<string> argument;
	vector<string> argval;
};

typedef struct msgtype {
	int src;
	int tag;
	int len;	
	char data[MSG_DATA_SIZE];
} msgtype;

vector<job*> parseJobsFile (string filename);

int schedulerInterpreter(int* argc, char*** argv);

void convertModelScanToJobs(model*& currentModel, scan* currentScan, vector<job*>& joblist);

const char* itoa(int inNum);

const char* dtoa(double inNum);

string getFileLine(ifstream &input);

vector<string>* stringToStrings(string fullString, const char* delim, bool treatConsecutiveDelimAsOne = true);

vector<string>* getStringsFileline(ifstream &input, const char* delim, bool treatConsecutiveDelimAsOne);

void findandreplace(string &source, string find, string replace);

void printFileLineOutput();

void push_stream(int rank, NFstream& strm);

void send_to_slave(int slave, int tag, int datalen, char *data);

void send_to_master(int myid, int tag, int datalen, char *data);

void recv_from_slave();

void recv_from_master();

void perr(const char*);

void printParallelJobOutput(vector<job*> jobQueue); 

void FinalizeMPI();

void InitializeMPI(int* argc, char*** argv,int& Size,int& Rank);

//Loads a file into a buffer
string load_to_buffer(string filename);

void DynamicParallel(map<string, string> argMap,int rank,int size);

void EmbarrassingParallel(map<string, string> argMap,int rank,int size);

string BroadcastString(int Rank,int From,string InBuffer);

string ConvergeAllData(int Rank,int Size,string Buffer);

void ConvertStringToBufferMap(map<string, map<int, string> >& FileMap,string ReportBuffer);

string ConvertBufferMapToString(map<string, map<int,string> >& FileBuffers);

char* ConvertStringToCString(string Buffer);

void PrintFileBuffer(map<string, map<int, string> > FileMap,vector<job*> JobQueue);

string getPath(string Filename);

#endif /* SCHEDULER_H_ */
