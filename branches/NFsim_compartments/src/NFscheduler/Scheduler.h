#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "mpi.h"

#include "../NFsim.hh"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define MASTER 0
#define TAG    99             // universal MPI tag
#define MSG_DATA_SIZE (1<<20)
#define MAX_MPI_SIZE  4096

enum {cmd_free=0, rpt_ready, rpt_done};

struct job {
	string filename;
	int processors;
	double start;
	double stop;
	double timestep;
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
	double start;
	double stop;
	double timestep;
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

void send_to_slave(int slave, int tag, int datalen, char *data);

void send_to_master(int myid, int tag, int datalen, char *data);

void recv_from_slave();

void recv_from_master();


#endif /* SCHEDULER_H_ */
