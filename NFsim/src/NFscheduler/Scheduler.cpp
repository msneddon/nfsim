#include "Scheduler.h"

#include "../NFsim.hh"

#include <iostream>
#include <map>

using namespace std;

msgtype msg;

job job1, job2, jobs[10];

// globals for slave
vector<string> slave_filenames;
vector<string> slave_buffers;

// globals for master
vector<char*> raw_buffers;
vector<int> incoming_sizes;
vector<job*> slave_assignment;

map<job*, string> filenames;
map<job*, string> buffers;

vector<job*> getTestJobs() {
	vector<job*> joblist;

	job1.filename = "motor.xml";
	job1.processors = 2;
	job1.argument.push_back("eqTime");
	job1.argval.push_back("0");
	job1.argument.push_back("sTime");
	job1.argval.push_back("10");
	job1.argument.push_back("oSteps");
	job1.argval.push_back("10");
	job1.parameters.push_back("motorCount");
	job1.values.push_back(12);
	job1.parameters.push_back("cellVolume");
	job1.values.push_back(1.72e-15);

	job2.filename = "simple_system.xml";
	job2.processors = 3;
	job2.argument.push_back("eqTime");
	job2.argval.push_back("10");
	job2.argument.push_back("sTime");
	job2.argval.push_back("20");
	job2.argument.push_back("oSteps");
	job2.argval.push_back("25");
	job2.parameters.push_back("kcat");
	job2.values.push_back(0.3);

	for (int i = 0; i < 5; ++i) {
	jobs[2*i] = job1;
	jobs[2*i+1] = job2;
	joblist.push_back(&jobs[2*i]);
	joblist.push_back(&jobs[2*i+1]);
	}

	joblist.push_back(&job1);
	joblist.push_back(&job2);
	
	return joblist;
}

vector<job*> parseJobsFile (string buffer) {
	//for test
	//return getTestJobs();

	vector<job*> joblist;
	vector<string>* lines = stringToStrings(buffer,"\n",true);

	bool inBlock;
	int currentJobID;
	scan* currentScan = NULL;
	model* currentModel = NULL;
	for (int m=0; m < int(lines->size()); m++) {
		vector<string>* strings = stringToStrings((*lines)[m],"\t ",true);
		for (int i=0; i < int(strings->size()); i++) {
			//checking for block start
			if ((*strings)[i].length() > 0 && (*strings)[i].substr(0,1).compare("<") == 0) {
				inBlock = true;
				//determining the block type
				if ((*strings)[i].length() >= 4 && (*strings)[i].substr(0,4).compare("<job") == 0) {
					currentJobID = 0;
				} else if ((*strings)[i].length() >= 6 && (*strings)[i].substr(0,6).compare("<model") == 0) {
					if (currentModel != NULL) {
						convertModelScanToJobs(currentModel,currentScan,joblist);
					}
					currentModel = new model;
					currentModel->processors = 1;
					currentModel->replicates = 1;
				} else if ((*strings)[i].length() >= 5 && (*strings)[i].substr(0,5).compare("<scan") == 0) {
					if (currentScan == NULL) {
						currentScan = new scan;
					}
				} else if ((*strings)[i].length() >= 6 && (*strings)[i].substr(0,6).compare("</job>") == 0) {					
					currentJobID = -1;
				} else if ((*strings)[i].length() >= 7 && (*strings)[i].substr(0,7).compare("</scan>") == 0) {
					if (currentScan->parameter.size() > 1) {
						currentScan->parameter.pop_back();
						currentScan->min.pop_back();
						currentScan->max.pop_back();
						currentScan->steps.pop_back();
					} else {
						delete currentScan;
						currentScan = NULL;
					}
				}
			} else if ((*strings)[i].find("=") != -1) {
				//checking for block stop
				if ((*strings)[i].substr(((*strings)[i].length()-1),1).compare(">") == 0) {
					inBlock = false;
					(*strings)[i] = (*strings)[i].substr(0,(*strings)[i].length()-1);
				}
						
				//breaking up the string into a parameter and value pair
				vector<string>* subStrings = stringToStrings((*strings)[i],"=");
				if (subStrings->size() >= 2) {
					if (currentScan != NULL && currentModel == NULL) {
						if ((*subStrings)[0].compare("param") == 0 && (*subStrings)[1].length() >= 3) {
							currentScan->parameter.push_back((*subStrings)[1].substr(1,(*subStrings)[1].length()-2));
						} else if ((*subStrings)[0].compare("min") == 0 && (*subStrings)[1].length() >= 3) {
							currentScan->min.push_back(atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data()));
						} else if ((*subStrings)[0].compare("max") == 0 && (*subStrings)[1].length() >= 3) {
							currentScan->max.push_back(atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data()));
						} else if ((*subStrings)[0].compare("steps") == 0 && (*subStrings)[1].length() >= 3) {
							int steps = atoi((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
							if (steps <= 2) {
							steps = 2;
							}
							currentScan->steps.push_back(steps);
						} else if ((*subStrings)[0].compare("stepsize") == 0 && (*subStrings)[1].length() >= 3) {
							int steps = 2;
							double stepsize = atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
							if (stepsize > 0) {
							steps = 1+int((currentScan->max[currentScan->max.size()-1] - currentScan->min[currentScan->min.size()-1])/stepsize);
							if (steps < 2) {
								steps = 2;
							}
							}
							currentScan->steps.push_back(steps);
						}
					} else if (currentModel != NULL) {
						if ((*subStrings)[0].compare("file") == 0 && (*subStrings)[1].length() >= 3) {
							currentModel->filename = (*subStrings)[1].substr(1,(*subStrings)[1].length()-2);
						} else if ((*subStrings)[0].compare("procs") == 0 && (*subStrings)[1].length() >= 3) {
							currentModel->processors = atoi((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
						} else if ((*subStrings)[0].compare("replicates") == 0 && (*subStrings)[1].length() >= 3) {
							currentModel->replicates = atoi((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
						} else {
							currentModel->argument.push_back((*subStrings)[0]);
							currentModel->argval.push_back((*subStrings)[1].substr(1,(*subStrings)[1].length()-2));
						}
					} else if (currentJobID != -1) {
						if ((*subStrings)[0].compare("id") == 0 && (*subStrings)[1].length() >= 3) {
							currentJobID = atoi((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
						}
					}
				}

				//if the block has ended we create the jobs in the job vector
				if  (!inBlock && currentModel != NULL) {
					convertModelScanToJobs(currentModel,currentScan,joblist);
				}
			} else if ((*strings)[i].compare(">") == 0) {
				inBlock = false;
				if  (!inBlock && currentModel != NULL) {
					convertModelScanToJobs(currentModel,currentScan,joblist);
				}
			}
		}
		delete strings;
	}

	return joblist;
}

void convertModelScanToJobs(model*& currentModel, scan* currentScan, vector<job*>& joblist) {
	vector<job*> currentJobVector;
	if (currentModel != NULL) {
		for (int i=0; i < currentModel->replicates; i++) {
			job* newJob = new job;
			newJob->filename = currentModel->filename;
			newJob->processors = currentModel->processors;
			newJob->argument = currentModel->argument;
			newJob->argval = currentModel->argval;
			currentJobVector.push_back(newJob);
		}
		if (currentScan != NULL) {
			for (int j=int(currentScan->parameter.size()-1); j >= 0; j--) {
				vector<job*> newJobVector;
				for (int i=0; i < currentScan->steps[j]; i++) {
					for (int k=0; k < int(currentJobVector.size()); k++) {
						job* newJob = new job;
						newJob->filename = currentJobVector[k]->filename;
						newJob->processors = currentJobVector[k]->processors;
						newJob->parameters = currentJobVector[k]->parameters;
						newJob->argument = currentJobVector[k]->argument;
						newJob->argval = currentJobVector[k]->argval;
						newJob->values = currentJobVector[k]->values;
						newJob->parameters.push_back(currentScan->parameter[j]);
						newJob->values.push_back(currentScan->min[j] + i*(currentScan->max[j]-currentScan->min[j])/(currentScan->steps[j]-1));
						newJobVector.push_back(newJob);
					}
				}
				for (int k=0; k < int(currentJobVector.size()); k++) {
					delete currentJobVector[k];
				}
				currentJobVector.clear();
				for (int i=0; i < int(newJobVector.size()); i++) {
					currentJobVector.push_back(newJobVector[i]);
				}
			}
		}
	}
	for (int k=0; k < int(currentJobVector.size()); k++) {
		joblist.push_back(currentJobVector[k]);
	}
	delete currentModel;
	currentModel = NULL;
}

string getFileLine(ifstream &input) {
	string buff; 
	getline( input, buff );
	return buff;
}

vector<string>* getStringsFileline(ifstream &input, const char* delim, bool treatConsecutiveDelimAsOne) {
	string buff = getFileLine(input);
	return stringToStrings(buff, delim, treatConsecutiveDelimAsOne);
}

vector<string>* stringToStrings(string fullString, const char* delim, bool treatConsecutiveDelimAsOne) {
	vector<string>* newVect = new vector<string>;
	string buff(fullString);

	int location;
	do {
	location = int(buff.find_first_of(delim));
	if (location != -1) {
		if (location == 0) {
		if (!treatConsecutiveDelimAsOne) {
			string newString;
			newVect->push_back(newString);
		}
		buff = buff.substr(location+1, buff.length()-(location+1));
		} else {
		string newString = buff.substr(0, location);
		newVect->push_back(newString);
		buff = buff.substr(location+1, buff.length()-(location+1));
		}
	}
	} while(location != -1);
	
	if (buff.length() != 0 || !treatConsecutiveDelimAsOne) {
	newVect->push_back(buff);
	}
	
	return newVect;
}

void findandreplace(string &source, string find, string replace) {
	size_t j;
	for (;(j = source.find( find )) != source.npos;) {
	source.replace( j, find.length(), replace );
	}
}

const char* itoa(int inNum) {
	stringstream strout;
	strout << inNum;
	strout.flush();
	return strout.str().data();
}

const char* dtoa(double inNum) {
	stringstream strout;
	strout << inNum;
	strout.flush();
	return strout.str().data();
}

// MPI communication routines
void send_to_slave(int slave, int tag, int datalen, char *data) {
	msg.src = MASTER;
	msg.tag = tag;
	msg.len = datalen;
	if (datalen > 0 && data != 0) memcpy(msg.data, data, datalen);
	int actlen = sizeof(msg.src) + sizeof(msg.tag) + sizeof(msg.len) + datalen;
#ifdef NF_MPI
	MPI_Send(&msg, actlen, MPI_CHAR, slave, TAG_MSG, MPI_COMM_WORLD);	
#endif
}

void send_to_master(int myid, int tag, int datalen, char *data) {
	if (tag == rpt_data) {
#ifdef NF_MPI
	MPI_Send(data, datalen, MPI_CHAR, MASTER, TAG_DATA, MPI_COMM_WORLD);
#endif	
	} else {
	msg.src = myid;
	msg.tag = tag;
	msg.len = datalen;
	if (datalen > 0 && data != 0) memcpy(msg.data, data, datalen);
	int actlen = sizeof(msg.src) + sizeof(msg.tag) + sizeof(msg.len) + datalen;
#ifdef NF_MPI
	MPI_Send(&msg, actlen, MPI_CHAR, MASTER, TAG_MSG, MPI_COMM_WORLD);	
#endif
	}
}

void recv_from_slave() {
#ifdef NF_MPI
	MPI_Status status;
	MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (status.MPI_TAG == TAG_DATA) {
	int id = status.MPI_SOURCE;
	msg.src = id;
	msg.tag = rpt_data;
	MPI_Recv(raw_buffers[id], incoming_sizes[id], MPI_CHAR, MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, &status);
	} else {
	MPI_Recv(&msg, sizeof(msg), MPI_CHAR, MPI_ANY_SOURCE, TAG_MSG, MPI_COMM_WORLD, &status);
	}
#endif
}

void recv_from_master() {
#ifdef NF_MPI
	MPI_Status status;
	MPI_Recv(&msg, sizeof(msg), MPI_CHAR, MASTER, TAG_MSG, MPI_COMM_WORLD, &status);
#endif
}

void job2str(job& j, char* p) {
	sprintf(p, "%s,", j.filename.c_str());	p += strlen(p);	
	sprintf(p, "%d,", j.processors);	    p += strlen(p);

	int argc = j.argument.size();
	sprintf(p, "%d,", argc); p += strlen(p);
	if (argc > 0) {
	for (int i = 0; i < argc; ++i) {
		sprintf(p, "%s,",  j.argument[i].c_str()); p += strlen(p);
		sprintf(p, "%s,",  j.argval[i].c_str()); p += strlen(p);
	}
	}

	int n = j.parameters.size();
	sprintf(p, "%d,", n); p += strlen(p);
	if (n > 0) {
	for (int i = 0; i < n; ++i) {
		sprintf(p, "%s,",  j.parameters[i].c_str()); p += strlen(p);
		sprintf(p, "%lg,", j.values[i]); p += strlen(p);
	}
	}
}

void str2job(char* str, job& jnow) {
	char *p = str;
	char *ch = strtok(str, ",");
	ch = strtok(0, ","); jnow.filename = string(p);    p = ch; 
	ch = strtok(0, ","); jnow.processors = atoi(p);    p = ch; 
	ch = strtok(0, ","); int argc        = atoi(p);    p = ch;
	for (int i = 0; i < argc; ++i) { 
	ch = strtok(0, ","); jnow.argument.push_back(string(p)); p = ch; 
	ch = strtok(0, ","); jnow.argval.push_back(string(p));   p = ch; 
	}
	ch = strtok(0, ","); int n = atoi(p); p = ch;
	for (int i = 0; i < n; ++i) { 
	ch = strtok(0, ","); jnow.parameters.push_back(string(p));        p = ch; 
	ch = strtok(0, ","); jnow.values.push_back(atof((const char*)p)); p = ch; 
	}
}

void push_stream(int rank, NFstream& strm) {
	string fname(strm.getStrName());
	slave_filenames.push_back(fname);
	slave_buffers.push_back(strm.str());
}

void perr(const char* message) {
	fprintf(stderr, "%s\n", message);
}

void clear_slave_data() {
	slave_filenames.clear();
	slave_buffers.clear();
}

void slave_work(int rank, job& jnow) {
	clear_slave_data();

	bool verbose = true;
	map<string, string> argMap;

	argMap["xml"] = jnow.filename;
	System *s = initSystemFromFlags(argMap, verbose);
//     System* s = NFinput::initializeFromXML(jnow.filename, true, true);

	for (int i = 0; i < jnow.argument.size(); ++i) {
		argMap[jnow.argument[i]] = jnow.argval[i];
	}
	for (int i = 0; i < jnow.parameters.size(); ++i) {
		s->addParameter(jnow.parameters[i], jnow.values[i]);
	}

	runFromArgs(s, argMap, verbose);

	NFstream& strm = s->getOutputFileStream();
	push_stream(rank, strm);
	
// 	s->prepareForSimulation();
// 	s->updateSystemWithNewparameters();	

// 	double eqTime = 0;
// 	double sTime  = 10;
// 	int oSteps = 10;
	
// 	s->equilibrate(eqTime);
// 	s->sim(sTime, oSteps);
}

void master_init(int size) {
	raw_buffers.resize(size);
	incoming_sizes.resize(size);
	slave_assignment.resize(size);
}

int schedulerInterpreter(int* argc, char*** argv) {
	//Initializing problem
	int rank = 0;
	int size = 1;

	InitializeMPI(argc,argv,size,rank);

	//Parsing arguments
	map<string, string> argMap;
	NFinput::parseArguments(*argc, const_cast<const char**>(*argv), argMap);

	//Checking that a job file has been provided
	if (argMap.count("jobfile") == 0) {
		FinalizeMPI();
		return 1;
	}
	
	//Calling the appropriate parallel processing algorithm
	if (size  == 1 || argMap.count("embarrassing") > 0) {
		EmbarrassingParallel(argMap,rank,size);
	} else {
		DynamicParallel(argMap,rank,size);
	}

	FinalizeMPI();
	return 0;
};

void printParallelJobOutput(vector<job*> jobQueue) {
	//First organizing all output buffers from all jobs by filename
	map<string, map<job*, string> > FileMap;
	for (vector<job*>::iterator it = jobQueue.begin(); it != jobQueue.end(); ++it) {
		FileMap[filenames[*it]][*it] = buffers[*it];
	}
	
}

void FinalizeMPI() {
	#ifdef NF_MPI
	MPI_Finalize();
	#endif
}

void InitializeMPI(int* argc, char*** argv,int& Size,int& Rank) {
	Size = 1;
	Rank = 0;

	#ifdef NF_MPI
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	#endif
}

string load_to_buffer(string filename) {
	string buffer;
	ifstream Input;
	Input.open(filename.data());
	if (!Input.is_open()) {
		cout << "Could not open " << filename << endl;
		return "";
	}
	while(!Input.eof()) {
	    if (buffer.length() > 0) {
			buffer.append("\n");
	    }
	    buffer.append(getFileLine(Input));
	}
	Input.close();
	return buffer;
}

void DynamicParallel (map<string, string> argMap,int rank,int size) {
	printf("Hello, I am %d of %d.\n", rank, size);

	vector<job*> jobQueue;
	if (rank == 0) {
		int CurrentJob = 0;
		master_init(size);
		// Calling the job file parser to get unrolled list of jobs
		jobQueue = parseJobsFile (argMap["jobfile"]);
		job* pjob;
		job jnow;
		int jcount = 0;
		bool done = false;
		bool slave_available = true;
		int  left = size - 1;	// # slaves still working
		while (left > 0) {
			if (!done && slave_available) {
				if (jobQueue.size() <= jcount) {
					done = true;
				}
				pjob = jobQueue[jcount];
				jnow = job(*pjob);
				if (!done) {
					++jcount;
					printf("master: fetched job #%d\n", jcount);
				}
			}
				recv_from_slave();
				slave_available = false;
				if (msg.tag == rpt_ready || msg.tag == rpt_done) {
				slave_available = true;
				if (!done) {
					printf("master: assigning work #%d to slave #%d \n", jcount, msg.src);
					char str[MSG_DATA_SIZE];
					job2str(jnow, str);
					slave_assignment[msg.src] = pjob;
					send_to_slave(msg.src, cmd_job, strlen(str)+1, str);
				} else {
					--left;
					send_to_slave(msg.src, cmd_free, 0, 0);
				}
				} else if (msg.tag == rpt_pre_data) { // msg.data = "data_size,filename"
				char *p = strchr(msg.data, ','); 
				*(p++) = 0; 
				job *j = slave_assignment[msg.src];
				filenames[j] = p;
				int data_size = atoi(msg.data);
				raw_buffers[msg.src] = (char*)malloc(data_size);
				incoming_sizes[msg.src] = data_size;
				send_to_slave(msg.src, cmd_pre_data_ack, 0, 0);
				} else if (msg.tag == rpt_data) {
				job* j = slave_assignment[msg.src];
				buffers[j] = raw_buffers[msg.src];
				free(raw_buffers[msg.src]);
				send_to_slave(msg.src, cmd_data_ack, 0, 0);
			}
		}
		// filenames & buffers ready to be processed
		// code for parsing output goes here
		cout << endl << "---------------- output (process me)---------------" << endl << endl;
		printParallelJobOutput(jobQueue);
	} else {
		// slave 
		send_to_master(rank, rpt_ready, 0, 0);
		while (1) {
			recv_from_master();
			if (msg.tag == cmd_free) {
				printf("slave #%d : free now\n", rank);
				break;
			}			
			job jnow;
			char str[MSG_DATA_SIZE];
			memcpy(str, msg.data, msg.len);	
			printf("slave #%d : got work (%s)\n", rank, str);
			str2job(str, jnow);
			slave_work(rank, jnow);

			for (int i = 0; i < slave_filenames.size(); ++i) {
				sprintf(str, "%d,%s", slave_buffers[i].length()+1, slave_filenames[i].c_str());
				send_to_master(rank, rpt_pre_data, strlen(str)+1, str);
				recv_from_master();
				if (msg.tag != cmd_pre_data_ack) perr("Error: expecting cmd_pre_data_ack");
				send_to_master(rank, rpt_data, slave_buffers[i].length()+1, const_cast<char*>(slave_buffers[i].c_str()));
				recv_from_master();
				if (msg.tag != cmd_data_ack) perr("Error: expecting cmd_data_ack");
			}
			send_to_master(rank, rpt_done, 0, 0);
		}
	}
}

void EmbarrassingParallel(map<string, string> argMap,int rank,int size) {
	//Master node reads job buffer and broadcasts it to all other nodes
	vector<job*> jobQueue;
	string JobBuffer;
	if (rank == 0) {
		JobBuffer = load_to_buffer(argMap["jobfile"]);
		if (size > 1) {
			BroadcastString(rank,0,JobBuffer);
		}
	} else {
		JobBuffer = BroadcastString(rank,0,"");
	}

	//Parsing the jobs buffer
	jobQueue = parseJobsFile(JobBuffer);
	
	//Every processor runs job based on rank
	int Processor = 0;
	map<string, map<int,string> > FileBuffers;
	for (int i=0; i < int(jobQueue.size()); i++) {
		if (rank == Processor) {
			bool verbose = false;
			map<string, string> CurrentArgs = argMap;
			CurrentArgs["xml"] = jobQueue[i]->filename;
			for (int j = 0; j < jobQueue[i]->argument.size(); ++j) {
				CurrentArgs[jobQueue[i]->argument[j]] = jobQueue[i]->argval[j];
			}
			System *s = initSystemFromFlags(CurrentArgs, verbose);
			for (int j = 0; j < jobQueue[i]->parameters.size(); ++j) {
				s->addParameter(jobQueue[i]->parameters[j], jobQueue[i]->values[j]);
			}
			s->getOutputFileStream().setUseFile(false);
			runFromArgs(s, CurrentArgs, verbose);
			FileBuffers[s->getOutputFileStream().getStrName()][i].assign(s->getOutputFileStream().str());
			delete s;
		}
		Processor++;
		if (Processor >= size) {
			Processor = 0;
		}
	}

	//Gathering all results from all nodes
	if (size > 1 && rank != 0) {
		//Converting buffer map to single string
		string ReportBuffer = ConvertBufferMapToString(FileBuffers);
		//Transmitting all buffers to the master node
		cout << "Sending results" << endl;
		ConvergeAllData(rank,size,ReportBuffer);
		cout << "Results sent" << endl;
	} else if (size > 1) {
		cout << "Receiving results" << endl;
		string ReportBuffer = ConvergeAllData(rank,size,"");
		ConvertStringToBufferMap(FileBuffers,ReportBuffer);
		cout << "Results received" << endl;
	}

	//Printing results from master node
	if (rank == 0) {
		PrintFileBuffer(FileBuffers,jobQueue);
	}
}

string BroadcastString(int Rank,int From,string InBuffer) {
	#ifdef NF_MPI
	int Length;
	if (Rank == From) {
		Length = InBuffer.length();
	}
	MPI_Bcast(&Length, 1, MPI_INT, From, MPI_COMM_WORLD);

	if (Length > 0) {
		char* Buffer;
		if (Rank == From) {
			Buffer = new char[InBuffer.length()+1];		
			strcpy(Buffer,InBuffer.data());
		} else {
			Buffer = new char[Length];
		}
		MPI_Bcast(Buffer, Length, MPI_CHAR, From, MPI_COMM_WORLD);
		InBuffer.assign(Buffer);
		delete [] Buffer;
	}
	#endif
	cout << InBuffer << endl;
	return InBuffer;
}

string ConvertBufferMapToString(map<string, map<int, string> >& FileMap) {
	stringstream Result;
	Result << FileMap.size();
	for (map<string, map<int, string> >::iterator MapIT = FileMap.begin(); MapIT != FileMap.end(); MapIT++) {
		Result << "`" << MapIT->first << "`" << MapIT->second.size();
		for (map<int, string>::iterator MapITT = FileMap[MapIT->first].begin(); MapITT != FileMap[MapIT->first].end(); MapITT++) {
			Result << "`" << MapITT->first << "`" << MapITT->second;
		}
	}
	Result.flush();
	return Result.str();
}

//Transmitting all buffers to the master node
string ConvergeAllData(int Rank,int Size,string Buffer) {
	int OtherNode = 1;
	int CurrentRank = Rank;
	int CurrentMessageSize = int(Buffer.length());
	char* CurrentMessage = ConvertStringToCString(Buffer);
	int CurrentMultiplier = 1;
	bool Done = false;
	while (!Done) {
		if (CurrentRank % 2 == 0) {
			OtherNode = CurrentMultiplier*(CurrentRank+1);
			if (OtherNode >= Size) {
				//Indicating when the master node is done
				if (CurrentRank == 0) {
					Done = true;
				}
			} else {			
				#ifdef NF_MPI
				MPI_Status status;
				int MessageSize;
				MPI_Recv(&MessageSize,1, MPI_INT, OtherNode, TAG_DATA, MPI_COMM_WORLD, &status);
				if (MessageSize > 0) {
					char* OldMessage = CurrentMessage;
					CurrentMessage = new char[CurrentMessageSize+MessageSize];
					MPI_Recv(CurrentMessage,MessageSize, MPI_CHAR, OtherNode, TAG_DATA, MPI_COMM_WORLD, &status);
					if (OldMessage != NULL) {
						for (int i=MessageSize; i < (MessageSize+CurrentMessageSize); i++) {
							CurrentMessage[i] = OldMessage[i-MessageSize];
						}
						delete [] OldMessage;
					}
					CurrentMessageSize += MessageSize;
				}
				#endif
			}
		} else {
			OtherNode = CurrentMultiplier*(CurrentRank-1);
			if (!Done) {
				#ifdef NF_MPI
				MPI_Send(&CurrentMessageSize,1, MPI_INT, OtherNode, TAG_DATA, MPI_COMM_WORLD);
				if (CurrentMessageSize > 0) {
					MPI_Send(CurrentMessage,CurrentMessageSize, MPI_CHAR, OtherNode, TAG_DATA, MPI_COMM_WORLD);
				}
				Done = true;
				delete [] CurrentMessage;
				CurrentMessage = NULL;
				#endif
			}
		}
		//Doubling the current multiplier
		CurrentMultiplier = CurrentMultiplier*2;
		//Changing the rank
		CurrentRank = CurrentRank/2;
	}

	if (Rank == 0) {
		Buffer.assign(CurrentMessage);
		delete [] CurrentMessage;
	}
	return Buffer;
}


void ConvertStringToBufferMap(map<string, map<int, string> >& FileMap,string ReportBuffer) {
	vector<string>* Strings = stringToStrings(ReportBuffer,"`",false);
	int CurrentIndex = 0;
	while (CurrentIndex <= int(Strings->size()-5)) {
		int NumFiles = atoi((*Strings)[CurrentIndex].data());
		CurrentIndex++;
		for (int i=0;i < NumFiles; i++) {
			string CurrentFile = (*Strings)[CurrentIndex];
			CurrentIndex++;
			int Jobs = atoi((*Strings)[CurrentIndex].data());
			CurrentIndex++;
			for (int j=0; j < Jobs; j++) {
				int CurrentJob = atoi((*Strings)[CurrentIndex].data());
				FileMap[CurrentFile][CurrentJob] = (*Strings)[CurrentIndex+1];
				CurrentIndex += 2;
			}
		}
	}
	delete Strings;
}

char* ConvertStringToCString(string Buffer) {
	char* CString = NULL;
	if (Buffer.length() > 0) {
		CString = new char[Buffer.length()+1];
		strcpy(CString,Buffer.data());
	}
	return CString;
}

void PrintFileBuffer(map<string, map<int, string> > FileMap,vector<job*> JobQueue) {
	//Now going through each distinct filename and printing output in a variety of formats
	for (map<string, map<int, string> >::iterator it = FileMap.begin(); it != FileMap.end(); ++it) {
		string Filename = getPath(JobQueue[it->second.begin()->first]->filename)+it->first;
		//if (Filename.length() <= 4 || Filename.substr(Filename.length()-4).compare(".gdat") != 0) {
			//Appending the output buffer from each job 
			ofstream Output;
			Output.open(Filename.data());
			for (map<int, string>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt) {
				//Printing job data
				job* CurrentJob = JobQueue[itt->first];
				Output << "**JOB SPECIFICATIONS**" << endl; 
				Output << "JOBOVERVIEW(MODEL:" << CurrentJob->filename;
				for (int i=0; i < int(CurrentJob->argument.size()); i++) {
					Output << "," << CurrentJob->argument[i] << ":" << CurrentJob->argval[i];
				}
				Output << ")" << endl;
				if (CurrentJob->parameters.size() > 0) {
					Output << "PARAMETER;VALUE" << endl;
					for (int i=0; i < int(CurrentJob->parameters.size()); ++i) {
						Output << CurrentJob->parameters[i] << ";" << CurrentJob->values[i] << endl; 
					}
				}
				Output << "**JOB OUTPUT**" << endl;
				//Printing file data
				Output << itt->second << endl << endl;
			}
			Output.close();
		//}
	}
}

string getPath(string InFilename) {
	int Position = int(InFilename.rfind("/"));
	if (Position < int(InFilename.rfind("\\"))) {
		Position = int(InFilename.rfind("\\"));
	}

	if (Position != -1) {
		return InFilename.substr(0,Position+1);
	}
	return "";
}
