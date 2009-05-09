#include "Scheduler.h"

#include "../NFsim.hh"

#include <iostream>
#include <map>

using namespace std;

msgtype msg;
vector<job*> jobQueue;

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

vector<job*> parseJobsFile (string filename) {
    //for test
    //return getTestJobs();

    vector<job*> joblist;
	
    ifstream input;
    input.open(filename.data());
    if (!input.is_open()) {
	cout << "could not open " << filename << endl;
	return joblist;
    }

    bool inBlock;
    int currentJobID;
    scan* currentScan = NULL;
    model* currentModel = NULL;
    while (!input.eof()) {
	vector<string>* strings = getStringsFileline(input, " ",true);
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
		    currentModel->start = 0;
		    currentModel->stop = 100;
		    currentModel->timestep = 0.01;
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
			} else if ((*subStrings)[0].compare("start") == 0 && (*subStrings)[1].length() >= 3) {
			    currentModel->start = atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
			} else if ((*subStrings)[0].compare("stop") == 0 && (*subStrings)[1].length() >= 3) {
			    currentModel->stop = atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
			} else if ((*subStrings)[0].compare("timestep") == 0 && (*subStrings)[1].length() >= 3) {
			    currentModel->timestep = atof((*subStrings)[1].substr(1,(*subStrings)[1].length()-2).data());
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
    input.close();

    return joblist;
}

void convertModelScanToJobs(model*& currentModel, scan* currentScan, vector<job*>& joblist) {
    vector<job*> currentJobVector;
    if (currentModel != NULL) {
	for (int i=0; i < currentModel->replicates; i++) {
	    job* newJob = new job;
	    newJob->filename = currentModel->filename;
	    newJob->processors = currentModel->processors;
	    newJob->start = currentModel->start;
	    newJob->stop = currentModel->stop;
	    newJob->timestep = currentModel->timestep;
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
			newJob->start = currentJobVector[k]->start;
			newJob->stop = currentJobVector[k]->stop;
			newJob->timestep = currentJobVector[k]->timestep;
			newJob->parameters = currentJobVector[k]->parameters;
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
    ostringstream strout;
    strout << inNum;
    return strout.str().data();
}

const char* dtoa(double inNum) {
    ostringstream strout;
    strout << inNum;
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

bool fetch_job(job* &pjob) {
    static vector<job*>::iterator it = jobQueue.begin();
    while (it != jobQueue.end()) {
	pjob = *it;
	++it;
	return true;
    }
    return false;
}

void job2str(job& j, char* p) {
    sprintf(p, "%s,", j.filename.c_str());	p += strlen(p);	
    sprintf(p, "%d,", j.processors);	    p += strlen(p);

    int argc = j.argument.size();
    sprintf(p, "%d,", argc); p += strlen(p);
    if (argc > 0) {
	for (int i = 0; i < argc; ++i) {
	    sprintf(p, "%s,",  j.argument[i].c_str()); p += strlen(p);
	    sprintf(p, "%s,",  j.argval[i]); p += strlen(p);
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
    string fname("slave." + string(itoa(rank)) + "." + strm.getStrName());
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
    int rank, size = 1;
#ifdef NF_MPI
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    if (size < 2) {
#ifdef NF_MPI
	cerr << "MPI: at least 2 nodes needed." << endl;
	MPI_Finalize();
#endif
	return -1;
    }

    map<string, string> argMap;
    NFinput::parseArguments(*argc, const_cast<const char**>(*argv), argMap);

    map<string,string>::iterator mapIt = argMap.find("jobfile");
    if (mapIt == argMap.end()) {
#ifdef NF_MPI
	MPI_Finalize();
#endif
	return 1;
    }
		
    printf("Hello, I am %d of %d.\n", rank, size);

    if (rank == 0) {
	master_init(size);
	// Calling the job file parser to get unrolled list of jobs
	jobQueue = parseJobsFile (mapIt->second);
	job* pjob;
	job jnow;
	int jcount = 0;
	bool done = false;
	bool slave_available = true;
	int  left = size - 1;	// # slaves still working
	while (left > 0) {
	    if (!done && slave_available) {
		done = !fetch_job(pjob);
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
	for (vector<job*>::iterator it = jobQueue.begin(); it != jobQueue.end(); ++it) {
	    cout << "filename = " << filenames[*it] << endl;
	    cout << "  output = " << endl << buffers[*it] << endl;
	}
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
	
#ifdef NF_MPI
    MPI_Finalize();
#endif
    return 0;
};

