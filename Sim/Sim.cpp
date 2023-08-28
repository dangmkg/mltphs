// Sim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <intrin.h>
#include <windows.h>

#include <cmath>

#include <map>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <thread>
#include <atomic>
#include <random>

import Conf;
import SimCore;
import RNG;
import Json;
import SortData;

class DataPoint {
public:
    double lambda[8];
    double mu[8];
    size_t queueSize[8];
    double minDelta;
    double meanServTime, meanTimeInSys;
    size_t dropPktCnt;
};

class Worker {
    std::string filename;
    time_t startTime, endTime;
    int nDataPoints, nNodes, nPackets, curPoint;
    double gamma;
    std::thread* threadObj;
    std::atomic_bool done;
    void generateInputData(size_t i, int limit, double scale, double shift, int limitDelta, double scaleDelta, double shiftDelta);
    void mainFunc();
public:
    std::vector<DataPoint> data;
    int id;

    Worker(int id, std::string filename, int nDataPoints, int nNodes, int nPackets, time_t startTime, double gamma);
    void setParams(std::string filename, int nDataPoints, int nNodes, int nPackets, time_t startTime, double gamma);
    int getCurrentProgress();
    void run();
    void writeToFile();
    void waitTillDone();
    bool isFree();
};

Worker::Worker(int id, std::string filename, int nDataPoints, int nNodes, int nPackets, time_t startTime, double gamma) {
    this->id = id;
    this->done.store(true);
    this->setParams(filename, nDataPoints, nNodes, nPackets, startTime, gamma);
}

int Worker::getCurrentProgress() {
    return curPoint;
}

void Worker::setParams(std::string filename, int nDataPoints, int nNodes, int nPackets, time_t startTime, double gamma) {
    if (filename != "") this->filename = filename;
    else {
        std::string rngtype;
        if (Conf::rngType == Conf::MersenneTwister) rngtype = "mt19937";
        else if (Conf::rngType == Conf::RDRAND) rngtype = "rdrand";
        this->filename = std::format(
            "{}_n{}_p{}_g{}_{}_w{}_d{}.json",
            rngtype,
            nNodes,
            nPackets,
            (int)(gamma * 10),
            startTime,
            this->id,
            nDataPoints
        );
    }
    this->nDataPoints = nDataPoints;
    this->nNodes = nNodes;
    this->nPackets = nPackets;
    this->startTime = startTime;
    this->gamma = gamma;
}

void Worker::run() {
    this->done.store(false);
    threadObj = new std::thread(&Worker::mainFunc, this);
}

void Worker::writeToFile() {
    std::string key;
    std::ofstream out;
    Json::Writer json;
    int i, j;

    out.open(Conf::workingPath / "result" / filename);
    json.setOutput(&out);
    json.openObject(1, "");
    json.writeInt("nData", nDataPoints);
    json.writeInt("nNodes", nNodes);
    json.writeInt("timeRunning", endTime - startTime);

    for (i = 0; i < nNodes; i++) {
        key = "lambda" + std::to_string(i);
        json.openArray(0, key);
        for (j = 0; j < nDataPoints; j++)
            json.writeDouble("", data[j].lambda[i]);
        json.closeArray();

        key = "mu" + std::to_string(i);
        json.openArray(0, key);
        for (j = 0; j < nDataPoints; j++)
            json.writeDouble("", data[j].mu[i]);
        json.closeArray();
    }

    json.openArray(0, "meanServTime");
    for (i = 0; i < nDataPoints; i++)
        json.writeDouble("", data[i].meanServTime);
    json.closeArray();

    json.openArray(0, "meanTimeInSys");
    for (i = 0; i < nDataPoints; i++)
        json.writeDouble("", data[i].meanTimeInSys);
    json.closeArray();

    json.closeObject();
    out.close();
}

unsigned __int64 generateUniform(unsigned __int64 limit) {
    unsigned __int64 result;
    _rdrand64_step(&result);
    return (result % limit);
}

void Worker::generateInputData(size_t i, int limit, double scale, double shift, int limitDelta, double scaleDelta, double shiftDelta) {
    size_t j, minIdx;
    std::vector<double> lambda, mu;
    double sumLambda, delta, minDelta;

    sumLambda = 0.0;
    lambda.resize(nNodes); mu.resize(nNodes);
    do minDelta = (double)generateUniform(limitDelta);
    while (minDelta == 0.0);
    minIdx = generateUniform(nNodes);
    for (j = 0; j < nNodes; j++) {
        if (j != minIdx) {
            do delta = (double)generateUniform(minDelta + 1);
            while (delta == 0.0);
        }
        else {
            delta = minDelta;
        }
        delta = -delta / double(limitDelta) * scaleDelta + shiftDelta;
        delta = exp(delta);

        do lambda[j] = (double)generateUniform(limit) + shift;
        while (lambda[j] == 0.0);
        sumLambda = sumLambda + lambda[j];

        mu[j] = round(sumLambda / (1 - delta));

        data[i].lambda[j] = lambda[j] / double(limit) * scale;
        data[i].mu[j] = mu[j] / double(limit) * scale;

        if (j == 0) data[i].minDelta = delta;
        else if (data[i].minDelta > delta)
            data[i].minDelta = delta;
    }
}

void Worker::mainFunc() {
    size_t i, j;
    double sumFactors, normFactor;
    std::vector<size_t> queueLen;
    std::string msg, bkfn;
    Json::Writer json;
    std::ofstream bk, bkresult;
    
    if (Conf::workingMode == "interactive") {
        msg = std::format("Worker {} starts\n", id);
        std::cout << msg;
    }

    /*
    bkfn = "bk_"; bkfn.append(filename);
    bk.open(Conf::workingPath / "backup" / bkfn);
    
    json.setOutput(&bk);
    json.openObject(1, "");
    json.writeInt("startTime", startTime);
    json.writeInt("nDataPoints", nDataPoints);
    json.writeInt("nNodes", nNodes);
    json.writeInt("nPackets", nPackets);
    json.closeObject();
    bk << "\n";
    */

    
    this->data.clear();
    this->data.resize(nDataPoints);
    sumFactors = 0.0;
    for (i = 0; i < nDataPoints; i++) {
        generateInputData(i, 16384, 16.0, 0, 16384, 2.0, -4.0);
        //generateInputData(i, 16384, 16.0, 16384, 16384, 4.0, -2.0);
        //generateInputData(i, 65536, 32.0);
        sumFactors = sumFactors + pow(gamma, -log(data[i].minDelta));

        /*
        json.openObject(0, "");
        json.writeArrayOfDouble("lambda", data[i].lambda, nNodes);
        json.writeArrayOfDouble("mu", data[i].mu, nNodes);
        json.writeDouble("minDelta", data[i].minDelta);
        json.closeObject();
        bk << "\n";
        */
    }
    normFactor = ((double)nDataPoints) / sumFactors;

    //bk.close();
    /*filename = "backup\\result";
    filename.append(std::to_string(startTime));
    filename.append("_w"); filename.append(std::to_string(id));
    filename.append(".json");
    bkresult.open(filename);
    json.setOutput(&bkresult);*/

    SimCore::Sim sim;
    queueLen.resize(nNodes);
    std::fill_n(queueLen.begin(), nNodes, nPackets);
    double servTime, timeInSys;

    for (i = 0; i < nDataPoints; i++) {
        if (i % 1000 == 0 && Conf::workingMode == "interactive") {
            msg = std::format("Worker {} is at iteration {}\n", id, i);
            std::cout << msg;
        }
        else if (Conf::workingMode == "batch") {
            this->curPoint = i;
        }

        sim.init(nNodes, data[i].lambda, data[i].mu, queueLen.data());
        double k = pow(gamma, -log(data[i].minDelta)) * normFactor;
        sim.runSim((size_t)(nPackets * k));

        servTime = 0.0;
        timeInSys = 0.0;
        for (j = 0; j < sim.dataCnt; j++) {
            servTime = servTime + sim.pkt[j].totalServTime;
            timeInSys = timeInSys + sim.pkt[j].totalServTime + sim.pkt[j].totalQueueTime;
        }
        data[i].meanServTime = servTime / sim.dataCnt;
        data[i].meanTimeInSys = timeInSys / sim.dataCnt;

        /*json.openScope(0, 0, "");
        json.writeArrayOfDouble("lambda", data[i].lambda, nNodes);
        json.writeArrayOfDouble("mu", data[i].mu, nNodes);
        json.writeDouble("minDelta", data[i].minDelta);
        json.writeDouble("dataCnt", sim.dataCnt);
        json.writeDouble("meanServTime", data[i].meanServTime);
        json.writeDouble("meanTimeInSys", data[i].meanTimeInSys);
        json.closeScope();
        bkresult << "\n";*/

        sim.cleanUp();
    }

    //bkresult.close();
    std::time(&this->endTime);
    writeToFile();
    this->done.store(true);
}

void Worker::waitTillDone() {
    threadObj->join();
}

bool Worker::isFree() {
    return this->done.load();
}

class JobInfo {
public:
    int nData;
    int nNodes;
    int nPackets;
    double gamma;
    bool inProcessing = false;
    Worker* worker = nullptr;

    int linenumber;

    void loadFromJson(Json::Object* job, bool& error) {
        if (job->numfields.contains("nData"))
            nData = job->numfields["nData"];
        else {
            error = true;
        }
        if (job->numfields.contains("nNodes"))
            nNodes = job->numfields["nNodes"];
        else {
            error = true;
        }
        if (job->numfields.contains("nPackets"))
            nPackets = job->numfields["nPackets"];
        else {
            error = true;
        }
        if (job->numfields.contains("gamma"))
            gamma = job->numfields["gamma"];
        else {
            gamma = Conf::gamma;
        }
    }
};

class BatchMode {
private:
    int currentLine;
    COORD oldCursor;

    void saveCursorPos() {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        CONSOLE_SCREEN_BUFFER_INFO scrBufInfo;
        GetConsoleScreenBufferInfo(hStdout, &scrBufInfo);
        oldCursor = scrBufInfo.dwCursorPosition;
    }

    void restoreCursorPos() {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleCursorPosition(hStdout, oldCursor);
    }

    void changeCursorPosTo(int linenumber) {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        CONSOLE_SCREEN_BUFFER_INFO scrBufInfo;
        GetConsoleScreenBufferInfo(hStdout, &scrBufInfo);
        scrBufInfo.dwCursorPosition.Y = scrBufInfo.dwCursorPosition.Y + linenumber - currentLine;
        scrBufInfo.dwCursorPosition.X = 0;
        SetConsoleCursorPosition(hStdout, scrBufInfo.dwCursorPosition);
    }

    void startbatchjobs(Json::Object* json) {
        int nWorkers = 0;
        size_t i, j, nJobsDone;
        std::vector<Worker*> worker;
        std::string fname, key, errormsg;
        Json::Object* jobs, * job;
        std::vector<JobInfo> jobinfo;
        bool error, foundFreeWorker;

        time_t startTime = std::time(0);

        if (json->objfields.contains("jobs") && json->objfields["jobs"]->type == Json::array) {
            if (json->numfields.contains("nWorkers"))
                nWorkers = json->numfields["nWorkers"];
            if (nWorkers <= 0) nWorkers = 1;
            worker.resize(nWorkers);
            for (i = 0; i < worker.size(); i++) {
                worker[i] = new Worker(i, "", 0, 0, 0, 0, 0);
            }

            jobs = json->objfields["jobs"];
            jobinfo.resize(jobs->idx.size());
            for (i = 0; i < jobs->idx.size(); i++) {
                error = false;
                errormsg = std::format("Job item {} is invalid\n", i + 1);
                key = std::to_string(i);
                if (jobs->objfields.contains(key)) {
                    job = jobs->objfields[key];
                    jobinfo[i].loadFromJson(job, error);
                    if (error) std::cout << errormsg;
                }
                else {
                    std::cout << errormsg; error = true;
                }
            }
        }
        else {
            jobs = nullptr; error = true;
            std::cout << "Invalid batch file\n";
        }

        nJobsDone = 0;
        while (!error && jobs != nullptr && nJobsDone < jobs->idx.size()) {
            for (i = 0; i < jobs->idx.size(); i++) {
                if (jobinfo[i].inProcessing) {
                    if (jobinfo[i].worker->isFree()) {
                        jobinfo[i].inProcessing = false;
                        nJobsDone++;
                    }
                    else {
                        saveCursorPos();
                        changeCursorPosTo(jobinfo[i].linenumber);
                        int curPoint = jobinfo[i].worker->getCurrentProgress();
                        std::cout << std::format("Job {} with d:{} n:{} p:{} g:{} - current progress {}\n",
                            i + 1, jobinfo[i].nData, jobinfo[i].nNodes, jobinfo[i].nPackets, (int)(jobinfo[i].gamma * 10), curPoint);
                        restoreCursorPos();
                    }
                }
                else if (jobinfo[i].worker == nullptr) {
                    j = 0; foundFreeWorker = false;
                    while (!foundFreeWorker && j < worker.size()) {
                        if (worker[j]->isFree()) foundFreeWorker = true;
                        else j++;
                    }
                    if (foundFreeWorker) {
                        fname = std::format("Job{}_", i + 1);
                        fname = std::format(
                            "{}_n{}_p{}_g{}_{}_Job{}_d{}.json",
                            "mt19937",
                            jobinfo[i].nNodes,
                            jobinfo[i].nPackets,
                            (int)(jobinfo[i].gamma * 10),
                            startTime,
                            i + 1,
                            jobinfo[i].nData
                        );

                        jobinfo[i].inProcessing = true;
                        jobinfo[i].worker = worker[j];
                        worker[j]->setParams(fname, jobinfo[i].nData, jobinfo[i].nNodes, jobinfo[i].nPackets, startTime, jobinfo[i].gamma);
                        worker[j]->run();
                        jobinfo[i].linenumber = currentLine; currentLine++;
                        std::cout << std::format("Job {} with d:{} n:{} p:{} g:{} - current progress 0\n",
                            i + 1, jobinfo[i].nData, jobinfo[i].nNodes, jobinfo[i].nPackets, (int)(jobinfo[i].gamma * 10));
                    }
                }              
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        }
    }
public:
    BatchMode(int argc, char* argv[]) {
        std::filesystem::path path;
        std::ifstream file;
        Json::Reader jsonReader;
        Json::Object* json = nullptr;
        bool error = false;

        if (argc >= 2) {
            path.assign(argv[1]);
        }
        else if (Conf::batchFile.compare("") != 0) {
            path.assign(Conf::batchFile);
        }
        else {
            std::cout << "No input parameter!\n";
            error = true;
        }
        if (!error) {
            if (std::filesystem::exists(path)) {
                file.open(path);
                jsonReader.setInput(&file);
                json = jsonReader.parse();
                file.close();
                if (json != nullptr) startbatchjobs(json);
                else std::cout << "Invalid json input file";
            }
            else std::cout << "Input file not found\n";
        }
    }

    BatchMode() {
        std::string fname;
        std::ifstream file;
        Json::Reader jsonReader;
        Json::Object* json = nullptr;

        std::cout << "Enter batch filename: ";
        std::cin >> fname;
        if (std::filesystem::exists(fname)) {
            file.open(fname);
            jsonReader.setInput(&file);
            json = jsonReader.parse();
            file.close();
            if (json != nullptr) startbatchjobs(json);
            else std::cout << "Invalid json input file";
        }
        else std::cout << "Input file not found\n";
    }
};

namespace SimMM1k {
    class DataPoint {
    public:
        double lambda[128];
        double mu[128];
        size_t k[128];

        int nNodes;
        unsigned __int64 seedsLambda[128];
        unsigned __int64 seedsMu[128];
    };

    class Worker {
        std::string filename;
        time_t startTime, endTime;
        int nPackets, curPoint;
        std::thread* threadObj;
        std::atomic_bool done;
        void mainFunc();
        double eps;
    public:
        DataPoint* data;
        size_t nData;
        int id;

        Worker(int id, DataPoint* data, size_t nData);
        int getCurrentProgress();
        void run(double eps, std::string fname);
        void writeToFile();
        void waitTillDone();
        bool isFree();
    };

    void generateInputDataWithFixedK(DataPoint* &data, int nNodes, int k, size_t nData) {
        size_t i, j, limit;
        double scale, lambda, rho;

        data = new DataPoint[nData];
        limit = 65536;
        scale = 10.0;

        for (i = 0; i < nData; i++) {
            lambda = (double)generateUniform(limit) + 1;
            for (j = 0; j < nNodes; j++) {
                rho = (double)generateUniform(limit) - (limit / 2);
                rho = rho / limit * 7.0;
                rho = exp(rho);

                data[i].mu[j] = lambda * rho / double(limit) * scale;
                data[i].lambda[j] = 0.0;
                data[i].k[j] = k;
            }
            data[i].lambda[0] = lambda / double(limit) * scale;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputDataWithFixedKAndRatio(DataPoint*& data, int nNodes, int k, double ratioK0, double rhoScale, size_t nData) {
        size_t i, j, limit;
        double scale, lambda, rho;

        data = new DataPoint[nData];
        limit = 65536;
        scale = 100.0;

        for (i = 0; i < nData; i++) {
            lambda = (double)generateUniform(limit) + 1;
            for (j = 0; j < nNodes; j++) {
                rho = (double)generateUniform(limit) - (limit / 2);
                rho = rho / limit * rhoScale * 2;
                rho = exp(rho);

                data[i].mu[j] = lambda * rho / double(limit) * scale;
                data[i].lambda[j] = 0.0;
                data[i].k[j] = k;
            }
            data[i].k[0] = round(k * ratioK0);
            data[i].lambda[0] = lambda / double(limit) * scale;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputData(DataPoint*& data, size_t nData) {
        size_t i, j, limit;
        double scale, lambda, rho;

        int nNodes = 3;
        int k = 10;

        data = new DataPoint[nData];

        for (i = 0; i < nData; i++) {
            for (j = 0; j < nNodes; j++) {
                data[i].k[j] = k;
            }
            data[i].lambda[0] = 6.09756;
            data[i].mu[0] = 181.953;
            data[i].mu[1] = 19.0229;
            data[i].mu[2] = 0.220227;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputData2(DataPoint*& data, size_t nNodes, double lambda, double mu, size_t nData) {
        size_t i, j, limit;

        data = new DataPoint[nData];

        for (i = 0; i < nData; i++) {
            for (j = 0; j < nNodes; j++) {
                data[i].k[j] = i;
                data[i].mu[j] = mu;
            }
            data[i].lambda[0] = lambda;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputData3(DataPoint*& data, size_t nNodes, double lambda, double mu, int k, size_t nData) {
        size_t i, j, limit;

        data = new DataPoint[nData];

        for (i = 0; i < nData; i++) {
            for (j = 0; j < nNodes; j++) {
                data[i].k[j] = k;
                data[i].mu[j] = mu;
            }
            data[i].k[0] = i;
            data[i].lambda[0] = lambda;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputData4(DataPoint*& data, size_t nNodes, double lambda, double mu, int k, int kOrbit, size_t nData) {
        size_t i, j, limit;

        data = new DataPoint[nData];

        for (i = 0; i < nData; i++) {
            for (j = 1; j < nNodes; j++) {
                data[i].k[j] = k;
                data[i].mu[j] = mu;
            }
            data[i].k[0] = kOrbit;
            data[i].mu[0] = i+1;
            data[i].lambda[0] = lambda;
            data[i].nNodes = nNodes;
        }
    }

    void generateInputData5(DataPoint*& data, size_t nNodes, double lambda, double mu, size_t nData) {
        size_t i, j, limit;

        data = new DataPoint[nData];

        for (i = 0; i < nData; i++) {
            for (j = 0; j < nNodes; j++) {
                data[i].k[j] = i;
                data[i].mu[j] = mu;
            }
            data[i].lambda[0] = lambda;
            data[i].mu[0] = 0;
            data[i].nNodes = nNodes;
        }
    }

    Worker::Worker(int id, DataPoint* data, size_t nData) {
        this->id = id;
        this->done.store(true);
        this->data = data;
        this->nData = nData;
    }

    int Worker::getCurrentProgress() {
        return curPoint;
    }

    void Worker::run(double eps, std::string fname) {
        this->eps = eps;
        this->filename = fname;
        this->done.store(false);
        threadObj = new std::thread(&Worker::mainFunc, this);
    }

    void Worker::mainFunc() {
        size_t i, j;
        std::string msg, bkfn;
        Json::Writer json;
        time_t start, end;
        SimCore::Sim sim;
        std::ofstream f;

        f.open(this->filename);
        json.setOutput(&f);
        json.openObject(1, "");
        json.openArray(1, "data");
        for (i = 0; i < nData; i++) {
            std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
            //sim.init2(data[i].nNodes, data[i].lambda, data[i].mu, data[i].k);
            sim.init3(data[i].nNodes, data[i].lambda[0], data[i].mu, data[i].k);
            sim.runSim2(this->eps);
            std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

            std::cout << "Datapoint " << i << " avg outflow: " << sim.nPktAtSink / sim.lastTimeLenChanged << "\n";
            json.openObject(1, "");
            json.writeInt("nNodes", data[i].nNodes);
            json.openArray(0, "k");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", data[i].k[j]);
            }
            json.closeArray();
            json.writeDouble("lambd", data[i].lambda[0]);
            json.openArray(0, "mu");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", data[i].mu[j]);
            }
            json.closeArray();
            json.writeDouble("avgTime", sim.avgTimeInSys);
            json.writeDouble("avgLen", sim.avgLen / sim.lastTimeLenChanged);
            json.writeDouble("dropRate", (double)sim.nDrop / (double)sim.nPkt);
            std::chrono::duration<double> diff = end - start;
            json.writeDouble("simTime", diff.count());
            json.writeDouble("simTotalPkt", sim.nPkt);

            json.openArray(0, "avgLenNode");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", sim.node[j].avgLen / sim.lastTimeLenChanged);
            }
            json.closeArray();
            json.openArray(0, "lambdOut");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", sim.node[j].outFlow);
            }
            json.closeArray();
            //json.writeDouble("avgRetryCnt", sim.avgRetryCnt/sim.nPktAtSink);
            json.openArray(0, "avgTimeNode");
            for (j = 0; j < data[i].nNodes; j++) {
                if ((double)sim.node[j].nOutPkt > 0)
                    json.writeDouble("", sim.node[j].totalPktTime / (double)sim.node[j].nOutPkt);
                else
                    json.writeDouble("", 0);
            }
            json.closeArray();
            /*json.openArray(0, "avgRetryCntNode");
            for (j = 0; j < data[i].nNodes; j++) {
                if ((double)sim.node[j].nOutPkt > 0)
                    json.writeDouble("", (double)sim.node[j].totalRetryCnt / (double)sim.node[j].nOutPkt);
                else
                    json.writeDouble("", 0);
            }
            json.closeArray();
            */

            json.closeObject();

            sim.cleanUp();
        }
        json.closeArray();
        json.closeObject();
        f.close();

        //bkresult.close();
        std::time(&this->endTime);
        this->done.store(true);
    }

    void Worker::waitTillDone() {
        threadObj->join();
    }

    void Worker::writeToFile() {
        return;
    }
}

void runMM1k() {
    SimMM1k::DataPoint *data;
    int nNodes;
    double eps;
    std::string fname;
    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "eps: ";
    std::cin >> eps;
    std::cout << "outfile: ";
    std::cin >> fname;
    int k = 5;
    size_t nData = 200;
    data = nullptr;
    SimMM1k::generateInputDataWithFixedK(data, nNodes, k, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

void runMM1k_1() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower;
    double eps, ratio, rhoScale;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "k: ";
    std::cin >> k;
    std::cout << "ratioK0: ";
    std::cin >> ratio;
    std::cout << "rhoScale: ";
    std::cin >> rhoScale;
    std::cout << "eps: ";
    std::cin >> epsPower;
    eps = pow(10, -epsPower);
    std::cout << "nData: ";
    std::cin >> nData;
    //std::cout << "outfile: ";
    //std::cin >> fname;

    time_t startTime = std::time(0);
    fname = std::format("n{}_k{}_rk{}_rho{}_eps{}_d{}_{}.json", nNodes, k, ratio, rhoScale, epsPower, nData, startTime);
    
    data = nullptr;
    //SimMM1k::generateInputDataWithFixedK(data, nNodes, k, nData);
    SimMM1k::generateInputDataWithFixedKAndRatio(data, nNodes, k, ratio, rhoScale, nData);
    //SimMM1k::generateInputData(data, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

void runMM1k_2() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower;
    double eps, ratio, rhoScale, lambda, mu;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "lambda: ";
    std::cin >> lambda;
    std::cout << "mu: ";
    std::cin >> mu;
    std::cout << "eps: ";
    std::cin >> epsPower;
    eps = pow(10, -epsPower);
    std::cout << "nData: ";
    std::cin >> nData;
    //std::cout << "outfile: ";
    //std::cin >> fname;

    time_t startTime = std::time(0);
    fname = std::format("n{}_lambda{}_mu{}_eps{}_d{}_{}.json", nNodes, lambda, mu, epsPower, nData, startTime);

    data = nullptr;
    //SimMM1k::generateInputDataWithFixedK(data, nNodes, k, nData);
    //SimMM1k::generateInputDataWithFixedKAndRatio(data, nNodes, k, ratio, rhoScale, nData);
    SimMM1k::generateInputData2(data, nNodes, lambda, mu, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

void runMM1k_3() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower;
    double eps, ratio, rhoScale, lambda, mu;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "lambda: ";
    std::cin >> lambda;
    std::cout << "mu: ";
    std::cin >> mu;
    std::cout << "k: ";
    std::cin >> k;
    std::cout << "eps: ";
    std::cin >> epsPower;
    eps = pow(10, -epsPower);
    std::cout << "nData: ";
    std::cin >> nData;
    //std::cout << "outfile: ";
    //std::cin >> fname;

    time_t startTime = std::time(0);
    fname = std::format("n{}_lambda{}_mu{}_k{}_eps{}_d{}_{}.json", nNodes, lambda, mu, k, epsPower, nData, startTime);

    data = nullptr;
    SimMM1k::generateInputData3(data, nNodes, lambda, mu, k, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

void runMM1k_4() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower, kOrbit;
    double eps, ratio, rhoScale, lambda, mu;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "lambda: ";
    std::cin >> lambda;
    std::cout << "mu: ";
    std::cin >> mu;
    std::cout << "k: ";
    std::cin >> k;
    std::cout << "kOrbit: ";
    std::cin >> kOrbit;
    std::cout << "eps: ";
    std::cin >> epsPower;
    eps = pow(10, -epsPower);
    std::cout << "nData: ";
    std::cin >> nData;
    //std::cout << "outfile: ";
    //std::cin >> fname;

    time_t startTime = std::time(0);
    fname = std::format("n{}_lambda{}_mu{}_k{}_kOrb{}_eps{}_d{}_{}.json", nNodes, lambda, mu, k, kOrbit, epsPower, nData, startTime);

    data = nullptr;
    SimMM1k::generateInputData4(data, nNodes, lambda, mu, k, kOrbit, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

void runMM1k_5() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower, kOrbit;
    double eps, ratio, rhoScale, lambda, mu;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "lambda: ";
    std::cin >> lambda;
    std::cout << "mu: ";
    std::cin >> mu;
    std::cout << "eps: ";
    std::cin >> epsPower;
    eps = pow(10, -epsPower);
    std::cout << "nData: ";
    std::cin >> nData;
    //std::cout << "outfile: ";
    //std::cin >> fname;

    time_t startTime = std::time(0);
    fname = std::format("n{}_lambda{}_mu{}_noOrb_eps{}_d{}_{}.json", nNodes, lambda, mu, epsPower, nData, startTime);

    data = nullptr;
    SimMM1k::generateInputData5(data, nNodes, lambda, mu, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.run(eps, fname);
    worker.waitTillDone();
}

int main(int argc, char* argv[])
{
    //runMM1k();
    char key;
    
    std::cout << "1. Run simulation with random generated data\n";
    std::cout << "2. Run simulation with fixed lambda and varied buffer size\n";
    std::cout << "3. Run simulation with fixed lambda, fixed buffer size but varied orbit buffer size\n";
    std::cout << "4. Run simulation with fixed lambda, fixed buffer size but varied orbit mu\n";
    std::cout << "5. Run simulation with fixed lambda, orbit turned off and fixed buffer size\n";
    std::cout << "Choice: ";
    std::cin >> key;

    if (key == '1') {
        runMM1k_1();
    }
    else if (key == '2') {
        runMM1k_2();
    }
    else if (key == '3') {
        runMM1k_3();
    }
    else if (key == '4') {
        runMM1k_4();
    }
    else if (key == '5') {
        runMM1k_5();
    }

    return 0;

    //int nWorkers = 1;
    //int nData = 10000;
    //int nNodes = 1;
    //int nPackets = 100000;
    //time_t startTime, endTime;
    //std::string str;

    //int i;
    //Worker** worker;

    //Conf::initProgram();

    //if (Conf::workingMode.compare("batch") == 0) {
    //    BatchMode batchmode(argc, argv);
    //}
    //else {
    //    std::cout << "1. Run simulation (default option)\n";
    //    std::cout << "2. Batch mode\n";
    //    std::cout << "3. Sort data\n";
    //    std::cout << "4. M/M/1/k sim\n";
    //    std::cout << "Choose: ";
    //    std::cin >> str;
    //    if (str == "2") {
    //        BatchMode batchmode();
    //        return 0;
    //    }
    //    if (str == "3") {
    //        SortData::sortByDelta();
    //        return 0;
    //    }
    //    if (str == "4") {
    //        runMM1k();
    //        return 0;
    //    }

    //    std::cout << "List of PRNG:\n";
    //    std::cout << "1. C++ Standard Library Mersenne Twister\n";
    //    std::cout << "2. RDRAND instruction\n";
    //    std::cout << "Choose: ";
    //    std::cin >> str;
    //    if (str.compare("2") == 0) {
    //        Conf::rngType = Conf::RDRAND;
    //        std::cout << "PRNG is set to RDRAND\n";
    //    }
    //    else {
    //        std::cout << "PRNG is set to Mersenne Twister\n";
    //    }

    //    std::cout << "Number of workers: ";
    //    std::cin >> nWorkers;
    //    std::cout << "Number of simulations each worker: ";
    //    std::cin >> nData;
    //    std::cout << "Number of nodes each simulation: ";
    //    std::cin >> nNodes;
    //    std::cout << "Number of packets each simulation: ";
    //    std::cin >> nPackets;

    //    std::time(&startTime);
    //    worker = new Worker * [nWorkers];
    //    for (i = 0; i < nWorkers; i++) {
    //        worker[i] = new Worker(i, "", nData, nNodes, nPackets, startTime, Conf::gamma);
    //        worker[i]->run();
    //    }
    //    for (i = 0; i < nWorkers; i++) {
    //        worker[i]->waitTillDone();
    //    }
    //    std::time(&endTime);
    //    std::cout << "Running time: " << endTime - startTime << " seconds";

    //    /*for (i = 0; i < nWorkers; i++) {
    //        worker[i]->writeToFile();
    //    }*/
    //}

    //return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
