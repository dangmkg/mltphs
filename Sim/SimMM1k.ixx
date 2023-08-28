module;

export module SimMM1k;

#include <intrin.h>
#include <string>
#include <thread>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>

import Conf;
import SimCore;
import RNG;
import Json;

unsigned __int64 generateUniform(unsigned __int64 limit) {
    unsigned __int64 result;
    _rdrand64_step(&result);
    return (result % limit);
}

export namespace SimMM1k {
    class DataPoint {
    public:
        double lambda[16];
        double mu[16];
        size_t k[16];

        int nNodes;
        double meanNumOfPktsAtNode[16];
        double meanPktTimeAtNode[16];
        size_t totalNumOfPktsAtNode[16];
        size_t numOfDroppedPktsAtNode[16];
        unsigned __int64 seedsLambda[16];
        unsigned __int64 seedsMu[16];
    };

    class Worker {
        std::string filename;
        time_t startTime, endTime;
        int nPackets, curPoint;
        std::thread* threadObj;
        std::atomic_bool done;
        void mainFunc();
    public:
        DataPoint *data;
        size_t nData;
        int id;

        Worker(int id, DataPoint *data, size_t nData);
        int getCurrentProgress();
        void run();
        void writeToFile();
        void waitTillDone();
        bool isFree();
    };

    void generateInputDataWithFixedK(DataPoint *data, int nNodes, int k, size_t nData) {
        size_t i, j, limit;
        double scale, lambda, rho;

        data = new DataPoint[nData];
        limit = 65536;
        scale = 10.0;

        for (i = 0; i < nData; i++) {
            lambda = (double)generateUniform(limit) + 1;
            for (j = 0; j < nNodes; j++) {
                rho = (double)generateUniform(limit) - (limit / 2);
                rho = rho / limit * 9.0;
                rho = exp(rho);

                data[i].mu[j] = lambda * rho / double(limit) * scale;
                data[i].lambda[j] = 0.0;
                data[i].k[j] = k;
            }
            data[i].lambda[0] = lambda / double(limit) * scale;
            data[i].nNodes = nNodes;
        }
    }

    Worker::Worker(int id, DataPoint *data, size_t nData) {
        this->id = id;
        this->done.store(true);
        this->data = data;
        this->nData = nData;
    }

    int Worker::getCurrentProgress() {
        return curPoint;
    }

    void Worker::run() {
        this->done.store(false);
        threadObj = new std::thread(&Worker::mainFunc, this);
    }

    void Worker::mainFunc() {
        size_t i, j;
        std::string msg, bkfn;
        Json::Writer json;
        time_t start, end;
        SimCore::Sim sim;

        for (i = 0; i < nData; i++) {
            std::time(&start);
            sim.init(data[i].nNodes, data[i].lambda, data[i].mu, data[i].k);
            sim.runSim2(1e-2);
            std::time(&end);

            std::cout << data[i].nNodes << ", ";
            std::cout << data[i].k << ", ";
            std::cout << data[i].lambda[0] << ", ";
            for (j = 0; j < data[i].nNodes; j++) {
                std::cout << data[i].mu[j] << ", ";
            }
            std::cout << sim.avgTimeInSys << ", ";
            std::cout << sim.nDrop / sim.nPkt << ", ";
            std::cout << end - start << "\n";

            sim.cleanUp();
        }

        //bkresult.close();
        std::time(&this->endTime);
        this->done.store(true);
    }

    void Worker::waitTillDone() {
        threadObj->join();
    }

    void Worker::writeToFile() {
        /*std::string key;
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
        out.close();*/
    }
}
