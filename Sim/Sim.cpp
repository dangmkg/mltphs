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

#include <Eigen/Dense>

#include "Typedefs.h"

import Conf;
import SimCore;
import RNG;
import Json;
import SortData;

unsigned __int64 generateUniform(unsigned __int64 limit) {
    unsigned __int64 result;
    _rdrand64_step(&result);
    return (result % limit);
}

namespace SimMM1k {
    class DataPoint {
    public:
        double lambda[128];
        double mu[128];
        size_t k[128];
        Matrix d0[128];
        Matrix d1[128];

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
        void mainFuncMAP();
        double eps;
    public:
        DataPoint* data;
        size_t nData;
        int id;

        Worker(int id, DataPoint* data, size_t nData);
        int getCurrentProgress();
        void run(double eps, std::string fname);
        void runMAP(double eps, std::string fname);
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

    double generateMatrix(Matrix& d0, Matrix& d1, size_t dim, double limit, double scale) {
        size_t i, j;
        double sum;
        Eigen::MatrixXd d0_(dim, dim), d1_(dim, dim);
        SimCore::resize_Matrix(d0, dim, dim);
        SimCore::resize_Matrix(d1, dim, dim);
        for (i = 0; i < d0.size(); i++) {
            for (j = 0; j < d0[i].size(); j++) {
                if (i != j) {
                    d0[i][j] = ((double)generateUniform(limit) + 1) / double(limit) * scale;
                    d1[i][j] = 0;
                }
                if (i == j) {
                    d0[i][j] = 0;
                    d1[i][j] = ((double)generateUniform(limit) + 1) / double(limit) * scale;
                }
            }
            sum = 0;
            for (j = 0; j < d0[i].size(); j++) {
                d0_(i, j) = d0[i][j];
                d1_(i, j) = d1[i][j];
                sum = sum + d0[i][j] + d1[i][j];
            }
            d0[i][i] = -sum;
            d0_(i, i) = -sum;
        }
        Eigen::MatrixXd a(dim, dim);
        a = (d0_ + d1_).transpose();
        for (i = 0; i < dim; i++) a(dim - 1, i) = 1.0;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(dim);
        b(dim - 1) = 1.0;
        Eigen::VectorXd x = a.colPivHouseholderQr().solve(b);
        for (i = 0; i < dim; i++) x(i) = x(i) * d1_(i,i);
        double lambda = x.sum();
        return lambda;
    }

    void generateMMPPWithFixedKAndRatio(DataPoint*& data, int nNodes, int k, double dim, double rhoScale, size_t nData) {
        size_t i, j, limit;
        double scale, lambda, rho;

        data = new DataPoint[nData];
        limit = 65536;
        scale = 100.0;

        for (i = 0; i < nData; i++) {
            data[i].lambda[0] = generateMatrix(data[i].d0[0], data[i].d1[0], dim, limit, scale);
            for (j = 0; j < nNodes; j++) {
                rho = (double)generateUniform(limit) - (limit / 2);
                rho = rho / limit * rhoScale * 2;
                rho = exp(rho);

                data[i].mu[j] = data[i].lambda[0] * rho;
                if (j != 0) data[i].lambda[j] = 0.0;
                data[i].k[j] = k;
            }
            //data[i].k[0] = round(k * ratioK0);
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

    void Worker::runMAP(double eps, std::string fname) {
        this->eps = eps;
        this->filename = fname;
        this->done.store(false);
        threadObj = new std::thread(&Worker::mainFuncMAP, this);
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
            std::chrono::duration<double> diff = end - start;

            std::cout << "Datapoint " << i << " avg outflow: " << sim.nPktAtSink / sim.lastTimeLenChanged << " simulation duration: " << diff.count() << "\n";
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

    void Worker::mainFuncMAP() {
        size_t i, j;
        std::string msg, bkfn;
        Json::Writer json;
        SimCore::Sim sim;
        std::ofstream f;

        f.open(this->filename);
        json.setOutput(&f);
        json.openObject(1, "");
        json.openArray(1, "data");
        for (i = 0; i < nData; i++) {
            std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
            sim.init4(data[i].nNodes, data[i].lambda[0], data[i].d0[0], data[i].d1[0], data[i].mu, data[i].k);
            sim.runSim3(this->eps);
            std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;

            std::cout << "Datapoint " << i << " avg outflow: " << sim.nPktAtSink / sim.lastTimeLenChanged << " simulation duration: " << diff.count() << "\n";
            json.openObject(1, "");
            json.writeInt("nNodes", data[i].nNodes);
            json.openArray(0, "k");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", data[i].k[j]);
            }
            json.closeArray();
            json.writeDouble("lambd", data[i].lambda[0]);
            json.openArray(0, "d0");
            for (j = 0; j < data[i].d0[0].size(); j++) {
                json.writeArrayOfDouble("", data[i].d0[0][j].data(), data[i].d0[0][j].size());
            }
            json.closeArray();
            json.openArray(0, "d1");
            for (j = 0; j < data[i].d1[0].size(); j++) {
                json.writeArrayOfDouble("", data[i].d1[0][j].data(), data[i].d1[0][j].size());
            }
            json.closeArray();
            json.openArray(0, "mu");
            for (j = 0; j < data[i].nNodes; j++) {
                json.writeDouble("", data[i].mu[j]);
            }
            json.closeArray();
            json.writeDouble("avgTime", sim.avgTimeInSys);
            json.writeDouble("avgLen", sim.avgLen / sim.lastTimeLenChanged);
            json.writeDouble("dropRate", (double)sim.nDrop / (double)sim.nPkt);
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

void runMM1k_6() {
    SimMM1k::DataPoint* data;
    int nNodes, k, epsPower;
    double eps, dim, rhoScale;
    std::string fname;
    size_t nData = 1000;

    std::cout << "nNodes: ";
    std::cin >> nNodes;
    std::cout << "k: ";
    std::cin >> k;
    std::cout << "dim of MAP Matrix: ";
    std::cin >> dim;
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
    fname = std::format("n{}_k{}_dim{}_rho{}_eps{}_d{}_{}.json", nNodes, k, dim, rhoScale, epsPower, nData, startTime);

    data = nullptr;
    SimMM1k::generateMMPPWithFixedKAndRatio(data, nNodes, k, dim, rhoScale, nData);
    SimMM1k::Worker worker(0, data, nData);
    worker.runMAP(eps, fname);
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
    std::cout << "6. Run simulation with MMPP(2) flow and random generated data\n";
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
    else if (key == '6') {
        runMM1k_6();
    }

    return 0;
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
