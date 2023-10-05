module;

#include <stdexcept>
#include <queue>
#include <algorithm>
#include <vector>
#include <list>
#include <iostream>

#include "Typedefs.h"

export module SimCore;
import RNG;

namespace SimCore {
	void assert(bool cond) {
		if (!cond)
			throw std::runtime_error("assertion failed");
	}

	/* ----------------------------------------------------------------------------- */
	/* Type definitions */

	class Sim;
	class Object;
	class Node;
	enum Message { arrive, done, gen };
	class Packet;

	class Event {
	public:
		double t;
		Object* receiver;
		Message action;
		Packet* pkt;

		Event();
		static void init(Event& ev, Object* rec, Message act, double t, Packet* p);

		bool operator<(const Event& Rhs) const {
			return this->t < Rhs.t;
		}

		bool operator>(const Event& Rhs) const {
			return this->t > Rhs.t;
		}

		bool operator==(const Event& Rhs) const {
			return this->t == Rhs.t;
		}
	};

	class EventList {
	private:
		std::vector<Event> heap;
		Event* events;
		size_t n, maxlen;
	public:
		EventList();
		void put(Event ev);
		bool get(Event& ev);
		bool notEmpty();
	};

	class Packet {
	public:
		bool dropped, inSystem, fromRetryQueue;
		size_t id;
		int genId, retryCnt;
		double creationTime, totalServTime, totalQueueTime;
		double arriveTime;

		Packet() {
			this->totalQueueTime = 0;
			this->totalServTime = 0;
			this->dropped = false;
			this->retryCnt = 0;
			this->inSystem = false;
			this->fromRetryQueue = false;
		};
	};
	
	class Object {
	protected:
		Sim* sim;
	public:
		virtual void handle(Message msg, double t, Packet* pkt) {}
	};

	class Node : public Object {

	};

	struct State {
		std::vector<double> tIntensity;
		std::vector<double> tFIntensity;
		std::vector<int> tDest;
		std::vector<bool> tType;
		double intensity;
	};

	export void resize_Matrix(Matrix& mat, size_t m, size_t n) {
		mat.resize(m);
		for (size_t i = 0; i < m; i++) {
			mat[i].resize(n);
		}
	}

	class Flow : public Object {
	private:
		int id;
		bool isMAPflow;

		RNG::Generator rng;
		double lambda;
		Matrix d0, d1;
		std::vector<State> states;
		int curState;

		double t;
		Node* outNode;

		double getInverval();
	public:
		Flow(Sim* sim, double lambda, int id);
		Flow(Sim* sim, double lambda, Matrix& d0, Matrix& d1, int id);
		void setOutput(Node* out);
		void generate();
		double getTime();
		double getLambda();
		const Matrix& getD(int i);
		void handle(Message msg, double t, Packet* pkt);
	};

	class SinkNode : public Node {
	public:
		void handle(Message msg, double t, Packet* pkt);
		SinkNode(Sim* sim);
	};

	class BasicNode : public Node {
	private:
		RNG::Generator rng;
		double mu;

		double t;
		double curLen, lastTimeLenChanged;

		std::queue<Packet*> queue;
		size_t queueLen, nInQueue, nDrop, maxRetry;
		bool busyState;
		Node* outNode;
		Node* retryQueue;

		double getServiceTime();
	public:
		double avgLen, totalPktTime;
		double outFlow;
		uint64_t nOutPkt, totalRetryCnt;

		BasicNode(Sim* sim, double mu, size_t queueLen);
		void setOutput(Node* out);
		void setRetryQueue(Node* rq);
		void handle(Message msg, double evTime, Packet* pkt);
		bool isFull();
		bool isRetryQueue;
	};

	export class Sim {
	private:
		void generate();
		void resetStateVars(size_t limit);
		void resetStateVars2();
	public:
		size_t nDrop;
		double avgTime;
		std::vector<Packet> pkt;
		std::list<Packet> tmpPkt;
		
		size_t dataCnt;
		size_t gid;
		std::vector<Flow> gen;
		std::vector<BasicNode> node;
		SinkNode* sink;
		int nGen, nNode;
		EventList evList;
		
		double totalTimeInSys, avgTimeInSys;
		double avgLen, lastTimeLenChanged, curLen;
		double nPktAtSink;
		double avgRetryCnt;
		size_t nPkt;

		Sim();
		void init3(int n, double lambda, double mu[], size_t queueLen[]);
		void init4(int n, double lambda, Matrix& d1, Matrix& d0, double mu[], size_t queueLen[]);
		void runSim(size_t limit);
		void runSim2(double eps);
		void runSim3(double eps);
		void cleanUp();
	};

	/* ----------------------------------------------------------------------------- */
	/* Functions of Event class */

	Event::Event() {
	}

	void Event::init(Event& ev, Object* rec, Message act, double t, Packet* p) {
		ev.receiver = rec;
		ev.action = act;
		ev.pkt = p;
		ev.t = t;
	}

	/* ----------------------------------------------------------------------------- */
	/* Functions of EventList class */

	EventList::EventList() {
		events = new Event[128];
		n = 0;
		maxlen = 128;
	}

	void EventList::put(Event ev) {
		heap.push_back(ev);
		std::push_heap(heap.begin(), heap.end(), std::greater<>{});
	}

	bool EventList::get(Event &ev) {
		bool success;
		if (heap.size() > 0) {
			ev = heap[0];
			std::pop_heap(heap.begin(), heap.end(), std::greater<>{});
			heap.pop_back();
			success = true;
		}
		else {
			success = false;
		}
		return success;
	}

	bool EventList::notEmpty() {
		return (!heap.empty());
	}

	/* ----------------------------------------------------------------------------- */
	/* Misc */

	void dropPacket(Sim* sim, Node* node, Packet* pkt, double t) {
		if (pkt->inSystem) {
			sim->avgLen = sim->avgLen + sim->curLen * (t - sim->lastTimeLenChanged);
			sim->lastTimeLenChanged = t;
			sim->curLen = sim->curLen - 1;
		}
		pkt->dropped = true;
		sim->nDrop++;
		delete pkt;
	}

	/* ----------------------------------------------------------------------------- */
	/* Functions of Flow class */

	double Flow::getInverval() {
		return this->rng.generateExponential(this->lambda);
	}

	Flow::Flow(Sim* sim, double lambda, int id) {
		this->id = id;
		this->lambda = lambda;
		this->sim = sim;
		this->t = 0;
		this->isMAPflow = false;
	}

	Flow::Flow(Sim* sim, double lambda, Matrix& d0, Matrix& d1, int id) {
		this->id = id;
		this->lambda = lambda;
		this->sim = sim;
		this->t = 0;
		this->isMAPflow = true;
		this->curState = 0;

		this->d0 = d0;
		this->d1 = d1;

		size_t i, j;
		double sum;
		/*
		states.resize(d0.size());
		for (i = 0; i < d0.size(); i++) {
			states[i].tIntensity.clear();
			states[i].tFIntensity.clear();
			states[i].tDest.clear();
			states[i].tType.clear();
			sum = 0;
			for (j = 0; j < d0[i].size(); j++) {
				if (i == j) {
					states[i].intensity = abs(d0[i][j]);
				}
				else if (d0[i][j] > 0) {
					states[i].tIntensity.push_back(d0[i][j]);
					sum = sum + d0[i][j];
					states[i].tFIntensity.push_back(sum);
					states[i].tDest.push_back(j);
					states[i].tType.push_back(false);
				}
				if (d1[i][j] > 0) {
					states[i].tIntensity.push_back(d1[i][j]);
					sum = sum + d1[i][j];
					states[i].tFIntensity.push_back(sum);
					states[i].tDest.push_back(j);
					states[i].tType.push_back(true);
				}
			}
		}
		*/
		for (i = 0; i < d0.size(); i++) {
			sum = 0;
			for (j = 0; j < d0[i].size(); j++) {
				sum = sum + d1[i][j];
				if (i != j) {
					sum = sum + d0[i][j];
				}
			}
			d0[i][i] = abs(d0[i][i] - sum);
		}
	}

	void Flow::setOutput(Node* out) {
		this->outNode = out;
	}

	void Flow::generate() {
		size_t pktIdx, i, iMin;
		double interval = 0;
		Packet* pkt;
		Event ev, ev2;
		bool newPacket = false;

		if (!this->isMAPflow && this->lambda > 0.0) {
			interval = this->getInverval();
			assert(interval >= 0);
			this->t = this->t + interval;
			newPacket = true;
		}
		else if (this->isMAPflow) {
			std::vector<double> intervals;
			intervals.resize(d0.size() * 2);
			for (i = 0; i < d0.size(); i++) {
				if (d0[curState][i] > 0) intervals[i * 2] = rng.generateExponential(d0[curState][i]);
				else intervals[i * 2] = std::numeric_limits<double>::infinity();
				if (d1[curState][i] > 0) intervals[i * 2 + 1] = rng.generateExponential(d1[curState][i]);
				else intervals[i * 2 + 1] = std::numeric_limits<double>::infinity();
			}
			iMin = 0;
			for (i = 1; i < intervals.size(); i++) {
				if (intervals[i] < intervals[iMin]) iMin = i;
			}
			interval = intervals[iMin];
			if (interval < std::numeric_limits<double>::infinity()) {
				assert(interval >= 0);
				this->t = this->t + interval;
				curState = iMin / 2;
				newPacket = (iMin % 2 != 0);
			}
			else {
				interval = 0;
			}
		}

		if (newPacket) {
			pkt = new Packet();
			pkt->id = sim->gid;
			sim->gid++;
			sim->nPkt++;
			pkt->genId = this->id;
			pkt->creationTime = t;

			Event::init(ev, this->outNode, arrive, t, pkt);
			sim->evList.put(ev);
		}

		if (interval > 0) {
			Event::init(ev2, this, gen, t, nullptr);
			sim->evList.put(ev2);
		}
	}

	double SimCore::Flow::getTime() {
		return this->t;
	}

	double SimCore::Flow::getLambda() {
		return this->lambda;
	}

	const Matrix& SimCore::Flow::getD(int i) {
		if (i==0) return this->d0;
		return this->d1;
	}

	void SimCore::Flow::handle(Message msg, double t, Packet* pkt) {
		assert(msg == gen);
		generate();
	}

	/* ----------------------------------------------------------------------------- */
	/* Functions of SinkNode class */

	SinkNode::SinkNode(Sim* sim) {
		this->sim = sim;
	}

	void SinkNode::handle(Message msg, double t, Packet* pkt) {
		assert(msg == arrive);

		sim->avgLen = sim->avgLen + sim->curLen * (t - sim->lastTimeLenChanged);
		sim->lastTimeLenChanged = t;
		sim->curLen = sim->curLen - 1;
		
		sim->nPktAtSink = sim->nPktAtSink + 1.0;
		double duration = t - pkt->creationTime;
		sim->totalTimeInSys = sim->totalTimeInSys + duration;
		sim->avgRetryCnt += pkt->retryCnt;
		delete pkt;
		//printf("Packet %d: created at %.4f arrived to sink at %.4f\n", pkt->id, pkt->creationTime, t);
	}

	/* ----------------------------------------------------------------------------- */
	/* Functions of BasicNode class */

	double BasicNode::getServiceTime() {
		return this->rng.generateExponential(this->mu);
	}

	BasicNode::BasicNode(Sim* sim, double mu, size_t queueLen) {
		this->mu = mu;
		this->sim = sim;
		this->t = 0;
		this->queueLen = queueLen;
		this->nInQueue = 0;
		this->busyState = false;
		this->nDrop = 0;
		this->maxRetry = 0;

		this->avgLen = 0;
		this->curLen = 0;
		this->lastTimeLenChanged = 0;

		this->nOutPkt = 0;
		this->outFlow = 0;
		this->totalPktTime = 0;
		this->totalRetryCnt = 0;

		this->outNode = nullptr;
		this->retryQueue = nullptr;
		this->isRetryQueue = false;
	}

	void BasicNode::setOutput(Node* out) {
		this->outNode = out;
	}

	void BasicNode::setRetryQueue(Node* rq) {
		this->retryQueue = rq;
	}

	void BasicNode::handle(Message msg, double evTime, Packet* pkt) {
		if (msg == arrive) {
			if (!this->busyState || this->nInQueue < this->queueLen) {
				if (!pkt->inSystem) {
					pkt->inSystem = true;
					this->sim->avgLen = this->sim->avgLen + this->sim->curLen * (evTime - this->sim->lastTimeLenChanged);
					this->sim->lastTimeLenChanged = evTime;
					this->sim->curLen = this->sim->curLen + 1;
				}

				this->avgLen = this->avgLen + this->curLen * (evTime - this->lastTimeLenChanged);
				this->lastTimeLenChanged = evTime;
				this->curLen = this->curLen + 1;

				pkt->arriveTime = evTime;
			}
			if (!this->busyState) {
				/* Available */
				if (this->mu > 0) {
					double serviceTime = this->getServiceTime();
					this->t = evTime + serviceTime;
					pkt->totalServTime = pkt->totalServTime + serviceTime;

					Event ev;
					Event::init(ev, this, done, this->t, pkt);
					this->sim->evList.put(ev);
				}
				this->busyState = true;
			}
			else if (this->nInQueue < this->queueLen) {
				/* Put into queue */
				this->queue.push(pkt);
				this->nInQueue++;
				pkt->totalQueueTime = pkt->totalQueueTime - evTime;
			}
			else {
				/* Queue is full */
				if (this->retryQueue == nullptr) {
					/* Drop packet */
					dropPacket(sim, this, pkt, evTime);
					this->nDrop++;
				}
				else {
					/* Try sending it to another node */
					/*
					pkt->retryCnt++;
					if (pkt->retryCnt <= this->maxRetry) {
						Event ev;
						Event::init(ev, this->retryQueue, arrive, evTime, pkt);
						sim->evList.put(ev);
					}
					else {
						dropPacket(sim, this, pkt, evTime);
						this->nDrop++;
					}
					*/
					this->nDrop++;
					if (!pkt->fromRetryQueue) {
						if (pkt->inSystem) pkt->retryCnt++;

						Event ev;
						Event::init(ev, this->retryQueue, arrive, evTime, pkt);
						sim->evList.put(ev);
					}
					else {
						dropPacket(sim, this, pkt, evTime);
					}
				}
			}
		}
		else if (msg == done) {
			assert(this->t == evTime);
			//if (!this->isRetryQueue || !(((BasicNode*)(this->outNode))->isFull())) {
			if (true) {
				pkt->fromRetryQueue = this->isRetryQueue;

				this->avgLen = this->avgLen + this->curLen * (evTime - this->lastTimeLenChanged);
				this->lastTimeLenChanged = evTime;
				this->curLen = this->curLen - 1;

				this->nOutPkt++;
				this->outFlow = (double)nOutPkt / evTime;
				this->totalPktTime += (evTime - pkt->arriveTime);
				this->totalRetryCnt += pkt->retryCnt;

				Event ev;
				Event::init(ev, this->outNode, arrive, evTime, pkt);
				sim->evList.put(ev);

				/* if queue is not empty */
				if (this->nInQueue > 0) {
					Packet* pkt2 = this->queue.front();
					this->queue.pop();
					this->nInQueue--;
					pkt2->totalQueueTime = pkt2->totalQueueTime + evTime;

					double serviceTime = this->getServiceTime();
					this->t = this->t + serviceTime;
					pkt2->totalServTime = pkt2->totalServTime + serviceTime;

					Event ev2;
					Event::init(ev2, this, done, this->t, pkt2);
					sim->evList.put(ev2);
				}
				else {
					/* node is available now */
					this->busyState = false;
				}
			}
			else {
				double serviceTime = this->getServiceTime();
				this->t = evTime + serviceTime;
				pkt->totalServTime = pkt->totalServTime + serviceTime;

				Event ev;
				Event::init(ev, this, done, this->t, pkt);
				sim->evList.put(ev);
			}
		}
		else
			assert(false);
	}

	bool BasicNode::isFull() {
		return (this->nInQueue == this->queueLen);
	}

	/* ----------------------------------------------------------------------------- */
	/* Functions of Sim class */

	Sim::Sim() {
		sink = nullptr;
	}

	void Sim::generate() {
		int i;
		double min = 0, max = 0;
		double t;

		for (i = 0; i < nGen; i++) {
			gen[i].generate();
		}
	}

	void Sim::resetStateVars(size_t limit)
	{
		nDrop = 0;
		avgTime = 0;
		dataCnt = limit;

		pkt.resize(limit);

		tmpPkt.clear();

		gid = 0;
	}

	void Sim::runSim(size_t limit) {
		Event ev;

		resetStateVars(limit);
		generate();
		do {
			evList.get(ev);
			ev.receiver->handle(ev.action, ev.t, ev.pkt);
		} while (gid < limit);

		/* process the remaining events */
		while (evList.notEmpty()) {
			evList.get(ev);
			//ev.receiver->handle(ev.action, ev.t, ev.pkt);
		}
	}

	void Sim::cleanUp() {
		gen.clear();
		node.clear();
		delete sink;
		sink = nullptr;
		pkt.clear(); pkt.shrink_to_fit();
		tmpPkt.clear();
	}

	void Sim::resetStateVars2()
	{
		nDrop = 0;
		avgTime = 0;

		totalTimeInSys = 0;
		nPktAtSink = 0;
		nPkt = 0;
		
		pkt.clear();
		tmpPkt.clear();

		gid = 0;

		avgLen = 0;
		lastTimeLenChanged = 0;
		curLen = 0;
		avgRetryCnt = 0;
	}

	void Sim::runSim2(double epsTarget) {
		#define WINDOW_SIZE 10000
		Event ev;
		double W[WINDOW_SIZE];
		double len[WINDOW_SIZE];
		double delta[WINDOW_SIZE];
		double eps = 1;
		int cnt = 0;
		int stage = 0;
		double nPktAtSink0, n;
		double totalTimeInSys0 = 0;
		size_t evCnt = 0;

		resetStateVars2();
		n = 0; nPktAtSink0 = 0;
		generate();
		do {
			evCnt++;
			/*if ((evCnt - 1) % 1000000 == 0) {
				std::cout << "[";
				for (int i = 0; i < this->nNode; i++) {
					std::cout << this->node[i].nOutPkt << "; ";
				}
				std::cout << "]\n";
			}*/
			evList.get(ev);
			ev.receiver->handle(ev.action, ev.t, ev.pkt);
			if (nPktAtSink > n) {
				/*if (stage == 0) {
					len[cnt] = avgLen / lastTimeLenChanged;
					if (cnt > 0) {
						delta[cnt] = abs(len[cnt] - len[0]);
					}
					else {
						delta[cnt] = 0;
					}
					cnt++;
					if (cnt == WINDOW_SIZE) {
						int i;
						double sumLoss, sumDelta;
						sumLoss = 0;
						sumDelta = 0;
						for (i = 0; i < WINDOW_SIZE; i++) {
							sumLoss = sumLoss + len[i];
							sumDelta = sumDelta + delta[i];
						}
						eps = sumDelta / sumLoss;
						cnt = 0;
						if (eps <= epsTarget) {
							stage = 1;
							totalTimeInSys0 = totalTimeInSys;
							nPktAtSink0 = nPktAtSink;
							std::cout << "stage 0 done\n";
						}
					}
				}
				else {*/
				W[cnt] = (totalTimeInSys - totalTimeInSys0) / (nPktAtSink - nPktAtSink0);
				if (cnt > 0) {
					delta[cnt] = abs(W[cnt] - W[0]);
				}
				else {
					delta[cnt] = 0;
				}
				cnt++;
				if (cnt == WINDOW_SIZE) {
					int i;
					double sumW, sumDelta;
					sumW = 0;
					sumDelta = 0;
					for (i = 0; i < WINDOW_SIZE; i++) {
						sumW = sumW + W[i];
						sumDelta = sumDelta + delta[i];
					}
					eps = sumDelta / sumW;
					cnt = 0;
					//std::cout << W[0] << "\n";
				}
				//std::cout << "Eps: " << eps << "\n";
				n = nPktAtSink;
				// std::cout << W[0] << "\n";
			}
		} while (nPktAtSink - nPktAtSink0 <= WINDOW_SIZE || eps > epsTarget);

		avgTimeInSys = W[0];

		/* clear the remaining events */
		while (evList.notEmpty()) {
			evList.get(ev);
		}
	}

	void Sim::runSim3(double epsTarget) {
		Event ev;
		double W[WINDOW_SIZE];
		double len[WINDOW_SIZE];
		double delta[WINDOW_SIZE];
		double eps = 1;
		double epsLambda = 1;
		int cnt = 0;
		double nPktAtSink0, n, t, lambda;
		double totalTimeInSys0 = 0;
		size_t evCnt = 0;

		resetStateVars2();
		n = 0; nPktAtSink0 = 0;
		generate();
		do {
			evCnt++;
			evList.get(ev);
			t = ev.t;
			ev.receiver->handle(ev.action, ev.t, ev.pkt);
			if (nPktAtSink > n) {
				W[cnt] = (totalTimeInSys - totalTimeInSys0) / (nPktAtSink - nPktAtSink0);
				if (cnt > 0) {
					delta[cnt] = abs(W[cnt] - W[0]);
				}
				else {
					delta[cnt] = 0;
				}
				cnt++;
				if (cnt == WINDOW_SIZE) {
					int i;
					double sumW, sumDelta;
					sumW = 0;
					sumDelta = 0;
					for (i = 0; i < WINDOW_SIZE; i++) {
						sumW = sumW + W[i];
						sumDelta = sumDelta + delta[i];
					}
					eps = sumDelta / sumW;
					lambda = gen[0].getLambda();
					epsLambda = abs(lambda - (double)nPkt / t) / lambda;
					cnt = 0;
					//std::cout << W[0] << "\n";
					//std::cout << lambda << " " << (double)nPkt / t << "\n";
				}
				//std::cout << "Eps: " << eps << "\n";
				n = nPktAtSink;
			}
		} while (nPktAtSink - nPktAtSink0 <= WINDOW_SIZE || eps > epsTarget || epsLambda > epsTarget);

		avgTimeInSys = W[0];

		/* clear the remaining events */
		while (evList.notEmpty()) {
			evList.get(ev);
		}
	}

	void Sim::init3(int n, double lambda, double mu[], size_t queueLen[]) {
		int i;

		nGen = 1; nNode = n;
		gen.reserve(1); node.reserve(n);
		gen.emplace_back(this, lambda, 0);
		for (i = 0; i < n; i++) {
			node.emplace_back(this, mu[i], queueLen[i]);
		}
		sink = new SinkNode(this);
		gen[0].setOutput(&node[1]);
		for (i = n - 1; i >= 1; i--) {
			if (i == n - 1) node[i].setOutput(sink);
			else node[i].setOutput(&node[i + 1]);
			node[i].setRetryQueue(&node[0]);
		}
		node[0].setOutput(&node[1]);
		node[0].isRetryQueue = true;
	}

	void Sim::init4(int n, double lambda, Matrix& d0, Matrix& d1, double mu[], size_t queueLen[]) {
		int i;

		nGen = 1; nNode = n;
		gen.reserve(1); node.reserve(n);
		gen.emplace_back(this, lambda, d0, d1, 0);
		for (i = 0; i < n; i++) {
			node.emplace_back(this, mu[i], queueLen[i]);
		}
		sink = new SinkNode(this);
		gen[0].setOutput(&node[1]);
		for (i = n - 1; i >= 1; i--) {
			if (i == n - 1) node[i].setOutput(sink);
			else node[i].setOutput(&node[i + 1]);
			node[i].setRetryQueue(&node[0]);
		}
		node[0].setOutput(&node[1]);
		node[0].isRetryQueue = true;
	}
}