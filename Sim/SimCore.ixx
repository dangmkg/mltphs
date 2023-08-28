module;

#include <stdexcept>
#include <queue>
#include <algorithm>
#include <vector>
#include <list>
#include <iostream>

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
	class Node;
	enum Message { arrive, done };
	class Packet;

	class Event {
	public:
		double t;
		Node* receiver;
		Message action;
		Packet* pkt;

		Event();
		static void init(Event& ev, Node* rec, Message act, double t, Packet* p);

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
	};
	
	class Node {
	protected:
		Sim* sim;
	public:
		virtual void handle(Message msg, double t, Packet* pkt) {}
	};

	class GenNode : public Node {
	private:
		int id;

		RNG::Generator rng;
		double lambda;

		double t;
		Node* outNode;

		double getInverval();
	public:
		GenNode(Sim* sim, double lambda, int id);
		void setOutput(Node* out);
		void generate();
		double getTime();
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
		double genMinT, genMaxT;

		void generateUntil(double threshold);
		void resetStateVars(size_t limit);
		void resetStateVars2();
	public:
		size_t nDrop;
		double avgTime;
		std::vector<Packet> pkt;
		std::list<Packet> tmpPkt;
		
		size_t dataCnt;
		size_t gid;
		std::vector<GenNode> gen;
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
		void init(int n, double lambda[], double mu[], size_t queueLen[]);
		void init2(int n, double lambda[], double mu[], size_t queueLen[]);
		void init3(int n, double lambda, double mu[], size_t queueLen[]);
		void runSim(size_t limit);
		void runSim2(double eps);
		void cleanUp();
	};

	/* ----------------------------------------------------------------------------- */
	/* Functions of Event class */

	Event::Event() {
	}

	void Event::init(Event& ev, Node* rec, Message act, double t, Packet* p) {
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
	/* Functions of GenNode class */

	double GenNode::getInverval() {
		return this->rng.generateExponential(this->lambda);
	}

	GenNode::GenNode(Sim* sim,  double lambda, int id) {
		this->id = id;
		this->lambda = lambda;
		this->sim = sim;
		this->t = 0;
	}

	void GenNode::setOutput(Node* out) {
		this->outNode = out;
	}

	void GenNode::generate() {
		size_t pktIdx;
		double interval;
		Packet* pkt;
		Event ev;

		if (this->lambda > 0.0) {
			interval = this->getInverval();
			assert(interval >= 0);
			this->t = this->t + interval;
			/*
			pktIdx = sim->gid;
			if (pktIdx < sim->dataCnt)
				pkt = &sim->pkt[pktIdx];
			else {
				pkt = &sim->tmpPkt.emplace_back();
			}
			*/
			pkt = new Packet();
			pkt->id = sim->gid;
			sim->gid++;
			sim->nPkt++;
			pkt->genId = this->id;
			pkt->creationTime = t;
			pkt->totalQueueTime = 0;
			pkt->totalServTime = 0;
			pkt->dropped = false;
			pkt->retryCnt = 0;
			pkt->inSystem = false;
			pkt->fromRetryQueue = false;

			Event::init(ev, this->outNode, arrive, t, pkt);
			sim->evList.put(ev);
		}
	}

	double SimCore::GenNode::getTime() {
		return this->t;
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

	/* POSTCONDITION: this procedure ensures all possible packets
	   before and at threshold time have been generated */
	void Sim::generateUntil(double threshold) {
		int i;
		double min = 0, max = 0;
		double t;

		for (i = 0; i < nGen; i++) {
			while (gen[i].getTime() <= threshold) {
				gen[i].generate();
			}
		}
		
		min = gen[0].getTime(); max = min;
		for (i = 1; i < nGen; i++) {
			t = gen[i].getTime();
			if (t < min) min = t;
			else if (t > max) max = t;
		}
		assert(min > genMinT); genMinT = min;
		assert(max >= genMaxT); genMaxT = max;
	}

	void Sim::resetStateVars(size_t limit)
	{
		nDrop = 0;
		avgTime = 0;
		dataCnt = limit;

		pkt.resize(limit);

		tmpPkt.clear();

		gid = 0;
		genMinT = 0;
		genMaxT = 0;
	}

	/* Init a linear network of Poisson flow generators and nodes with exponential service time */
	void Sim::init(int n, double lambda[], double mu[], size_t queueLen[]) {
		int i;

		nGen = n; nNode = n;
		gen.reserve(n); node.reserve(n);
		for (i = 0; i < n; i++) {
			gen.emplace_back(this, lambda[i], i);
			node.emplace_back(this, mu[i], queueLen[i]);
		}
		sink = new SinkNode(this);
		for (i = n-1; i >= 0; i--) {
			gen[i].setOutput(&node[i]);
			if (i == n - 1) node[i].setOutput(sink);
			else node[i].setOutput(&node[i + 1]);
		}
	}

	void Sim::runSim(size_t limit) {
		Event ev;

		resetStateVars(limit);
		generateUntil(0);
		do {
			evList.get(ev);
			if (ev.t <= genMinT)
				ev.receiver->handle(ev.action, ev.t, ev.pkt);
			else {
				evList.put(ev);
				double threshold = ev.t;
				generateUntil(threshold);	
			}
		} while (gid < limit);

		/* process the remaining events */
		while (evList.notEmpty()) {
			evList.get(ev);
			ev.receiver->handle(ev.action, ev.t, ev.pkt);
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

	/* M/M/1/k */
	void Sim::init2(int n, double lambda[], double mu[], size_t queueLen[]) {
		int i;

		nGen = 1; nNode = n;
		gen.reserve(1); node.reserve(n);
		gen.emplace_back(this, lambda[0], 0);
		for (i = 0; i < n; i++) {
			node.emplace_back(this, mu[i], queueLen[i]);
		}
		sink = new SinkNode(this);
		gen[0].setOutput(&node[0]);
		for (i = n - 1; i >= 0; i--) {
			if (i == n - 1) node[i].setOutput(sink);
			else node[i].setOutput(&node[i + 1]);
		}
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
		genMinT = 0;
		genMaxT = 0;

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
		generateUntil(0);
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
			if (ev.t <= genMinT)
				ev.receiver->handle(ev.action, ev.t, ev.pkt);
			else {
				evList.put(ev);
				double threshold = ev.t;
				generateUntil(threshold);
			}
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
					std::cout << W[0] << "\n";
				}
				//std::cout << "Eps: " << eps << "\n";
				n = nPktAtSink;
				//std::cout << W[0] << "\n";
			}
		} while (nPktAtSink - nPktAtSink0 <= WINDOW_SIZE || eps > epsTarget);

		avgTimeInSys = W[0];

		/* process the remaining events */
		while (evList.notEmpty()) {
			evList.get(ev);
			ev.receiver->handle(ev.action, ev.t, ev.pkt);
		}
		//avgTimeInSys = totalTimeInSys / nPktAtSink;
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

	class StateSimEvent {
	public:
		int sender, recv;
		double time;
		int msg;
		#define MSG_NEW 1
		#define MSG_DEL 0

		bool operator<(const StateSimEvent& Rhs) const {
			return this->time < Rhs.time;
		}

		bool operator>(const StateSimEvent& Rhs) const {
			return this->time > Rhs.time;
		}

		bool operator==(const StateSimEvent& Rhs) const {
			return this->time == Rhs.time;
		}
	};

	export class StateSimMM1k {
	private:
		size_t stepCnt;
		size_t pktId;
		std::vector<StateSimEvent> evQueue;
		std::vector<size_t> state;
		double t, time;
		RNG::Generator rand;

		size_t genIdx(int senderId) {
			assert(senderId < 0);
			return -senderId - 1;
		};
		int senderId(size_t idx, bool isGenIdx) {
			if (isGenIdx) {
				assert(idx < -INT_MIN);
				return -(int)idx - 1;
			}
			else {
				return idx;
			}
		};

		void gen(size_t idx);
		void serv(size_t idx);
		void doStep();
	public:
		int nNodes;
		std::vector<double> lambda, mu;
		std::vector<size_t> k;
		std::vector<int> outNodeL, outNodeM;
		std::vector<double> Pi;
		std::vector<size_t> nPkt, nPktDrop;

		size_t stateIdx(const std::vector<size_t>& state, const std::vector<size_t>& k) {
			size_t idx = state[0];
			size_t arrSize = 1;
			for (size_t i = 1; i < state.size(); i++) {
				arrSize = arrSize * k[i - 1];
				idx = idx + state[i] * arrSize;
			}
			return idx;
		};

		void init(std::vector<double>& lambda0, std::vector<double>& mu0, std::vector<size_t>& k0);
		void runSteps(size_t maxStep);
		void runEps(double eps, std::vector<double> L0);
	};

	void StateSimMM1k::init(std::vector<double>& lambda0, std::vector<double>& mu0, std::vector<size_t>& k0) {
		lambda = lambda0;
		mu = mu0;
		k = k0;

		assert(lambda.size() == mu.size());
		assert(lambda.size() == k.size());

		size_t i;
		size_t nStates = 1;
		for (i = 0; i < k0.size(); i++) nStates = nStates * k0[i];
		Pi.resize(nStates);
		nPkt.resize(mu.size());
		nPktDrop.resize(mu.size());
		state.resize(mu.size());

		outNodeL.resize(lambda.size());
		for (i = 0; i < lambda.size(); i++) {
			outNodeL[i] = i;
		}
		outNodeM.resize(mu.size());
		for (i = 0; i < mu.size(); i++) {
			outNodeM[i] = i + 1;
		}
		outNodeM[mu.size() - 1] = -1; // sink
	}

	void StateSimMM1k::gen(size_t idx) {
		double interval, evTime;
		if (this->lambda[idx] > 0) {
			interval = rand.generateExponential(this->lambda[idx]);
			StateSimEvent ev;
			ev.msg = MSG_NEW;
			ev.time = this->time + interval;
			ev.sender = senderId(idx, true);
			ev.recv = senderId(this->outNodeL[idx], false);

			evQueue.push_back(ev);
			std::push_heap(evQueue.begin(), evQueue.end(), std::greater<>{});
		}
	}

	void StateSimMM1k::serv(size_t idx) {
		double interval, evTime;
		if (this->mu[idx] > 0) {
			interval = rand.generateExponential(this->lambda[idx]);
			StateSimEvent ev;
			ev.msg = MSG_DEL;
			ev.time = this->time + interval;
			ev.sender = senderId(idx, false);
			ev.recv = senderId(idx, false);

			evQueue.push_back(ev);
			std::push_heap(evQueue.begin(), evQueue.end(), std::greater<>{});
		}
	}

	void StateSimMM1k::doStep() {
		int i;

		stepCnt++;
		StateSimEvent ev = evQueue[0];
		std::pop_heap(evQueue.begin(), evQueue.end(), std::greater<>{});
		evQueue.pop_back();

		i = ev.recv;
		this->time = ev.time;
		if (ev.msg == MSG_NEW) {
			this->nPkt[i]++;
			if (ev.sender < 0) {
				this->gen(genIdx(ev.sender));
				pktId++;
			}
			if (this->state[i] == this->k[i] - 1) {
				this->nPktDrop[i]++;
			}
			else {
				double interval = ev.time - this->t;
				this->t = ev.time;
				size_t st = stateIdx(this->state, this->k);
				this->Pi[st] = this->Pi[st] + interval;
				if (this->state[i] == 0) this->serv(i);
				state[i]++;
			}
		}
		else if (ev.msg == MSG_DEL) {
			double interval = ev.time - this->t;
			this->t = ev.time;
			size_t st = stateIdx(this->state, this->k);
			this->Pi[st] = this->Pi[st] + interval;
			state[i]--;
			if (this->state[i] > 0) this->serv(i);
			if (this->outNodeM[i] >= 0) {
				StateSimEvent newEv;
				newEv.recv = senderId(this->outNodeM[i], false);
				newEv.sender = i;
				newEv.time = ev.time;
				newEv.msg = MSG_NEW;
				evQueue.push_back(newEv);
				std::push_heap(evQueue.begin(), evQueue.end(), std::greater<>{});
			}
		}
	}

	void StateSimMM1k::runSteps(size_t maxStep) {
		size_t i;

		std::fill(Pi.begin(), Pi.end(), 0);
		std::fill(nPkt.begin(), nPkt.end(), 0);
		std::fill(nPktDrop.begin(), nPktDrop.end(), 0);
		std::fill(state.begin(), state.end(), 0);
		stepCnt = 0;
		pktId = 0;
		t = 0;
		evQueue.clear();

		for (i = 0; i < lambda.size(); i++) gen(i);
		for (i = 0; i < maxStep; i++) {
			doStep();
		}
	}

	void StateSimMM1k::runEps(double eps, std::vector<double> L0) {
		std::fill(Pi.begin(), Pi.end(), 0);
	}
}