module;

#include <random>
#include <thread>
#include <intrin.h>

export module RNG;

import Conf;

namespace RNG {
	export class Generator {
	private:
		std::mt19937_64* genMt19937;
		double expLambda;
	public:
		unsigned __int64 seed;
		Generator();
		double generateExponential(double lambda);
		double generateUniform();
		~Generator();
	};

	unsigned __int64 rdrand() {
		unsigned __int64 result;
		_rdrand64_step(&result);
		return result;
	}

	unsigned __int64 rdseed() {
		unsigned __int64 result;
		_rdseed64_step(&result);
		return result;
	}

	Generator::Generator() {
		seed = rdseed();
		genMt19937 = new std::mt19937_64(seed);
		rdseed();
	}

	Generator::~Generator() {
		delete genMt19937;
	}

	double Generator::generateExponential(double lambda) {
		if (Conf::rngType == Conf::RDRAND) {
			double x;
			x = rdrand(); x = x / 18446744073709551615.0;
			x = std::log(1 - x) / (-expLambda);
			return x;
		}
		else {
			std::exponential_distribution<double> expDist(lambda);
			return expDist(*genMt19937);
		}
	}

	double Generator::generateUniform() {
		if (Conf::rngType == Conf::RDRAND) {
			double x;
			x = rdrand(); x = x / 18446744073709551615.0;
			x = std::log(1 - x) / (-expLambda);
			return x;
		}
		else {
			std::uniform_real_distribution<> uniformDist(0, 1);
			return uniformDist(*genMt19937);
		}
	}
}