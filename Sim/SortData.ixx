module;

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <map>
#include <string>

export module SortData;

import Json;
import Conf;

bool getArrayOfDouble(Json::Object* json, std::string fieldname, std::vector<double>& res) {
	bool error = false;
	size_t i, len;
	Json::Object* array = nullptr;

	if (json->objfields.contains(fieldname)) array = json->objfields[fieldname];
	else error = true;
	if (!error) {
		len = array->idx.size();
		error = (len != res.size());
	}
	if (!error) {
		for (i = 0; i < len; i++) {
			if (array->idx[i] == Json::Type::number)
				res[i] = array->numfields[std::to_string(i)];
			else error = true;
		}
	}
	return error;
}

double calcDelta(std::vector<std::vector<double>>& lambda, std::vector<std::vector<double>>& mu, size_t idx) {
	size_t nNodes, i;
	double sumLambda, value, t;

	nNodes = lambda.size();
	value = (mu[0][idx] - lambda[0][idx]) / mu[0][idx];
	sumLambda = lambda[0][idx];
	i = 1;
	while (i < nNodes) {
		sumLambda = sumLambda + lambda[i][idx];
		t = (mu[i][idx] - sumLambda) / mu[i][idx];
		if (t < value) value = t;
		i++;
	}
	return value;
}

void swap(
	std::vector<std::vector<double>>& lambda,
	std::vector<std::vector<double>>& mu,
	std::vector<double>& delta,
	std::vector<double>& meanServTime,
	std::vector<double>& meanTimeInSys,
	size_t x, size_t y
) {
	double t3, t4, t5;
	size_t i;
	for (i = 0; i < lambda.size(); i++) {
		double t, t2;
		t = lambda[i][x]; lambda[i][x] = lambda[i][y]; lambda[i][y] = t;
		t2 = mu[i][x]; mu[i][x] = mu[i][y]; mu[i][y] = t2;
	}
	t3 = delta[x]; delta[x] = delta[y]; delta[y] = t3;
	t4 = meanServTime[x]; meanServTime[x] = meanServTime[y]; meanServTime[y] = t4;
	t5 = meanTimeInSys[x]; meanTimeInSys[x] = meanTimeInSys[y]; meanTimeInSys[y] = t5;
}

void myqsort(
	std::vector<std::vector<double>>& lambda,
	std::vector<std::vector<double>>& mu,
	std::vector<double>& delta,
	std::vector<double>& meanServTime,
	std::vector<double>& meanTimeInSys,
	size_t beg, size_t end
) {
	size_t pivot, i, j;
	//rng = np.random.default_rng()
	//pivot = rng.integers(beg, end + 1)
	pivot = (beg + end) / 2;
	i = beg; j = end;
	while (i < pivot || j > pivot) {
		if (i < pivot) {
			while (i < pivot) {
				if (delta[i] > delta[pivot]) {
					swap(lambda, mu, delta, meanServTime, meanTimeInSys, i, pivot);
					pivot = i;
				}
				else i++;
			}
		}
		else {
			while (j > pivot) {
				if (delta[j] < delta[pivot]) {
					swap(lambda, mu, delta, meanServTime, meanTimeInSys, j, pivot);
					pivot = j;
				}
				else j--;
			}
		}
	}

	if (pivot - beg > 2) {
		myqsort(lambda, mu, delta, meanServTime, meanTimeInSys, beg, pivot - 1);
	}
	else if (pivot - beg == 2) {
		if (delta[beg] > delta[beg + 1])
			swap(lambda, mu, delta, meanServTime, meanTimeInSys, beg, beg + 1);
	}

	if (end - pivot > 2) {
		myqsort(lambda, mu, delta, meanServTime, meanTimeInSys, pivot + 1, end);
	}
	else if (end - pivot == 2) {
		if (delta[end - 1] > delta[end])
			swap(lambda, mu, delta, meanServTime, meanTimeInSys, end - 1, end);
	}
}

export namespace SortData {
	void sortByDeltaJson(Json::Object* json, std::string fname) {
		size_t nData, nNodes, i;
		std::time_t startTime, endTime;
		bool error = false;
		std::vector<std::vector<double>> lambda, mu;
		std::vector<double> meanServTime, meanTimeInSys, delta;
		Json::Writer writer;
		std::ofstream out;

		if (json->numfields.contains("nData")) nData = json->numfields["nData"];
		else error = true;
		if (json->numfields.contains("nNodes")) nNodes = json->numfields["nNodes"];
		else error = true;
		
		if (!error) {
			lambda.resize(nNodes); mu.resize(nNodes);
			delta.resize(nData);
			meanServTime.resize(nData); meanTimeInSys.resize(nData);

			for (i = 0; i < nNodes; i++) {
				lambda[i].resize(nData); mu[i].resize(nData);
				error = getArrayOfDouble(json, "lambda"+std::to_string(i), lambda[i]);
				error = getArrayOfDouble(json, "mu"+std::to_string(i), mu[i]);
			}

			error = getArrayOfDouble(json, "meanServTime", meanServTime);
			error = getArrayOfDouble(json, "meanTimeInSys", meanTimeInSys);
		}
		if (!error) {
			for (i = 0; i < nData; i++) delta[i] = calcDelta(lambda, mu, i);
			std::cout << std::format("Sorting {} ... ", fname);
			std::time(&startTime);
			myqsort(lambda, mu, delta, meanServTime, meanTimeInSys, 0, nData - 1);
			std::time(&endTime);
			std::cout << std::format("in {}s\n", endTime - startTime);

			out.open(Conf::workingPath / "sorted result" / fname);
			writer.setOutput(&out);
			writer.openObject(1, "");
			writer.writeDouble("nData", nData);
			writer.writeDouble("nNodes", nNodes);
			for (i = 0; i < nNodes; i++) {
				writer.writeArrayOfDouble("lambda" + std::to_string(i), lambda[i].data(), nData);
				writer.writeArrayOfDouble("mu" + std::to_string(i), mu[i].data(), nData);
			}
			writer.writeArrayOfDouble("delta", delta.data(), nData);
			writer.writeArrayOfDouble("meanServTime", meanServTime.data(), nData);
			writer.writeArrayOfDouble("meanTimeInSys", meanTimeInSys.data(), nData);
			writer.closeObject();
			out.close();
		}

		delete json;
	}

	void sortByDelta() {
		std::filesystem::path path(Conf::workingPath / "result");
		std::filesystem::directory_entry entry;
		std::filesystem::directory_iterator iter(path);

		while (!iter._At_end()) {
			entry = iter.operator*();
			if (entry.is_regular_file() && entry.path().extension() == ".json") {
				std::ifstream file;
				Json::Object* json;
				Json::Reader jsonReader;

				file.open(entry.path());
				jsonReader.setInput(&file);
				json = jsonReader.parse();
				file.close();

				if (json != nullptr) sortByDeltaJson(json, entry.path().filename().string());
			}
			iter++;
		}
	}
}