module;

#include <iostream>
#include <fstream>
#include <filesystem>

export module Conf;

import Json;

void createDirIfNotExisted(std::filesystem::path const& path) {
	if (!std::filesystem::exists(path))
		std::filesystem::create_directory(path);
	else if (!std::filesystem::is_directory(path)) {
		std::cout << "Error: " << path << " existed but isn't a directory\n";
		std::cout << "Ensure that " << path << " is directory before running\n";
		exit(-1);
	}
}

export namespace Conf {
	enum RNGType {
		MersenneTwister, RDRAND
	};

	RNGType rngType = MersenneTwister;
	bool reuseGeneratedNumbers = false;

	std::filesystem::path workingPath = ".";
	std::string workingMode = "interactive";
	std::string batchFile = "";

	double gamma = 1.0;

	void initProgram() {
		std::string str;
		std::filesystem::path filename;
		std::ifstream conf;
		Json::Reader jsonReader;
		Json::Object* json;

		filename = "config.json";
		conf.open(filename);
		jsonReader.setInput(&conf);
		json = jsonReader.parse();
		if (json != nullptr) {
			if (json->strfields.contains("workingPath")) {
				str = json->strfields["workingPath"];
				if (std::filesystem::exists(str) && std::filesystem::is_directory(str)) {
					workingPath = str;
				}
				else std::cout << "Working path is not valid in config.json file\n";
			}
			if (json->strfields.contains("workingMode")) {
				workingMode = json->strfields["workingMode"];
			}
			if (json->strfields.contains("batchFile")) {
				batchFile = json->strfields["batchFile"];
			}
		}
		delete json;
		conf.close();

		createDirIfNotExisted(workingPath / "result");
		createDirIfNotExisted(workingPath / "backup");
		createDirIfNotExisted(workingPath / "sorted result");
	}
}