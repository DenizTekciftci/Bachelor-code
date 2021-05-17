#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

void writeList(std::vector<std::vector<double>> vector, std::string variableName, std::string fileName) {
	std::ofstream fileHandler;
	fileHandler.open(fileName);

	fileHandler << variableName << " = [ ";

	for (std::vector<std::vector<double>>::iterator i = vector.begin(); i != vector.end(); i++){
		fileHandler << "[";
		for (std::vector<double>::iterator j = (*i).begin(); j+1 != (*i).end(); j++) {
			fileHandler << *j << ",";
		}
		fileHandler << (*i)[(*i).size()-1] << "]";
		if (i + 1 != vector.end()) fileHandler << ",";
	}
	fileHandler << "]";
	fileHandler.close();

}

void writeList(std::vector<double> vector, std::string variableName, std::string fileName) {
	std::ofstream fileHandler;
	fileHandler.open(fileName);

	fileHandler << variableName << " = [ ";
	for (std::vector<double>::iterator i = vector.begin(); i + 1 != vector.end(); i++) {
		fileHandler << *i << ",";
	}
	fileHandler << vector[vector.size() - 1] << "]";
	fileHandler.close();

}