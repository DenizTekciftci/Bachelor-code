#pragma once
#include <vector>
#include "math.h"

double Mean(std::vector<double> vector) {
	return accumulate(vector.begin(), vector.end(), 0.0) / vector.size();
}

double Variance(std::vector<double> vector, double mean) {
	double SSE = 0;
	for (std::vector<double>::iterator i = vector.begin(); i != vector.end(); ++i) {
		SSE += (*i - mean) * (*i - mean);
	}
	return SSE / (double)vector.size();
}

double Std(std::vector<double> vector, double mean) {
	return sqrt(Variance(vector, mean));
}

std::vector<double> Summary(std::vector<double> vector) {
	double mean = Mean(vector);
	double std = Std(vector, mean);
	return 	std::vector<double>  {(double)mean, std, (double)mean - 1.96 * std / sqrt(vector.size()), (double)mean + 1.96 * std / sqrt(vector.size())};
}

void LogSummary(std::vector<double> vector) {
	double mean = Mean(vector);
	double std = Std(vector, mean);
	std::cout << "Mean: " << mean << ", Std: " << std << ", 95%CI: [" << mean - 1.96 * std / sqrt(vector.size()) << ", " << mean + 1.96 * std / sqrt(vector.size()) << "]" << std::endl;

}

std::vector<double> CI(std::vector<double> vector, double standardDeviations) {
	double mean = Mean(vector);
	double std = Std(vector, mean);
	return 	std::vector<double>  {mean - standardDeviations * std / sqrt(vector.size()), mean + standardDeviations * std / sqrt(vector.size())};
}