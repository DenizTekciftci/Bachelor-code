#include <iostream>
#include <chrono>
#include "QuestionWrapper.h"


int main() {
	// Timing
	auto start = std::chrono::steady_clock::now();
	//Q1();
	//Q2();
	//Q3();
	//Q4();
	//Q5();
	//Q5_2();
	//Q6();
	//Q7();
	//Q8();
	Q8_1();
	//Q8_2();
	//ER1();
	//ER2();
	//ER3();
	//ER4();

	// End timing
	std::cout << "--------------------" << std::endl;
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	log(std::chrono::duration <double, std::milli>(diff).count(), " ms\n");
	return 0;
}