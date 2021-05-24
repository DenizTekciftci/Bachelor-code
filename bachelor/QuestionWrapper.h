#include <iostream>
#include <iomanip>
#include <random> 
#include "CashOption.h"
#include <numeric>
#include "RNG.h"
#include "math.h"
#include "write.h"
#include <chrono>
#include "stats.h"
#include "integral.h"
#include "simulations.h"
#include <sstream>
#include <tuple>

void log() {}

template<typename First, typename ...Rest>
void log(First&& first, Rest && ...rest)
{
	std::cout << std::forward<First>(first);
	log(std::forward<Rest>(rest)...);
}

// Timing
auto start = std::chrono::steady_clock::now();

// Initial values (Vasicek parameters from Brøgger & Andreasen)
double r0 = -0.005, th = 0.0291, s = 0.0069, k = 0.0325, b_star = 0.0259, n = 50, T = 30, dt = 1 / 12.0, l = 0.02;
double b_hi = 1, b_lo = 0, eps = 0.00001, target = 1;
const int paths = 100000;

// Random number generator
RNG<std::mt19937> rng(time(NULL));

// Cash Option object
CashOption cashOption(th, s, k);

// Calculate ZCB-prices
void Q1() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	int maturities[6] = { 1, 5, 10, 20, 30, 50 };  // Used maturities
	for (int t : maturities) log("P(0, ", t, "): ", cashOption.P(r0, t), "\n"); // Print prices to console

	std::cout << std::endl;
}


// Calculate ZCB-prices using MC and Antithetic sampling
void Q2() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	const int numOfMats = 50;
	std::vector<std::vector<double>> zcb(numOfMats+1, std::vector<double>(paths));
	std::vector<std::vector<double>> summaries1;
	std::vector<std::vector<double>> summaries2;
	std::vector<double> times1;
	std::vector<double> times2;
	// Normal sampling
	for (int mat = 0; mat <= numOfMats; ++mat) {
		// Start timer
		auto t0 = std::chrono::steady_clock::now();

		// Generate normals
		std::vector<double> W(paths);
		rng.genNormal(W);

		// Calculations
		double ekt = exp(-k * mat);
		double mean = mat * th + (r0 - th) / k * (1 - ekt);
		double var = s * s / (k * k * k) * (k * mat + 2 * ekt - 0.5 * exp(-2 * k * mat) - 1.5);

		for (int i = 0; i < paths; ++i) {
			// Reg sampling
			zcb[mat][i] = exp(-(mean + sqrt(var) * W[i]));
		}

		// End timer
		auto t1 = std::chrono::steady_clock::now();
		auto diff = t1 - t0;
		times1.push_back(std::chrono::duration <double, std::milli>(diff).count());

		// Summary
		std::vector<double> summary = Summary(zcb[mat]);
		summaries1.push_back(summary);
		log("P(0, ", mat, "), mean: ", summary[0], ", 95%CI: [", summary[2], ",", summary[3], "]\n");
	}

	// Write to files
	writeList(summaries1, "zcbs", "out21.py");
	writeList(times1, "times", "out212.py");

	// Antithetic sampling
	for (int mat = 0; mat <= numOfMats; ++mat) {
		// Start timer
		auto t0 = std::chrono::steady_clock::now();

		// Generate normals
		std::vector<double> W(paths);
		rng.genNormal(W);

		// Calculations
		double ekt = exp(-k * mat);
		double mean = mat * th + (r0 - th) / k * (1 - ekt);
		double var = s * s / (k * k * k) * (k * mat + 2 * ekt - 0.5 * exp(-2 * k * mat) - 1.5);

		for (int i = 0; i < paths; ++i) {
			// AT sampling
			zcb[mat][i] = (exp(-(mean + sqrt(var) * W[i])) + exp(-(mean - sqrt(var) * W[i]))) / 2;

		}
		// End timer
		auto t1 = std::chrono::steady_clock::now();
		auto diff = t1 - t0;
		times2.push_back(std::chrono::duration <double, std::milli>(diff).count());

		// Summary
		std::vector<double> summary = Summary(zcb[mat]);
		summaries2.push_back(summary);
		log("P(0, ", mat, "), mean: ", summary[0], ", 95%CI: [", summary[2], ",", summary[3], "]\n");
	}

	// Write to files
	writeList(summaries2, "zcbs", "out22.py");
	writeList(times2, "times", "out222.py");

	std::cout << std::endl;
}


// Replicate Figure 3.1              !!!NOT DONE!!!
void Q3() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	std::vector<std::vector<double>> volSurface;
	for (double sigma = 0.00025; sigma <= 0.0075; sigma += 0.000125) {
		// New sigma, new option, new b_star
		//double s2 = 0.0075;
		CashOption cashOption2(th, sigma, k);
		std::cout << sigma << std::endl;

		double b_star2 = findB(1, 0.0, 0.5, 0.00001, cashOption2);

		// Find integral for different interest rates
		for (double i = -0.05; i < 0.025; i += 0.0025) {
			volSurface.push_back({ i, sigma, (C_integral_inf(cashOption2, i, b_star2) - 1) * 100 });
			std::cout << i << ", " << (C_integral_inf(cashOption2, i, b_star2) - 1) * 100 << std::endl;
		}
	}

	writeList(volSurface, "volSurface", "fig31.py");
	std::cout << std::endl;
}

// Replicate Figure 3.2
void Q4() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	double maxSigma = 0.0077;
	double minSigma = 0.0001;
	int numOfPoints = 50;

	std::vector<double> res(numOfPoints);
	std::vector<double>::iterator it = res.begin();

	for (double i = minSigma; i < maxSigma; i += (maxSigma - minSigma) / numOfPoints) {
		CashOption cashOption3(th, i, k);
		double result = 100 * (C_integral_inf(cashOption3, 0, findB(1, -0.02, 0.05, 0.001, cashOption3)) - 1);
		std::cout << i << ": " << result << std::endl;
		*it = result;
		if (!(it == res.end())) {
			it++;
		}
	}

	writeList(res, "test", "out4.py");

	std::cout << std::endl;
}


// Simulate the interest rate, until we hit barrier    !!!NOT DONE!!! missing 95CI
void Q5() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	dt = 1 / 12.0;
	std::tuple<std::vector<double>, std::vector<double>> res = r_to_b_t(cashOption, rng, r0, b_star, paths, dt);
	std::vector<double> integrals = std::get<0>(res);
	std::vector<double> taus = std::get<1>(res);

	writeList(taus, "taus", "out6.py");

	std::cout << std::endl;
}

// Simulate the interest rate, until T
void Q5_2() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	dt = 1 / 12.0;
	// low sigma
	std::vector<std::vector<double>> res = simulate_r_to_t(cashOption, rng, r0, T+20, 5, dt);
	writeList(res, "sims", "out52.py");

	// high sigma
	double s2 = 0.003;
	CashOption cashOption2(th, s2, k);
	std::vector<std::vector<double>> res2 = simulate_r_to_t(cashOption2, rng, r0, T+20, 5, dt);
	writeList(res2, "sims", "out53.py");

	std::cout << std::endl;
}

// Printing in python
void Q6() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;
	std::cout << "Printing in python" << std::endl;

	std::cout << std::endl;
}

// Control variable
void Q7() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	// Start time for regular sampling, only interested in the sample time
	auto t0 = std::chrono::steady_clock::now();
	double time = 200;

	// Sample
	std::vector<double> integrals = r_to_b(cashOption, rng, r0, b_star, paths, dt);
	std::cout << Mean(integrals) << std::endl;
	// Stop time
	auto t1 = std::chrono::steady_clock::now();
	auto diff1 = t1 - t0;

	auto t2 = std::chrono::steady_clock::now();
	double numerator = 0;
	double denominator = 0;

	// Find integrals, taus and control variables
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> cvs = r_to_b_cv_t(cashOption, rng, r0, b_star, paths, dt, time);
	integrals = std::get<0>(cvs);
	std::vector<double> control = std::get<2>(cvs);

	// Means
	double integralMean = Mean(integrals);
	double tauMean = Mean(std::get<1>(cvs));
	double controlMean = Mean(control);
	std::cout << integralMean << ", " << tauMean  <<", "<< controlMean << std::endl;

	// Calculate num and den
	for (int i = 0; i < integrals.size(); i++) {
		numerator += (integrals[i] - integralMean) * (control[i] - controlMean);
		denominator += (control[i] - controlMean) * (control[i] - controlMean);
	}

	// Parameters
	double betaHat = numerator / denominator;
	double analytical = cashOption.P(r0, time);
	std::cout << "Analytical: " << analytical << std::endl;
	// Calculate P_A^CV
	std::vector<double> pacv(paths);
	for (int i = 0; i < paths; i++) {
		pacv[i] = integrals[i] + betaHat * (analytical - control[i]);
	}

	// Stop time for CV
	auto t3 = std::chrono::steady_clock::now();
	auto diff2 = t3 - t2;

	// Print results
	std::cout << "T = " << time << std::endl;
	LogSummary(pacv);
	double pacvMean = Mean(pacv);

	std::cout << "Variance without control variable: " << Variance(integrals, integralMean) / sqrt(paths) << std::endl
		<< "Variance with control variable: " << Variance(pacv, pacvMean) / sqrt(paths) << std::endl
		<< "Efficiency ratio: " << Variance(pacv, pacvMean) / Variance(integrals, integralMean) << std::endl
		<< "T1: " << std::chrono::duration <double, std::milli>(diff1).count()  << ", T2: " << std::chrono::duration <double, std::milli>(diff2).count() << ", ratio: "
		<< std::chrono::duration <double, std::milli>(diff1).count() / std::chrono::duration <double, std::milli>(diff2).count() << std::endl;


	std::cout << std::endl;
}


// Calculate b(t)
void Q8() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	std::vector<std::vector<double>> B;
	std::vector<double> Cash;

	for(double t = 5; t <= 100; t+= 5){
		std::cout << t << std::endl;
		std::vector<double> bs;
		std::vector<double> times;
		for (double u = t; u > 0; u -= dt) {
			double b = findB(1, 0, .1, 0.000001, cashOption, 0, u, n);
			bs.push_back( b );
			times.push_back(t - u);
		}
		double price = cashOption.P(r0, t) + C_integral(cashOption, r0, bs, times, t, n);
		std::cout << C_integral_inf(cashOption, r0, b_star) << std::endl;
		Cash.push_back(price);
		B.push_back(times);
		B.push_back(bs);
	}

	writeList(B, "B", "out8.py");
	writeList(Cash, "C", "out81.py");

}

// Calculate b(t)
void Q8_1() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	std::vector<double> Cash;

	for (double t = 10; t <= 300; t += 20) {
		std::cout << t << std::endl;
		std::vector<double> bs;
		std::vector<double> times;
		for (double u = t; u > 0; u -= dt) {
			double b = findB(1, 0, .1, 0.000001, cashOption, 0, u, n);
			bs.push_back(b);
			times.push_back(t - u);
		}
		double price = cashOption.P(r0, t) + C_integral(cashOption, r0, bs, times, t, n);
		std::cout << price << std::endl;
		Cash.push_back(price);
	}

	writeList(Cash, "C", "out83.py");

}


// Calculate b* 
void Q8_2() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	std::vector<std::vector<double>> B;

	for (double sigma = 0.000125; sigma <= 0.00776; sigma += 0.000125) {
		CashOption cashOption2(th, sigma, k);
		double b = findB(1, 0, .1, 0.000001, cashOption2);
		std::cout << sigma << ", " << b << std::endl;
		B.push_back({ b, sigma });
	}

	writeList(B, "B", "out82.py");
}


// Find cash price with early redemption risk
void ER1() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	double b1 = findB_early_redemption(1, -1, 1, 0.00001, l, cashOption);
	double b2 = findB(1, -1, 1, 0.00001, cashOption);

	std::cout <<"Early redemption (l = " << l * 100 << "%), C(0) = " << C_integral_early_redemption(cashOption, r0, b1, l) << std::endl;

	std::cout << std::endl;
}


// Create data for figure 3.3: The term structure of cash
void ER2() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	// Initialize parameters
	int mat_lo = 1, mat_hi = 51, mat_hi2 = 501, mat_steps = 50, stepSize = (mat_hi - mat_lo) / mat_steps, stepSize2 = (mat_hi2 - mat_lo) / mat_steps;
	std::vector<double> rs = { -0.01, -0.025, -0.05 };

	for (double r : rs) {
		// Data holding matrix
		std::vector<std::vector<double>> prem(mat_steps);
		std::vector<std::vector<double>>::iterator prem_it = prem.begin();

		for (int mat = mat_lo; mat < mat_hi; mat += stepSize) {
			// Lambda
			double l = 1 / (double)mat;

			// Find b and premium
			double b = findB_early_redemption(target, b_lo, b_hi, eps, l, cashOption);
			double prem = 100 * (C_integral_early_redemption(cashOption, r, b, l) - 1);

			// Add slice to matrix and go to next pointer
			*prem_it = { (double)mat, prem };
			prem_it++;
		}

		// Create file name and write to file
		std::stringstream str1;
		str1 << "ER" << std::to_string(int(r * -100)) << ".py";
		writeList(prem, "prem", str1.str());
		std::cout << str1.str() << std::endl;


		// Expected mat between 1 and 500
		prem_it = prem.begin();
		for (int mat = mat_lo; mat < mat_hi2; mat += stepSize2) {
			// Lambda
			double l = 1 / (double)mat;

			// Find b and premium
			double b = findB_early_redemption(target, b_lo, b_hi, eps, l, cashOption);
			double prem = 100 * (C_integral_early_redemption(cashOption, r, b, l) - 1);

			// Add slice to matrix and go to next pointer
			*prem_it = { (double)mat, prem };
			prem_it++;
		}

		// Create file name and write to file
		std::stringstream str2;
		str2 << "ER" << std::to_string(int(r * -100)) << "2.py";
		std::cout << str2.str() << std::endl;
		writeList(prem, "prem", str2.str());
	}

	std::cout << std::endl;
} 


// Creates data for figure 4.1: Cash in fixed supply versus long-dated Treasury STRIPs.
void ER3() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	// Find b's
	double b1 = findB(target, b_lo, b_hi, eps, cashOption);
	double b2 = findB_early_redemption(target, b_lo, b_hi, eps, l, cashOption);
	double b3 = findB_early_redemption(target, b_lo, b_hi, eps, 0.04, cashOption);
	double b4 = findB_early_redemption(target, b_lo, b_hi, eps, 0.08, cashOption);
	double b5 = findB_early_redemption(target, b_lo, b_hi, eps, 0.16, cashOption);

	// Find base prices at r = 0
	double baseP = cashOption.P(0, 30);
	double baseC1 = C_integral_inf(cashOption, 0, b1);
	double baseC2 = C_integral_early_redemption(cashOption, 0, b2, l);
	double baseC3 = C_integral_early_redemption(cashOption, 0, b3, 0.04);
	double baseC4 = C_integral_early_redemption(cashOption, 0, b4, 0.08);
	double baseC5 = C_integral_early_redemption(cashOption, 0, b5, 0.16);

	// Matrix for data
	std::vector<std::vector<double>> res(7);

	for (double i = -0.05; i < 0.051; i += 0.001) {
		// Add slice to matrix
		res[0].push_back(i);
		res[1].push_back((cashOption.P(i, 30) / baseP) * 100);
		res[2].push_back((C_integral_inf(cashOption, i, b1) / baseC1) * 100);
		res[3].push_back((C_integral_early_redemption(cashOption, i, b2, l) / baseC2) * 100);
		res[4].push_back((C_integral_early_redemption(cashOption, i, b3, 0.04) / baseC3) * 100);
		res[5].push_back((C_integral_early_redemption(cashOption, i, b4, 0.08) / baseC4) * 100);
		res[6].push_back((C_integral_early_redemption(cashOption, i, b5, 0.16) / baseC5) * 100);
	}

	// Write matrix to .py
	writeList(res, "res", "fig41.py");

	std::cout << std::endl;
}

void ER4() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	// Initialize parameters
	double b_ZCB = (1 - exp(-k * 30)) / k;
	double dr = 0.00000001;

	// Find b's
	double b = findB(target, b_lo, b_hi, eps, cashOption);
	double bER1 = findB_early_redemption(target, b_lo, b_hi, eps, l, cashOption);
	double bER2 = findB_early_redemption(target, b_lo, b_hi, eps, 0.04, cashOption);
	double bER3 = findB_early_redemption(target, b_lo, b_hi, eps, 0.08, cashOption);
	double bER4 = findB_early_redemption(target, b_lo, b_hi, eps, 0.16, cashOption);

	// Data matrix
	std::vector<std::vector<double>> durs(7);

	for (double r = -0.05; r < 0.051; r += 0.0002) {
		// No early redemption
		std::cout << r << std::endl;

		// No early redemption
		double cash = C_integral_inf(cashOption, r, b);
		double der = (C_integral_inf(cashOption, r + dr, b) - cash) / dr;
		double dur = -1 / cash * der;

		// Early redemption
		double cashER1 = C_integral_early_redemption(cashOption, r, bER1, l);
		double cashER2 = C_integral_early_redemption(cashOption, r, bER1, 0.04);
		double cashER3 = C_integral_early_redemption(cashOption, r, bER1, 0.08);
		double cashER4 = C_integral_early_redemption(cashOption, r, bER1, 0.16);
		double derER1 = (C_integral_early_redemption(cashOption, r + dr, bER1, l) - cashER1) / dr;
		double derER2 = (C_integral_early_redemption(cashOption, r + dr, bER1, 0.04) - cashER2) / dr;
		double derER3 = (C_integral_early_redemption(cashOption, r + dr, bER1, 0.08) - cashER3) / dr;
		double derER4 = (C_integral_early_redemption(cashOption, r + dr, bER1, 0.16) - cashER4) / dr;
		double durER1 = -1 / cashER1 * derER1;
		double durER2 = -1 / cashER2 * derER2;
		double durER3 = -1 / cashER3 * derER3;
		double durER4 = -1 / cashER4 * derER4;

		// Add slice to matrix
		durs[0].push_back(r);
		durs[1].push_back(b_ZCB);
		durs[2].push_back(dur);
		durs[3].push_back(durER1);
		durs[4].push_back(durER2);
		durs[5].push_back(durER3);
		durs[6].push_back(durER4);
	}

	// Write matrix to py
	writeList(durs, "durs", "fig42.py");

	std::cout << std::endl;
}