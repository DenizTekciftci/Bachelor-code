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
double r0 = -0.005, th = 0.0291, s = 0.0069, k = 0.0325, b_star = 0.0259, n = 50, T = 30, dt = 1 / 52.0, l = 0.02;
double b_hi = 1, b_lo = -1, eps = 0.00001, target = 1;
const int paths = 100000;

// Random number generator
RNG<std::mt19937> rng(time(NULL));

// Cash Option object
CashOption cashOption(th, s, k);

// Swap prices from 07/05/2021; https://fred.stlouisfed.org/release/tables?rid=444&eid=783790&od=2021-05-07
std::vector<int> maturities = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30 };
std::vector<double> swapPrices = { 0.00189, 0.00243, 0.00413, 0.00639, 0.00852, 0.01045, 0.01209, 0.01341, 0.01448, 0.01538,
								  0.01813, 0.01926, 0.01981 };

double temp_th = th;
double temp_s = s;
double temp_k = k;

// Calibration
void calibrate(){
	double minSSE = 100;
	double r0_star = 0;
	double th_star = 0;
	double s_star = 0;
	double k_star = 0;
	double resolution = 0.0001;

	for (double t_r0 = -0.02; t_r0 < 0.01; t_r0 += resolution) {
		for (double t_th = -0.01; t_th < 0.04; t_th += resolution) {
			for (double t_k = 0.01; t_k < 0.045; t_k += resolution) {
				for (double t_s = 0.00000001; t_s < 0.02; t_s += resolution) {
					
					CashOption cashOption2(t_th, t_s, t_k);
					double SSE = 0;
					std::vector<double> Ps; // ZCB-prices

					for (int t = 1; t < T + 1; t++) {
						Ps.push_back(cashOption2.P(t_r0, t));
					}

					for (int i = 0; i < maturities.size(); i++) {
						double sum = 0;
						// Because index 0 of Ps = 1, j = 0
						for (int j = 0; j < maturities[i]; j++) {
							// Fixed leg is discounted
							sum += Ps[j];
						}
						SSE += (swapPrices[i] * sum / maturities[i] - (1 - Ps[maturities[i] - 1])) * (swapPrices[i] * sum / maturities[i] - (1 - Ps[maturities[i] - 1]));
					}
					if (SSE < minSSE) {
						minSSE = SSE;
						r0_star = t_r0;
						th_star = t_th;
						s_star = t_s;
						k_star = t_k;
					}
				}
			}
		}
	}


	std::cout << "SSE: " << minSSE << ", r0 = " << r0_star << ", th = " << th_star << ", s = " << s_star << ", k = " << k_star << std::endl;

	CashOption cashOption2(th_star, s_star, k_star);
	std::vector<double> Ps; // ZCB-prices

	for (int t = 1; t < T + 1; t++) {
		Ps.push_back(cashOption2.P(r0_star, t));
	}

	for (int i = 0; i < maturities.size(); i++) {
		double sum = 0;
		// Because index 0 of Ps = 1, j = 0
		for (int j = 0; j < maturities[i]; j++) {
			// Fixed leg is discounted
			sum += Ps[j];
		}
		std::cout << swapPrices[i] * sum / maturities[i] << ", " << 1 - Ps[maturities[i] - 1] << std::endl;
	}

	// returns 
	// r0 = 0.002
	// th = -0.0006
	// s = 0.0018
	// k = 0.0438
}

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
	std::vector<std::vector<double>> zcb(numOfMats, std::vector<double>(paths));

	// Normal sampling
	for (int mat = 1; mat < numOfMats + 1; ++mat) {
		std::vector<double> W(paths);
		rng.genNormal(W);

		double ekt = exp(-k * mat);
		double mean = mat * th + (r0 - th) / k * (1 - ekt);
		double var = s * s / (k * k * k) * (k * mat + 2 * ekt - 0.5 * exp(-2 * k * mat) - 1.5);

		for (int i = 0; i < paths; ++i) {
			zcb[mat - 1][i] = exp(-(mean + sqrt(var) * W[i]));
		}
		std::vector<double> summary = Summary(zcb[mat - 1]);
		log("P(0, ", mat, "), mean: ", summary[0], ", 95%CI: [", summary[2], ",", summary[3], "]\n");

	}

	// Antithetic sampling
	for (int mat = 1; mat < numOfMats + 1; ++mat) {
		std::vector<double> W(paths);
		rng.genNormal(W);

		double ekt = exp(-k * mat);
		double mean = mat * th + (r0 - th) / k * (1 - ekt);
		double var = s * s / (k * k * k) * (k * mat + 2 * ekt - 0.5 * exp(-2 * k * mat) - 1.5);

		for (int i = 0; i < paths; ++i) {
			zcb[mat - 1][i] = (exp(-(mean + sqrt(var) * W[i])) + exp(-(mean - sqrt(var) * W[i]))) / 2;

		}
		std::vector<double> summary = Summary(zcb[mat - 1]);
		log("P(0, ", mat, "), mean: ", summary[0], ", 95%CI: [", summary[2], ",", summary[3], "]\n");

	}

	writeList(zcb, "test", "out2.py");

	std::cout << std::endl;
}


// Replicate Figure 3.1              !!!NOT DONE!!!
void Q3() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	// New sigma, new option, new b_star
	double s2 = 0.0075;
	CashOption cashOption2(th, s2, k);
	double b_star2 = findB(1, -0.02, 0.05, 0.00001, cashOption);

	// Find integral for different interest rates
	for (double i = -0.1; i < 0.09; i += 0.02) {
		std::cout << i << " " << (C_integral_inf(cashOption, i, b_star2) - 1) * 100 << std::endl;
	}

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
	std::tuple<std::vector<double>, std::vector<double>> res = r_to_b(cashOption, rng, r0, b_star, paths, dt);
	std::vector<double> integrals = std::get<0>(res);
	std::vector<double> taus = std::get<1>(res);

	writeList(taus, "taus", "out6.py");

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

	// Initialize parameters
	double time = 200;
	double numerator = 0;
	double denominator = 0;

	// Find integrals, taus and control variables
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> res = r_to_b_cv(cashOption, rng, r0, b_star, paths, dt, time);
	std::vector<double> integrals = std::get<0>(res);
	std::vector<double> taus = std::get<1>(res);
	std::vector<double> cvs = std::get<2>(res);

	// Means
	double integralMean = Mean(integrals);
	double controlMean = Mean(cvs);

	// Calculate num and den
	for (int i = 0; i < integrals.size(); i++) {
		numerator += (integrals[i] - integralMean) * (cvs[i] - controlMean);
		denominator += (cvs[i] - controlMean) * (cvs[i] - controlMean);
	}

	// Parameters
	double betaHat = numerator / denominator;
	double analytical = cashOption.P(r0, time);
	std::vector<double> pcv(paths);

	// Calculate control variables
	for (int i = 0; i < pcv.size(); i++) {
		pcv[i] = integrals[i] + betaHat * (analytical - cvs[i]);
	}

	// Print results
	std::cout << "r(0) = " << r0 << ": ";
	LogSummary(pcv);
	double pcvMean = Mean(pcv);
	std::cout << "Variance without control variable: " << Variance(integrals, integralMean) << std::endl
		<< "Variance with control variable: " << Variance(pcv, pcvMean) << std::endl
		<< "Efficiency ratio: " << Variance(integrals, integralMean) / Variance(pcv, pcvMean) << std::endl;

	std::cout << std::endl;
}


// Calculate b(t) then find C(0, T)
void Q8() {
	std::cout << "--------- " << __FUNCTION__ << " ---------" << std::endl;

	std::vector<std::vector<double>> bs;
	double r = -0.05;
	T = 100;
	for (double u = T; u > 0; u -= dt*12) {
		double b = findB(1, 0, .1, 0.000001, cashOption, 0, u, n*10);
		std::cout << u << ", " << b << std::endl;
		bs.push_back({ T-u , b });
	}

	writeList(bs, "bs", "out8.py");

	std::cout << "C(0, " << T << ") = " << cashOption.P(r, T) + C_integral(cashOption, r, bs, T, n) << std::endl;
	std::cout << "Cinf: " << C_integral_inf(cashOption, r, findB(1, 0, 0.1, 0.000001, cashOption)) << std::endl;

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

	// Find base prices at r = 0
	double baseP = cashOption.P(0, 30);
	double baseC1 = C_integral_inf(cashOption, 0, b1);
	double baseC2 = C_integral_early_redemption(cashOption, 0, b2, l);

	// Matrix for data
	std::vector<std::vector<double>> res(4);

	for (double i = -0.05; i < 0.051; i += 0.001) {
		// Add slice to matrix
		res[0].push_back(i);
		res[1].push_back((cashOption.P(i, 30) / baseP) * 100);
		res[2].push_back((C_integral_inf(cashOption, i, b1) / baseC1) * 100);
		res[3].push_back((C_integral_early_redemption(cashOption, i, b2, l) / baseC2) * 100);
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
	double bER = findB_early_redemption(target, b_lo, b_hi, eps, l, cashOption);

	// Data matrix
	std::vector<std::vector<double>> durs(4);

	for (double r = -0.05; r < 0.051; r += 0.0005) {
		// No early redemption
		double cash = C_integral_inf(cashOption, r, b);
		double der = (C_integral_inf(cashOption, r + dr, b) - cash) / dr;
		double dur = -1 / cash * der;

		// Early redemption
		double cashER = C_integral_early_redemption(cashOption, r, bER, l);
		double derER = (C_integral_early_redemption(cashOption, r + dr, bER, l) - cashER) / dr;
		double durER = -1 / cashER * derER;

		// Add slice to matrix
		durs[0].push_back(r);
		durs[1].push_back(b_ZCB);
		durs[2].push_back(dur);
		durs[3].push_back(durER);
	}

	// Write matrix to py
	writeList(durs, "durs", "fig42.py");

	std::cout << std::endl;
}