#ifndef integral_h
#define integral_h

#pragma once
#include "cashOption.h"
#include "gaussian.h"
#include "math.h"
#include <tuple>
#include <time.h>

const double d = 0.02;
const double limit = 10000.0;

double C(CashOption cashOption, const double r, const double b, const double t) {
	double temp_v = cashOption.v(t);
	double temp_f = cashOption.f(r, t);
	double quantile = (temp_f - b) / temp_v;
	return cashOption.P(r, t) * (temp_f * normalCDF(quantile) + temp_v * normalPDF(quantile));
}

double C2(CashOption cashOption, const double r, const double l, const double t) {
	return l * exp(-l * t) * cashOption.P(r, t);
}

double C(CashOption cashOption, const double r, const double b, const double t, const double T) {
	double temp_v = cashOption.v(t, T);
	double temp_f = cashOption.f(r, t, T);
	double quantile = (temp_f - b) / temp_v;
	return cashOption.P(r, t, T) * (temp_f * normalCDF(quantile) + temp_v * normalPDF(quantile));
}


// Calculates the C-integral from 0-inf (10000 in practice)
double C_integral_inf(CashOption cashOption, const double r, const double b) {
	//Initialzie
	double integral = 0;
	double prev = C(cashOption, r, b, d);

	// Loop to find trapezoidal integral
	for (long double i = d; i < limit; i += d) {
		double curr = C(cashOption, r, b, i);
		integral += (prev + curr) / 2 * d;
		prev = curr;
	}

	return integral;
}

// Calculates the C-integral from u-T
double C_integral(CashOption cashOption, const double r, const double b, const double u, const double T, const int n) {
	//Initialize
	double integral = 0;
	double t0 = (u > 0) ? u : 0.000001;
	double prev = C(cashOption, r, b, t0, T);

	// Loop to find trapezoidal integral
	double dt = 1 / (double)n;
	for (double i = t0 + dt; i < T; i += dt) {
		double curr = C(cashOption, r, b, t0, i);
		integral += (prev + curr) / 2 * (1 / (double)n);
		prev = curr;
	}

	return integral;
}


//// Calculates the C-integral from u-T
//double C_integral(CashOption cashOption, const double r, const double b, const double u, const double T, const int n) {
//	//Initialize
//	double integral = 0;
//	double t0 = (u > 0) ? u : 0.000001;
//	double prev = C(cashOption, r, b, t0, T);
//
//	// Loop to find trapezoidal integral
//	for (double i = t0; i < T; i += 1/(double)n) {
//		double curr = C(cashOption, r, b, t0, i);
//		integral += (prev + curr) / 2 * (1 / (double)n);
//		prev = curr;
//	}
//
//	return integral;
//}

// Calculates the C-integral from u-T
double C_integral(CashOption cashOption, const double r, std::vector<std::vector<double>> b_n_i, const double T, const int n) {
	//Initialize
	std::vector<std::vector<double>>::iterator it = b_n_i.begin();
	double integral = 0;
	double t0 = ((*it)[0] > 0) ? (*it)[0] : 0.000001;
	double prev = C(cashOption, r, (*it)[1], t0, T);

	// Loop to find trapezoidal integral
	for (double i = t0; i < T; i += 1 / (double)n) {
		if (((it + 1) != b_n_i.end()) && i >= (*(it + 1))[0]) it++;
		double curr = C(cashOption, r, (*it)[1], t0, i);
		integral += (prev + curr) / 2 * (1 / (double)n);
		prev = curr;
	}

	return integral;
}

// Calculates the C-integral from 0-inf (10000 in practice)
double C_integral_early_redemption(CashOption cashOption, const double r, const double b, const double l) {
	//Initialzie
	double integral1 = 0;
	double integral2 = 0;
	double prev1 = exp(-l * d) * C(cashOption, r, b, d);
	double prev2 = C2(cashOption, r, l, d);

	// Loop to find trapezoidal integral
	for (long double i = d; i < limit; i += d) {
		double curr1 = exp(-l * i) * C(cashOption, r, b, i);
		double curr2 = C2(cashOption, r, l, i);
		integral1 += (prev1 + curr1) / 2 * d;
		integral2 += (prev2 + curr2) / 2 * d;
		prev1 = curr1;
		prev2 = curr2;
	}

	return integral1 + integral2;
}

// Finds b* given th, s and k using bisection
double findB(
	double y_target,         // Target y value
	double a,                // Lower bound
	double b,                // Upper bound
	double epsilon,          // Tolerance
	CashOption cashOption,   // Cash Option object holding parameters
	const double t,
	const double T,
	const int n)
{

	// Choose midpoint and initial y
	double x = 0.5 * (a + b);
	double y = cashOption.P(x,T) + C_integral(cashOption, x, x, 0, T, n);
	//std::cout << y << " " << " " << cashOption.P(x, T) << " " << C_integral(cashOption, x, x, 0, T, n) << std::endl;
	do {
		if (y > y_target) {
			a = x;
		}
		else if (y < y_target) {
			b = x;
		}

		if (x == 0.5 * (a + b)) epsilon = 0.001;
		x = 0.5 * (a + b);
		y = cashOption.P(x, T) + C_integral(cashOption, x, x, 0, T, n);
		//std::cout << x << " " << y << std::endl;
		//std::cout << x << " " <<fabs(y - y_target) << std::endl
	} while (fabs(y - y_target) > epsilon);
	return x;
}

//// Finds b* given th, s and k using bisection
//double findB(
//	double y_target,         // Target y value
//	double a,                // Lower bound
//	double b,                // Upper bound
//	double epsilon,          // Tolerance
//	CashOption cashOption,   // Cash Option object holding parameters
//	const double r0,
//	const double t, 
//	const double T,
//	const int n)   
//	{
//
//	// Choose midpoint and initial y
//	double x = 0.5 * (a + b);
//	double y = C_integral(cashOption, r0, x, t, T, n);
//
//	do {
//		if (y > y_target) {
//			a = x;
//		} else if(y < y_target) {
//			b = x;
//		}
//
//		if (x == 0.5 * (a + b)) epsilon = 0.001;
//		x = 0.5 * (a + b);
//		y = C_integral(cashOption, r0, x, t, T, n);
//		//std::cout << x << " " << y << std::endl;
//		//std::cout << x << " " <<fabs(y - y_target) << std::endl
//	} while (fabs(y - y_target) > epsilon);
//	return x;
//}

// Finds b* given th, s and k using bisection
double findB(
	double y_target,         // Target y value
	double m,                // Lower bound
	double n,                // Upper bound
	double epsilon,          // Tolerance
	CashOption cashOption)   // Cash Option object holding parameters
{

	// Choose midpoint and initial y
	double x = 0.5 * (m + n);
	double y = C_integral_inf(cashOption, x, x);

	do {
		if (y > y_target) {
			m = x;
		}
		else if (y < y_target) {
			n = x;
		}

		x = 0.5 * (m + n);
		y = C_integral_inf(cashOption, x, x);
	} while (fabs(y - y_target) > epsilon);

	return x;
}

// Finds b* given th, s and k using bisection
double findB_early_redemption(
	double y_target,         // Target y value
	double m,                // Lower bound
	double n,                // Upper bound
	double epsilon,          // Tolerance
	double l,          // Tolerance
	CashOption cashOption)   // Cash Option object holding parameters
{

	// Choose midpoint and initial y
	double x = 0.5 * (m + n);
	double y = C_integral_early_redemption(cashOption, x, x, l);

	do {
		if (y > y_target) {
			m = x;
		}
		else if (y < y_target) {
			n = x;
		}

		x = 0.5 * (m + n);
		y = C_integral_early_redemption(cashOption, x, x, l);
	} while (fabs(y - y_target) > epsilon);

	return x;
}

std::tuple<std::vector<double>, std::vector<double>> r_to_b(CashOption cashOption, RNG<std::mt19937> rng, double r0, double b, int paths, double dt) {
	std::vector<double> rs;
	std::vector<double> ts;

	for (int i = 0; i < paths; i++) {
		double r = r0;
		double t = 0;
		double integral = 0;

		while (r < b) {
			double dr = cashOption.get_k() * (cashOption.get_th() - r) * dt + cashOption.get_s() * rng.genNormal(0, sqrt(dt));
			r += dr;
			t += dt;
			integral += r * dt;
		}
		ts.push_back(t);
		rs.push_back(exp(-integral));
	}

	return { rs, ts };
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> r_to_b_cv(CashOption cashOption, RNG<std::mt19937> rng, double r0, double b, int paths, double dt, double time) {
	std::vector<double> rs;
	std::vector<double> ts;
	std::vector<double> cvs;

	for (int i = 0; i < paths; i++) {
		double r = r0;
		double t = 0;
		double integral = 0;
		bool hitBarrier = false;
		bool hitTime = false;

		while (!(hitBarrier && hitTime)) {
			double dr = cashOption.get_k() * (cashOption.get_th() - r) * dt + cashOption.get_s() * rng.genNormal(0, sqrt(dt));
			r += dr;
			t += dt;
			integral += r * dt;

			if (r >= b && !hitBarrier) {
				hitBarrier = true;
				ts.push_back(t);
				rs.push_back(exp(-integral));
			}

			if (t >= time && !hitTime) {
				hitTime = true;
				double cv = exp(-integral);
				if (isnan(cv)) {
					cvs.push_back(0);
				}
				else {
					cvs.push_back(exp(-integral));
				}
			}
		}
	}
	return { rs, ts, cvs };
}


#endif
