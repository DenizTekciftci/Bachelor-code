#include "cashOption.h"
#include "gaussian.h"
#include "math.h"
#include <iostream>


CashOption::CashOption(const double th_, const double s_, const double k_)
	:th(th_), s(s_), k(k_) {
}


double CashOption::B(const double t) {
	return (1 - exp(-k * t)) / k;
}


double CashOption::A(const double t) {
	return (th - s * s / (2 * k * k)) * (B(t) - t) - s * s / (4 * k) * B(t) * B(t);
}


double CashOption::P(const double r, const double t) {
	return exp(A(t) - B(t) * r);
}


double CashOption::v(const double t) {
	return sqrt(s * s / (2 * k) * (1 - exp(-2 * k * t)));
}


double CashOption::f(const double r, const double t) {
	return exp(-t * k) * r - (-(exp(-t * k) * (1 - exp(-t * k)) * s * s) / (2 * k * k) + (-1 + exp(-t * k)) * (th - s * s / (2 * k * k)));
}


double CashOption::B(const double t, const double T) {
	return (1 - exp(-k * (T - t))) / k;
}


double CashOption::A(const double t, const double T) {
	return (th - s * s / (2 * k * k)) * (B(t, T) - T + t) - s * s / (4 * k) * B(t, T) * B(t, T);
}


double CashOption::P(const double r, const double t, const double T) {
	return exp(A(t, T) - B(t, T) * r);
}


double CashOption::v(const double t, const double T) {
	return sqrt(s * s / (2 * k) * (1 - exp(-2 * k * (T - t))));
}


double CashOption::f(const double r, const double t, const double T) {
	return exp(-(T - t) * k) * r - (-(exp(-(T - t) * k) * (1 - exp(-(T - t) * k)) * s * s) / (2 * k * k) + (-1 + exp(-(T - t) * k)) * (th - s * s / (2 * k * k)));
}


double CashOption::get_th() const{
	return th;
}


double CashOption::get_s() const {
	return s;
}


double CashOption::get_k() const {
	return k;
}