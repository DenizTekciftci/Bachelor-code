#pragma once

class CashOption {
public:
	CashOption(const double th, const double s, const double k);
	double A(const double t);
	double B(const double t);
	double P(const double r, const double t);
	double v(const double t);
	double f(const double r, const double t);

	// Overloaded to take start and end time
	double A(const double t, const double T);
	double B(const double t, const double T);
	double P(const double r, const double t, const double T);
	double v(const double t, const double T);
	double f(const double r, const double t, const double T);

	// Return privates
	double get_th() const;
	double get_s()  const;
	double get_k()  const;


private:
	double th;    // theta
	double s;     // sigma
	double k;     // kappa
};
