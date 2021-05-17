#pragma once
#include "gaussian.h"

template <typename T>
class RNG {
public:

    // Constructor seeding with unsigned int
    RNG(const unsigned int seedValue) {
        rng.seed(seedValue);
    }

    // Overloaded constructor
    RNG(const T& generator) {
        rng = generator;
    }

    // Gen uniform realization
    double genUniform() {
        return rng() / (double)rng.max();
    }

    auto gen() {
        return rng();
    }

    double genNormal() {
        return invNormalCDF(genUniform());
    }

    // Overloaded genNormal to produce non-standard normally distributed random realizations
    double genNormal(const double mean, const double std) {
        return invNormalCDF(genUniform()) * std + mean;
    }

    void genNormal(std::vector<double> &realizations) {
        for (std::vector<double>::iterator i = realizations.begin(); i != realizations.end(); i++) {
            *i = genNormal();
        }
    }

    void genNormal(std::vector<double>& realizations, double mean, double std) {
        for (std::vector<double>::iterator i = realizations.begin(); i != realizations.end(); i++) {
            *i = genNormal(mean, std);
        }
    }

    auto max() {
        return rng.max();
    }

    auto min() {
        return rng.min();
    }

    // Seed the generator
    void setSeed(const unsigned int seedValue) {
        rng.seed(seedValue);
    }

private:
    T rng;
};