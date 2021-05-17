#ifndef gaussian_h
#define gaussian_h

#include <math.h>
#include <iostream>

inline double normalPDF(double x) {
    static constexpr double pi = 3.1415926535897932;
    return 1 / sqrt(2 * pi) * exp(-1.0 / 2.0 * x * x);
}

inline double normalCDF(double x)
{
    return 0.5 * erfc(-x * sqrt(0.5));
}

//    Inverse Normal CDF approximation
inline double invNormalCDF(double u)
{
    // Values used for calculating the result
    static constexpr double a0 = 2.50662823884;
    static constexpr double a1 = -18.61500062529;
    static constexpr double a2 = 41.39119773534;
    static constexpr double a3 = -25.44106049637;

    static constexpr double b0 = -8.47351093090;
    static constexpr double b1 = 23.08336743743;
    static constexpr double b2 = -21.06224101826;
    static constexpr double b3 = 3.13082909833;

    static constexpr double c0 = 0.3374754822726147;
    static constexpr double c1 = 0.9761690190917186;
    static constexpr double c2 = 0.1607979714918209;
    static constexpr double c3 = 0.0276438810333863;
    static constexpr double c4 = 0.0038405729373609;
    static constexpr double c5 = 0.0003951896511919;
    static constexpr double c6 = 0.0000321767881768;
    static constexpr double c7 = 0.0000002888167364;
    static constexpr double c8 = 0.0000003960315187;

    double r, x;
    double y = u - 0.5;
    if (abs(y) < 0.42) {
        r = y * y;
        x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1);
    }
    else {
        r = u;
        if (y > 0) r = 1 - u;
        r = log(-log(r));

        x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
        if (y < 0) x *= -1;
    }

    return x;
}

#endif
