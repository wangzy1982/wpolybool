#include "WGIntegration.h"
#include <cmath>

double gauss_high(const WGIntegrand1V* integrand, double a, double b) {
    double d = (b - a) * 0.5;
    double s = (b + a) * 0.5;
    double n0 = d * 0 + s;
    double n1 = d * 0.5384693101056831 + s;
    double n2 = d * -0.5384693101056831 + s;
    double n3 = d * 0.9061798459386640 + s;
    double n4 = d * -0.9061798459386640 + s;
    double w0 = d * 0.5688888888888889;
    double w1 = d * 0.4786286704993665;
    double w2 = w1;
    double w3 = d * 0.2369268850561891;
    double w4 = w3;
    return w0 * integrand->Calculate(n0) +
        w1 * integrand->Calculate(n1) +
        w2 * integrand->Calculate(n2) +
        w3 * integrand->Calculate(n3) +
        w4 * integrand->Calculate(n4);
}

double gauss_low(const WGIntegrand1V* integrand, double a, double b) {
    double d = (b - a) * 0.5;
    double s = (b + a) * 0.5;
    double n0 = d * 0 + s;
    double n1 = d * 0.7745966692414834 + s;
    double n2 = d * -0.7745966692414834 + s;
    double w0 = d * 0.8888888888888888;
    double w1 = d * 0.5555555555555556;
    double w2 = w1;
    return w0 * integrand->Calculate(n0) +
        w1 * integrand->Calculate(n1) +
        w2 * integrand->Calculate(n2);
}

double recursive_integrate(const WGIntegrand1V* integrand, double a, double b, int depth, double epsilon) {
    if (depth >= 32) {
        return gauss_high(integrand, a, b);
    }
    double high = gauss_high(integrand, a, b);
    double low = gauss_low(integrand, a, b);
    if (abs(high - low) <= epsilon) {
        return high;
    }
    double m = (a + b) * 0.5;
    epsilon *= 0.5;
    ++depth;
    return recursive_integrate(integrand, a, m, depth, epsilon) + recursive_integrate(integrand, m, b, depth, epsilon);
}

double integrate(const WGIntegrand1V* integrand, double a, double b, double epsilon) {
    return recursive_integrate(integrand, a, b, 0, epsilon);
}