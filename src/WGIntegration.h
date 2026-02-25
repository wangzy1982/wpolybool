#ifndef _WG_INTEGRATION_
#define _WG_INTEGRATION_

class WGIntegrand1V {
public:
    virtual ~WGIntegrand1V() {}
    virtual double Calculate(double t) const = 0;
};

double integrate(const WGIntegrand1V* integrand, double a, double b, double epsilon);

#endif
