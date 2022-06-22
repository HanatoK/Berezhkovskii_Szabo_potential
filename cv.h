#ifndef CV_H
#define CV_H

#include <vector>

double Ux(double x, double omega, double x0, double beta);
double Uxy(double x, double y, double z = 0);
double dUx_x(double x, double omega, double x0, double beta);
double getPotential(double x, double y, double z = 0);
std::vector<double> getGradients(double x, double y, double z = 0);
std::vector<double> getNumericalGradient(double x, double y, double z = 0, double epsilon = 0.00001);
std::vector<double> getForces(double x, double y, double z = 0);

#endif // CV_H
