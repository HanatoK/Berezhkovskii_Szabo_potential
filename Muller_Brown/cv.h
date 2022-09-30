#ifndef CV_H
#define CV_H

#include <vector>

double subterm(double x, double y, double A, double a, double b, double c, double x_0, double y_0);
void ds_dxdy(double x, double y, double A, double a, double b, double c, double x_0, double y_0, double& dx, double& dy);
double getPotential(double x, double y, double z = 0);
std::vector<double> getGradients(double x, double y, double z = 0);
std::vector<double> getNumericalGradient(double x, double y, double z = 0, double epsilon = 0.00001);
std::vector<double> getForces(double x, double y, double z = 0);

#endif // CV_H
