#include "cv.h"
#include <cmath>
#include <vector>


double Ux(double x, double omega, double x0, double beta) {
  const double Delta = omega * omega * x0 * x0 / 4.0;
  if (x <= -0.5 * x0) {
    return (-Delta + omega * omega * (x + x0) * (x + x0) / 2.0) / beta;
  } else if ((x > -0.5 * x0) && (x < 0.5 * x0)) {
    return (-omega * omega * x * x / 2.0) / beta;
  } else {
    return (-Delta + omega * omega * (x - x0) * (x - x0) / 2.0) / beta;
  }
}

double Uxy(double x, double y, double z) {
  const double omega = 2.0;
  const double big_omega_square = 1.01 * omega * omega;
  const double x0 = 2.2;
  const double beta = 1.0 / (300.0 * 0.0019872041);
  const double ux = Ux(x, omega, x0, beta);
  return ux + (big_omega_square * (x - y) * (x - y) / 2.0) / beta;
}

double dUx_x(double x, double omega, double x0, double beta) {
  if (x <= -0.5 * x0) {
    return (omega * omega * (x + x0)) / beta;
  } else if ((x > -0.5 * x0) && (x < 0.5 * x0)) {
    return (-omega * omega * x) / beta;
  } else {
    return (omega * omega * (x - x0)) / beta;
  }
}

double getPotential(double x, double y, double z) {
    return Uxy(x, y, z);
}

std::vector<double> getGradients(double x, double y, double z) {
  std::vector<double> grad(3);
  const double omega = 2.0;
  const double big_omega_square = 1.01 * omega * omega;
  const double x0 = 2.2;
  const double beta = 1.0 / (300.0 * 0.0019872041);
  const double dUxy_dx = dUx_x(x, omega, x0, beta) + (big_omega_square * (x - y)) / beta;
  const double dUxy_dy = -(big_omega_square * (x - y)) / beta;
  grad[0] = dUxy_dx;
  grad[1] = dUxy_dy;
  grad[2] = 0;
  return grad;
}

std::vector<double> getNumericalGradient(double x, double y, double z, double epsilon) {
  const double dV_dx = (getPotential(x + epsilon, y, z) - getPotential(x - epsilon, y, z)) / (2.0 * epsilon);
  const double dV_dy = (getPotential(x, y + epsilon, z) - getPotential(x, y - epsilon, z)) / (2.0 * epsilon);
  const double dV_dz = (getPotential(x, y, z + epsilon) - getPotential(x, y, z - epsilon)) / (2.0 * epsilon);
  return std::vector<double>({dV_dx, dV_dy, dV_dz});
}

std::vector<double> getForces(double x, double y, double z) {
  auto grad = getGradients(x, y, z);
  grad[0] *= -1.0;
  grad[1] *= -1.0;
  grad[2] *= -1.0;
  return grad;
}
