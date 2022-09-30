#include "cv.h"
#include <cmath>
#include <vector>

const double factor = 0.025;
const double param_A[] = {-200, -100, -170, 15};
const double param_a[] = {-1, -1, -6.5, 0.7};
const double param_b[] = {0, 0, 11, 0.6};
const double param_c[] = {-10, -10, -6.5, 0.7};
const double param_x_0[] = {1, 0, -0.5, -1};
const double param_y_0[] = {0, 0.5, 1.5, 1};

double subterm(double x, double y, double A, double a, double b, double c, double x_0, double y_0) {
  const double diff_x = x - x_0;
  const double diff_y = y - y_0;
  const double tmp = a * diff_x * diff_x + b * diff_x * diff_y + c * diff_y * diff_y;
  return A * std::exp(tmp);
}

void ds_dxdy(double x, double y, double A, double a, double b, double c, double x_0, double y_0, double& dx, double& dy) {
  const double diff_x = x - x_0;
  const double diff_y = y - y_0;
  const double tmp = a * diff_x * diff_x + b * diff_x * diff_y + c * diff_y * diff_y;
  const double dtmp_dx = 2.0 * a * diff_x + b * diff_y;
  const double dtmp_dy = b * diff_x + 2.0 * c * diff_y;
  dx = A * std::exp(tmp) * dtmp_dx;
  dy = A * std::exp(tmp) * dtmp_dy;
}

double getPotential(double x, double y, double z) {
  double p = 0;
  for (int i = 0; i < 4; ++i) {
    p += subterm(x, y, param_A[i], param_a[i], param_b[i], param_c[i], param_x_0[i], param_y_0[i]);
  }
  return factor * p;
}

std::vector<double> getGradients(double x, double y, double z) {
  std::vector<double> grad(3);
  double dx = 0;
  double dy = 0;
  for (int i = 0; i < 4; ++i) {
    ds_dxdy(x, y, param_A[i], param_a[i], param_b[i], param_c[i], param_x_0[i], param_y_0[i], dx, dy);
    grad[0] += dx;
    grad[1] += dy;
  }
  grad[0] *= factor;
  grad[1] *= factor;
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
