#include "cv.h"

#include <iostream>
#include <fstream>
#include <cmath>

double rmsd(const std::vector<double>& x, const std::vector<double>& y) {
  double sum = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
      sum += (x[i] - y[i]) * (x[i] - y[i]);
  }
  sum /= (double)x.size();
  return std::sqrt(sum);
}

void test_grad(const std::vector<double>& p) {
  const auto analytical_diff = getGradients(p[0], p[1], p[2]);
  double epsilon = 0.1;
  std::cout << "Test gradient at " << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
  for (int i = 0; i < 5; ++i) {
    const auto numerical_diff = getNumericalGradient(p[0], p[1], p[2], epsilon);
    std::cout << "Epsilon = " << epsilon << " ; ";
    std::cout << "Error = " << rmsd(analytical_diff, numerical_diff) << std::endl;
    epsilon /= 10.0;
  }
}

void dump_potential() {
  std::ofstream ofs("potential.dat");
  double x_lower = -2;
  double x_upper = 1;
  double x_width = 0.01;
  size_t nx = std::nearbyint((x_upper - x_lower) / x_width);
  double y_lower = -1;
  double y_upper = 2;
  double y_width = 0.01;
  size_t ny = std::nearbyint((y_upper - y_lower) / y_width);
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      double x = x_lower + (i + 0.5) * x_width;
      double y = y_lower + (j + 0.5) * y_width;
      double potential = getPotential(x, y);
      ofs << x << " " << y << " " << potential << std::endl;
    }
  }
}

int main() {
  std::vector<double> p{0.0, -0.05, 0};
  test_grad(p);
  p[0] = -1.0;
  test_grad(p);
  p[1] = 1.0;
  test_grad(p);
  dump_potential();
}
