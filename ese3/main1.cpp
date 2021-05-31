#include "uniconst.hpp"
#include "vector_help.hpp"
#include "HF.hpp"
#include "MF.hpp"
#include "gsl_help.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

// #include <eigen3/Eigen/Dense>
// #define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

double E_c(double r_s)
{
	constexpr double p = 1.0;
	constexpr double A = 0.031091;
	constexpr double alpha_1 = 0.21370;
	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};
	double DEN = 2 * A * (beta[0] * std::pow(r_s, 0.5) + beta[1] * r_s + beta[2] * std::pow(r_s, 3.0 / 2) + beta[3] * std::pow(r_s, p + 1));
	return -2 * A * (1 + alpha_1 * r_s) * std::log(1 + 1.0 / DEN);
}

double v_ext(const double r_s){
	double rho_b = 1.0 / (4.0 / 3 * M_PI * r_s);
	return 2*M_PI*1;
}
int main()
{
	std::cout << "DFT eseseses" << std::endl;

	return 0;
}
