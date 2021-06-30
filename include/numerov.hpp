#pragma once

// Numerov solver
#include "findzero.hpp"

#include <cstdint>
#include <functional>
#include <vector>
#include <cmath>


// #include <iostream>

void numerov_integrate(
	double y[],
	const double a, const double h, 
	const uint64_t N, std::function<double(double)> V,
	const double E);

double numerov_integrate_yxmax(
	double y0, double y1,
	const double a, const double h,
	const uint64_t N, std::function<double(double)> V,
	const double E);

double numerov_find_energy(
	double y0, double y1,
	const double a, const double h,
	const uint64_t N, std::function<double(double)> V,
	const double Ea, const double Eb, const double Eh);