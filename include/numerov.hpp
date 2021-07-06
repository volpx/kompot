#pragma once

// Numerov solver
#include "findzero.hpp"

#include <cstdint>
#include <functional>
#include <vector>
#include <cmath>


// #include <iostream>

// V as function versions
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

// V as computed vector
void numerov_integrate(
	double y[],
	const double h, 
	const uint64_t N, const double V[],
	const double E);

double numerov_integrate_yxmax(
	double y_ii,double y_i,
	const double h,
	const uint64_t N, const double V[],
	const double E);

double numerov_find_energy(
	double y_ii,double y_i,
	const double h,
	const uint64_t N, const double V[],
	const double Ea, const double Eb, const double Eh);