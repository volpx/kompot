#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <fstream>

#include "gsl/gsl_sf_bessel.h"

double j_m1(double x)
{
	return std::cos(x) / x;
}

double j_0(double x)
{
	return std::sin(x) / x;
}

double n_m1(double x)
{
	return std::sin(x) / x;
}

double n_0(double x)
{
	return -std::cos(x) / x;
}

int main()
{
	constexpr uint64_t N = 1e5;
	constexpr double a = 0. + 1. / N;
	constexpr double b = 10.;
	constexpr double h = (b - a) / N;
	constexpr uint64_t lmax = 3;

	std::vector<double> x(N);
	map(x, [a, h](uint64_t i) -> double { return a + i * h; });

	std::vector<double> j((lmax + 2) * N);
	std::vector<double> n((lmax + 2) * N);
	std::vector<double> j_gsl((lmax + 2) * N);
	std::vector<double> n_gsl((lmax + 2) * N);

	// Initialize l=-1
	for (uint64_t i = 0; i < N; i++)
	{
		j[0 * N + i] = j_m1(x[i]);
		// j_gsl[0 * N + i] = gsl_sf_bessel_jl(-1, x[i]);
		n[0 * N + i] = n_m1(x[i]);
		// n_gsl[0 * N + i] = gsl_sf_bessel_yl(-1, x[i]);
	}
	// Initialize l=0
	for (uint64_t i = 0; i < N; i++)
	{
		j[1 * N + i] = j_0(x[i]);
		j_gsl[1 * N + i] = gsl_sf_bessel_jl(0, x[i]);
		n[1 * N + i] = n_0(x[i]);
		n_gsl[1 * N + i] = gsl_sf_bessel_yl(0, x[i]);
	}

	for (uint64_t k = 2; k <= lmax + 1; k++)
	{
		// Computing this l
		uint64_t l = k - 1;
		for (uint64_t i = 0; i < N; i++)
		{
			// Compute j-functions
			j[k * N + i] = (2 * (l - 1) + 1) / x[i] * j[(k - 1) * N + i] - j[(k - 2) * N + i];
			j_gsl[k * N + i] = gsl_sf_bessel_jl(l, x[i]);

			// Compute n-functions
			n[k * N + i] = (2 * (l - 1) + 1) / x[i] * n[(k - 1) * N + i] - n[(k - 2) * N + i];
			n_gsl[k * N + i] = gsl_sf_bessel_yl(l, x[i]);
		}
	}

	std::ofstream file{"data/bessel.dat"};
	file << "#x y0 y_gsl0 n0 n_gsl0 ...\n";
	// Write to file
	for (uint64_t i = 0; i < N; i++)
	{
		file << x[i];
		for (uint64_t k = 1; k <= lmax + 1; k++)
		{
			file << ' ' << j[k * N + i] << ' ' << j_gsl[k * N + i]
				 << ' ' << n[k * N + i] << ' ' << n_gsl[k * N + i];
		}
		file << '\n';
	}
	file.close();

	return 0;
}