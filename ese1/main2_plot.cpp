#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_hermite.h>

double V(double x)
{
	return 0.5 * x * x;
}
double V_eff(double x, double l)
{
	return 0.5 * x * x + 0.5 * (l * (l + 1)) / (x * x);
}

double y_true(double x, int n)
{
	return gsl_sf_hermite(n, x) * std::exp(-0.5 * x * x);
}

void save_function(const double E, const int l, const double a, const double b, const uint64_t N, std::ofstream &file)
{
	const double h = (b - a) / N;
	const double h2 = h * h;

	double y = 0, y_i, y_ii;
	double k2, k_i2, k_ii2;

	y_ii = std::pow(a, (l + 1));
	y_i = std::pow(a + h, (l + 1));
	k_ii2 = 2. * (E - V_eff(a, l));
	k_i2 = 2. * (E - V_eff(a + h, l));

	file << a << ' ' << y_ii << '\n';
	file << a + h << ' ' << y_i << '\n';

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V_eff(a + i * h, l));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);

		file << a + i * h << ' ' << y << '\n';

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
}

int main()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t N = 3e5; // number of mesh points
	constexpr double a = 0.01;    // start included
	constexpr double b = 10.;     // stop excluded

	constexpr double h = (b - a) / N; // Mesh spacing
	constexpr double h2 = h * h;

	constexpr uint32_t nmax=3;

	std::cout << "Problem constants:"
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << std::endl;

	// Create an output of the wavefunction
	char filename[50];
	sprintf(filename, "data/eigen_mesh2e2.dat");
	std::ofstream file{filename};
	save_function(1.5, 0, a, b, 200, file);
	file.close();

	return 0;
}