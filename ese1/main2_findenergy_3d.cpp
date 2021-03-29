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

int main()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t M = 1e3; // Number of energies
	constexpr double Ea = 1.5-0.25;
	constexpr double Eb = 0.5 + 10 * 1 + 0.25;
	constexpr uint64_t lmax=3; // Stop value of l

	constexpr uint64_t N = 3e5; // number of mesh points
	constexpr double a = 0.;    // start included
	constexpr double b = 10.;     // stop excluded

	constexpr double h = (b - a) / N; // Mesh spacing
	constexpr double h2 = h * h;

	std::cout << "Problem constants:"
			  << "\nM: " << M
			  << "\nEa: " << Ea
			  << "\nEb: " << Eb
			  << "\nEh: " << (Eb - Ea) / (M - 1)
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << std::endl;

	std::vector<double> vE(M);
	linspace(vE, Ea, Eb);

	std::vector<double> y(N);
	std::vector<double> vy_xmax((lmax+1) * M);
	
	std::vector<double> x(N);
	for (uint64_t i = 0; i < N; i++)
	{
		x[i] = a + i * h;
	}

	for (uint64_t l = 0; l <= lmax; l++)
	{
		// Numerov algorithm
		for (uint64_t j = 0; j < M; j++)
		{

			if (j % 100 == 0)
			{
				std::cout << "l: "<<l<<" P: " << 1.0 * j / M << std::endl;
			}

			// Initialization
			y[0] = 0;
			y[1] = std::pow(x[1], (l + 1));

			double E = vE[j];
			double k_ii2 = 2. * (E - V_eff(x[1], l));
			double k_i2 = 2. * (E - V_eff(x[1], l));
			double k2;

			for (uint64_t i = 2; i < N; i++)
			{

				// New k
				k2 = 2. * (E - V_eff(x[i], l));

				y[i] = (2. * y[i - 1] * (1. - 5. / 12 * h2 * k_i2) - 
					y[i - 2] * (1. + 1. / 12 * h2 * k_ii2)) 
					/ (1. + h2 / 12.0 * k2);

				k_ii2 = k_i2;
				k_i2 = k2;
			}

			vy_xmax[l * M + j] = y[N - 1];
		}
	}

	std::ofstream file{"data/energy2.dat"};
	file << "#E y_xmax(l)\n";
	for (uint64_t i = 0; i < M; i++)
	{
		// Save the result
		file << vE[i] ;
		for (uint64_t l = 0; l <= lmax; l++)
		{
			file << ' ' << vy_xmax[l * M + i];
		}
		file << '\n';
	}
	file.close();

	// Print best estimate of Energy eigenvalues
	for (uint64_t l = 0; l <= lmax; l++)
	{
		uint64_t ii = 0;
		for (uint64_t i = 1; i < M; i++)
		{
			if (vy_xmax[l * M + i] * vy_xmax[l * M + i - 1] < 0)
			{
				std::cout << "Eigenvalue n " << ii
						  << " l: " << l << ": "
						  << (vE[i] - (vE[i] - vE[i - 1]) /
							(vy_xmax[l * M + i] - vy_xmax[l * M + i - 1]))
						  << " \t(" << (2 * ii + l + 1.5) << ")" << std::endl;
				ii++;
			}
		}
	}

	return 0;
}