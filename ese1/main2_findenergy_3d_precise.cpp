#include "vector_help.hpp"
#include "findzero.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>

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

// Function that returns the value of y(xmax) for given E and l
double F(const double E, const int l,const double a, const double b, const uint64_t N)
{
	const double h=(b-a)/N;
	const double h2=h*h;

	double y=0,y_i,y_ii;
	double k2,k_i2,k_ii2;

	y_ii = std::pow(a , (l + 1));
	y_i = std::pow(a + h, (l + 1));
	k_ii2 = 2. * (E - V_eff(a , l));
	k_i2 = 2. * (E - V_eff(a + h, l));

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V_eff(a+i*h, l));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
				y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			   (1. + h2 / 12.0 * k2);

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii=y_i;
		y_i=y;
	}
	return y;
}

void save_function(const double E, const int l,const double a, const double b,const uint64_t N, std::ofstream &file){
	const double h = (b - a) / N;
	const double h2 = h * h;

	double y = 0, y_i, y_ii;
	double k2, k_i2, k_ii2;

	y_ii = std::pow(a , (l + 1));
	y_i = std::pow(a + h, (l + 1));
	k_ii2 = 2. * (E - V_eff(a , l));
	k_i2 = 2. * (E - V_eff(a + h, l));

	file << a << ' ' << y_ii<<'\n';
	file << a+h << ' ' << y_i << '\n';

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V_eff(a + i * h, l));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);

		file << a + i*h << ' ' << y << '\n';

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
}

int main()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t M = 1e3; // Number of energies
	constexpr double Ea = 1.5-0.25;
	constexpr double Eb = 0.5 + 10 * 1 + 0.25;
	constexpr uint64_t lmax=3; // Stop value of l

	constexpr uint64_t N = 3e5; // number of mesh points
	constexpr double a = 0.+1./N;    // start included
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
			// Search for a change of sign
			if (vy_xmax[l * M + i] * vy_xmax[l * M + i - 1] < 0)
			{
				// Get a better value by secants method
				double E_after_secants=findzero_secants_xeps(
					[a,b,N,l](double E)->double{return F(E,l,a,b,N);},
					vE[i],vE[i-1],1e-7,vE[i-1],vE[i]);
				double E_true=(2 * ii + l + 1.5);
				std::cout << "Eigenvalue k " << ii
						  << " l: " << l << ": "
						  << (vE[i] - (vE[i] - vE[i - 1]) /
										  (vy_xmax[l * M + i] - vy_xmax[l * M + i - 1]))
						  << " -> "
						  << E_after_secants
						  << " \t( " << E_true << " , " << E_after_secants / E_true - 1 << " err )"
						  << "   " << F(E_after_secants, l, a, b, N)
						  << std::endl;

				// Create an output of the wavefunction
				char filename[50];
				sprintf(filename,"data/eigen_k%ld_l%ld.dat",ii,l);
				std::ofstream file{filename};
				save_function(E_true, l, a, b, 3e6, file);
				file.close();
				ii++;
			}
		}
	}

	return 0;
}