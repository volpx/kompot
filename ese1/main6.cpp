#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "gsl/gsl_sf_bessel.h"

double V_LJ(const double r)
{
	return (4.0 * (std::pow(r,-12.0)-std::pow(r,-6.0)));
}

double V_eff(const double r, const int l)
{
	return (V_LJ(r) + 0.5 * l * (l + 1) / (r * r));
}


int main()
{
	std::cout << "Ciao ese6!" << std::endl;

	constexpr double eps = 5.9e-3; // eV
	constexpr double sigma = 3.18; // Ams

	constexpr double B = 2.0 / 5;

	constexpr uint64_t N = 1e4; // number of mesh points
	constexpr double a = 0.5;	// start included
	constexpr double b = 50.;	// stop excluded
	constexpr double rmax= 5.0;

	constexpr double h = (b - a) / N; // Mesh spacing
	constexpr double h2 = h * h;

	constexpr double r1 = 25;
	constexpr uint64_t ir1 = (uint64_t)((r1 - a) / h);
	constexpr double r2 = b - 0.1;
	constexpr uint64_t ir2 = (uint64_t)((r2 - a) / h);
	static_assert(ir1 < N, "ir1 must be lower than N");
	static_assert(ir2 < N, "ir2 must be lower than N");

	constexpr uint64_t M = 1e3;
	constexpr double E_i = 0.1e-3 / eps;
	constexpr double E_f = 3.5e-3 / eps;
	constexpr double h_E = (E_f - E_i) / M;
	std::vector<double> vE(M);
	std::vector<double> sigma_tot(M);
	arange(vE,E_i,h_E);

	std::cout << "Probelm constants:"
			  << "\neps:\t" << eps
			  << "\nsigma:\t" << sigma
			  << "\nN:\t" << N
			  << "\na:\t" << a
			  << "\nb:\t" << b
			  << "\nM:\t" << M
			  << "\nE_i:\t" << E_i
			  << "\nE_f:\t" << E_f
			  << "\nh_E:\t" << h_E
			  << "\nr1:\t" << r1
			  << "\nir1:\t" << ir1
			  << "\nr2:\t" << r2
			  << "\nir2:\t" << ir2
			  << std::endl;

	constexpr int lmax = 6;
	std::vector<double> dl(lmax + 1);
	std::vector<double> y(N);


	std::vector<double> x(N);
	for (uint64_t i = 0; i < N; i++)
	{
		x[i] = a + i * h;
	}


	double E=0.3;
	double k_out = std::sqrt(E);

	for (int l = 0; l <= lmax; l++)
	{
		// Initial conditions
		y[0] = std::exp(-std::pow(B / x[0], 5.0));
		y[1] = std::exp(-std::pow(B / x[1], 5.0));

		// Numerov algorithm
		double k_ii2 =(E - V_eff(x[0], l));
		double k_i2 = (E - V_eff(x[1], l));
		double k2;
		// std::cout << y[0] << " " << y[1] << " " << std::pow(B / x[1], 5.0) << std::endl;

		for (uint64_t i = 2; i < N; i++)
		{

			// New k
			k2 =   (E - V_eff(x[i], l));

			y[i] = (2. * y[i - 1] * (1. - 5. / 12 * h2 * k_i2) - y[i - 2] *
						(1. + 1. / 12 * h2 * k_ii2)) /
					(1. + h2 / 12.0 * k2);

			// Move down the k for next iteration
			k_ii2 = k_i2;
			k_i2 = k2;

		}

		// Phase shift
		double K = y[ir1] * r2 / (y[ir2] * r1);
		// std::cout << l<< k_out *r2 <<std::endl;
		double tan_dl = (K * gsl_sf_bessel_jl(l, k_out * r2) - gsl_sf_bessel_jl(l, k_out * r1)) /
						(K * gsl_sf_bessel_yl(l, k_out * r2) - gsl_sf_bessel_yl(l, k_out * r1));

		dl[l] = std::atan(tan_dl);
		// std::cout << l << " " << tan_dl << " " << dl[l] << " " << y[(uint64_t)(r2 / h)] << std::endl;
	}

	// sigma_tot[m] = 0;
	// for (int l = 0; l <= lmax; l++)
	// {
	// 	sigma_tot[m] += (2 * l + 1) * std::pow(std::sin(dl[l]), 2.0);
	// }
	// sigma_tot[m] *= 4 * M_PI / (k_out * k_out);
	
	// if (m % 100 == 0)
	// {
	// 	std::cout << "m: " << m << " P: " << 1.0 * m / M << std::endl;
	// }
	

	std::ofstream file{"data/Scatter.dat"};
	file << "#x y(l)\n";
	for (uint64_t i = 0; i < N; i++)
	{
		// Save the result
		file << x[i] << ' ' << y[i] << '\n';
	}
	file.close();

	// for (int l=0;l<=lmax;l++)
	// 	std::cout << "l: " << l << "\tdl: " << dl[l] << std::endl;

	// Save sigma_tot
	// std::ofstream file{"data/sigma_tot.dat"};
	// file << "#E sigma_tot\n";
	// for (uint64_t m = 0; m < M; m++)
	// {
	// 	// Save the result
	// 	file << vE[m] << ' ' << sigma_tot[m] << '\n';
	// }
	// file.close();

	return 0;
}