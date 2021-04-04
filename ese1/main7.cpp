#include "uniconst.hpp"
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

double V_rot(const double r, const int l)
{
	return (l * (l + 1) / (r * r));
}

double u(const double r,const double B,const int l){
	return std::exp(-B * std::pow(r, -5.0));
}

int main()
{
	std::cout << "Ciao ese6!" << std::endl;

	// Physical quantities
	constexpr double eps_r = 5.9e-3; // eV
	constexpr double sigma_r = 3.18; // Ams
	constexpr double m_r =  (1. / (1 / (83.798) + 1 / (1.008)))*uni::Da_to_kg ; // kg

	// Reduced quantities
	constexpr double hbar2_2m = (uni::hbar_r * uni::hbar_r) / (2 * m_r * (eps_r * uni::eV_to_J) * (sigma_r * uni::Ams_to_m) * (sigma_r * uni::Ams_to_m));
	// const double B = 2. / 5. * std::sqrt(2 * m_r * (eps_r * uni::eV_to_J)) * (sigma_r * uni::Ams_to_m)/uni::hbar_r;
	constexpr double B = 2.1314279016703604310513399109722040174046218740519483004507;

	constexpr uint64_t N = 1e5; // number of mesh points
	constexpr double a = 0.8;	// start included
	constexpr double b = 15.;  // stop excluded
	// After this the particle is free, the wavefunction is a plane wave
	constexpr double rmax= 5;

	constexpr double h = (b - a) / N; // Mesh spacing
	constexpr double h2 = h * h;

	constexpr double r1 = 12;
	constexpr uint64_t ir1 = (uint64_t)((r1 - a) / h);
	constexpr uint64_t ir2 = N-1;
	constexpr double r2 = a + ir2 * h;
	static_assert(ir1 < N, "ir1 must be lower than N");
	static_assert(ir2 < N, "ir2 must be lower than N");

	constexpr uint64_t M = 300;
	constexpr double E_i = 0.1e-3 / eps_r;
	constexpr double E_f = 4e-3 / eps_r;
	constexpr double h_E = (E_f - E_i) / M;
	std::vector<double> vE(M);
	std::vector<double> sigma_tot(M);
	arange(vE,E_i,h_E);

	constexpr int lmax = 6;

	std::cout << "Probelm constants:"
			  << "\neps_r:   " << eps_r << " eV"
			  << "\nsigma_r: " << sigma_r << " Ams"
			  << "\nm_r:     " << m_r << " kg"
			  << "\nB:       " << B
			  << "\nhbar2_2m:" << hbar2_2m
			  << "\nN:       " << N
			  << "\na:       " << a
			  << "\nb:       " << b
			  << "\nM:       " << M
			  << "\nE_i:     " << E_i
			  << "\nE_f:     " << E_f
			  << "\nh_E:     " << h_E
			  << "\n"
			  << "\nrmax:    " << rmax
			  << "\nr1:      " << r1
			  << "\nr2:      " << r2
			  << "\n"
			  << std::endl;

	std::vector<double> dl(lmax + 1);
	std::vector<double> y(N);


	std::vector<double> x(N);
	for (uint64_t i = 0; i < N; i++)
	{
		x[i] = a + i * h;
	}

	// Loop over energies
	// for (uint64_t m = 0; m<M; m++)
	uint64_t m=M/2;
	{
		double E=vE[m];
		double k_out = std::sqrt(E / hbar2_2m);

		// Loop over l values
		// for (int l = 0; l <= lmax; l++)
		int l=5;
		{
			// Initial conditions
			y[0] = u(x[0],B,l);
			y[1] = u(x[1],B,l);

			// Numerov algorithm
			double k_ii2 = E - (V_LJ(x[0]) + hbar2_2m * V_rot(x[0],l));
			double k_i2 = E - (V_LJ(x[1]) + hbar2_2m * V_rot(x[1], l));
			double k2;
			// std::cout << y[0] << " " << y[1] << " " << std::pow(B / x[1], 5.0) << std::endl;

			// Loop over space
			for (uint64_t i = 2; i < N; i++)
			{

				// New k
				k2 = E - ((x[i] < rmax) ? V_LJ(x[i]) + hbar2_2m * V_rot(x[i], l) : hbar2_2m * V_rot(x[i], l));

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

			// std::cout << l << " " << tan_dl << " " << dl[l] << " " << std::endl;

			std::ofstream file{"data/Scatter.dat"};
			file << "#x y(l)\n";
			for (uint64_t i = 0; i < N; i++)
			{
				// Save the result
				file << x[i] << ' ' << y[i] << '\n';
			}
			file.close();

		}

		sigma_tot[m] = 0;
		for (int l = 0; l <= lmax; l++)
		{
			// 	std::cout << "l: " << l << "\tdl: " << dl[l] << std::endl;
			sigma_tot[m] += (2 * l + 1) * std::pow(std::sin(dl[l]), 2.0);
		}
		sigma_tot[m] *= 4 * M_PI / (k_out * k_out);
		
		if (m % 100 == 0)
		{
			std::cout << "m: " << m << " P: " << 1.0 * m / M << std::endl;
		}
	}

	

	// Save sigma_tot
	std::ofstream file{"data/sigma_tot.dat"};
	file << "#E sigma_tot\n";
	for (uint64_t m = 0; m < M; m++)
	{
		// Save the result
		file << vE[m] << ' ' << sigma_tot[m] << '\n';
	}
	file.close();

	return 0;
}