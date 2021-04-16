#include "vector_help.hpp"
#include "uniconst.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include <gsl/gsl_sf_bessel.h>
#include <static_math/cmath.h>

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
	constexpr double B = 2. / 5. * smath::sqrt(1./hbar2_2m);

	constexpr uint64_t N = 1e5; // number of mesh points
	constexpr double a = 0.8;	// start included

	// After this the particle is free, the wavefunction is a plane wave
	constexpr double rmax= 5.;
	constexpr double lam_wait_1=2.;
	constexpr double lam_wait_2=3./4;

	constexpr uint64_t M = 300;
	constexpr double E = 0.3;
	constexpr double lambda_max = 2.*M_PI/smath::sqrt(E/hbar2_2m);

	constexpr int lmax = 8;

	std::cout << "Probelm constants:"
			  << "\neps_r:      " << eps_r << " eV"
			  << "\nsigma_r:    " << sigma_r << " Ams"
			  << "\nm_r:        " << m_r << " kg"
			  << "\nB:          " << B
			  << "\nhbar2_2m:   " << hbar2_2m
			  << "\nN:          " << N
			  << "\na:          " << a
			  << "\nM:          " << M
			  << "\nE:          " << E
			  << "\nlmax:       " << lmax
			  << "\nrmax:       " << rmax
			  << "\nlam_wait_1: " << lam_wait_1
			  << "\nlam_wait_2: " << lam_wait_2
			  << "\nPotential values"
			  << "\nV_LJ(rmax):       " << V_LJ(rmax)
			  << "\nV_rot(rmax,1):    " << hbar2_2m * V_rot(rmax, 1)
			  << "\nV_rot(rmax,lmax): " << hbar2_2m * V_rot(rmax, lmax)
			  << "\n"
			  << std::endl;

	std::vector<double> dl(lmax + 1);
	std::vector<double> y(N);

	const double k_out = std::sqrt(E / hbar2_2m);
	const double lambda=2.*M_PI/k_out;

	// Redo the mesh
	const double r1=rmax+lam_wait_1*lambda;
	const double r2=r1+lam_wait_2*lambda;
	const double h = (r2 - a) / N; // Mesh spacing
	const double h2 = h * h;
	const uint64_t ir1 = (r1-a) / h; 
	const uint64_t ir2 = N-1;

	
	std::cout<<rmax<<' '<<r1<<' '<<r2<<std::endl;

	// Loop over l values
	for (int l = 0; l <= lmax; l++)
	// int l=0;
	{
		// Initial conditions
		y[0] = u(a+0*h,B,l);
		y[1] = u(a+1*h,B,l);

		// Numerov algorithm
		double k_ii2 = 1. / hbar2_2m * (E - (V_LJ(a + 0 * h) + hbar2_2m * V_rot(a + 0 * h, l)));
		double k_i2 = 1. / hbar2_2m * (E - (V_LJ(a + 1 * h) + hbar2_2m * V_rot(a + 1 * h, l)));
		double k2;
		// std::cout << y[0] << " " << y[1] << " " << std::pow(B / x[1], 5.0) << std::endl;

		// Loop over space (can stop at ir2)
		double x=a+2*h;
		for (uint64_t i = 2; i < N && i <= ir2 ; i++)
		{

			// New k
			k2 = 1. / hbar2_2m * 
				(E - (
					(x < rmax) ?
					V_LJ(x) + hbar2_2m * V_rot(x, l) : 
					hbar2_2m * V_rot(x, l))
					);

			y[i] = (2. * y[i - 1] * (1. - 5. / 12 * h2 * k_i2) - y[i - 2] *
						(1. + 1. / 12 * h2 * k_ii2)) /
					(1. + h2 / 12.0 * k2);

			// Move down the k for next iteration
			k_ii2 = k_i2;
			k_i2 = k2;
			x+=h;
		}

		// Phase shift
		double K = y[ir1] * r2 / (y[ir2] * r1);
		// std::cout << l<< k_out *r2 <<std::endl;
		double tan_dl = (K * gsl_sf_bessel_jl(l, k_out * r2) - gsl_sf_bessel_jl(l, k_out * r1)) /
						(K * gsl_sf_bessel_yl(l, k_out * r2) - gsl_sf_bessel_yl(l, k_out * r1));

		dl[l] = std::atan(tan_dl);

		// std::cout
		// 	<< l
		// 	<< "\t" << tan_dl
		// 	<< "\t" << dl[l]
		// 	<< "\t" << (2 * l + 1) * std::pow(std::sin(dl[l]), 2.0)
		// 	<< std::endl;


	}

	std::ofstream file{"data/dla0.8.dat"};
	file<<"#l dl (2l+1)sin^2(dl)\n";
	for (int l=0 ; l<=lmax;l++){
		file
			<< l << ' '
			<< dl[l] << ' '
			<< (2 * l + 1) * std::pow(std::sin(dl[l]), 2.0)
			<< '\n';
	}
	file.close();

	// std::ofstream file{"data/Emaxl0.dat"};
	// file << "#x y(0)\n";
	// for (uint64_t i = 0; i < N; i++)
	// {
	// 	// Save the result
	// 	file << a+i*h << ' ' << y[i] << '\n';
	// }
	// file.close();

	return 0;
}