#include "uniconst.hpp"
#include "vector_help.hpp"
#include "HF.hpp"
#include "MF.hpp"
#include "gsl_help.hpp"
#include "numerov.hpp"
#include "DFT.hpp"
#include "differentiate.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <static_math/cmath.h>

// Correlation energy parameters
constexpr double p = 1.0;
constexpr double A = 0.031091;
constexpr double alpha_1 = 0.21370;
constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};


// Correlation energy
double E_c(const double rs)
{
	const double DEN = 2 * A * (beta[0] * std::pow(rs, 0.5) + beta[1] * rs + beta[2] * std::pow(rs, 3.0 / 2) + beta[3] * std::pow(rs, p + 1));
	return -2 * A * (1 + alpha_1 * rs) * std::log(1 + 1.0 / DEN);
}

// Exchange energy
double E_x(const double rho)
{
	return -3. / 4 * std::pow(3 * rho / M_PI, 1. / 3);
}

// de_x/drho Derivative of exchange
double dE_x(const double rho)
{
	return -1. / 4 * std::pow(3. / M_PI, 1. / 3) * std::pow(rho, -2. / 3);
}

// de_c/drho Derivative of correlation
double dE_c(const double rho)
{
	const double rs = std::pow(3 / (4 * M_PI * rho), 1. / 3);
	const double dr_s = -1. / 3 * std::pow(3. / 4 / M_PI, 1. / 3) * std::pow(rho, -4. / 3);
	const double arg1 = 2 * A * (beta[0] * std::pow(rs, 0.5) + beta[1] * rs + beta[2] * std::pow(rs, 3. / 2) + beta[3] * std::pow(rs, p + 1));
	const double arg2 = 1 + 1. / arg1;
	return dr_s *(
		-2*A*alpha_1*rs*std::log(arg2)
		-2*A*(1+alpha_1*rs)*1./(arg2)*(-(2*A*
			(beta[0]*0.5/std::sqrt(rs)+beta[1]+beta[2]*3./2*std::sqrt(rs)+
			beta[3]*(p+1)*std::pow(rs,p)))/(arg1*arg1))
	);
}

// External potential
double V_ext(const double r, const double rhoB, const double Rc)
{
	double v_ext = 2 * M_PI * rhoB * 
	((r > Rc) ? 
		(-2. / 3 * std::pow(Rc, 3.0) / r): 
		1. / 3 * r * r - Rc * Rc);
	
	return v_ext;
}

double V_eff_0(const double r, const double rhoB, const double Rc, const int l)
{
	return V_ext(r, rhoB, Rc) + 0.5 * l * (l + 1) / (r * r);
}


int main()
{
	std::cout << "DFT: Independent Electron Model" << std::endl;
	
	// Problem constants
	constexpr double r_s_Na = 3.93;
	constexpr double r_s_K = 4.86;
	constexpr double rs = r_s_K;
	constexpr uint32_t N = 8;
	constexpr int Nlevels = 2; // N = 8
	// constexpr int Nlevels = 4; // N = 20
	// constexpr int Nlevels = 6; // N = 40

	// Angular momentum
	constexpr int lmax = 1;
	// Occupation order
	constexpr int is[Nlevels] = {0, 1}; // N = 8 lmax=1
	// constexpr int is[Nlevels] = {0, 3, 1, 2}; // N = 20 lmax=2
	// constexpr int is[Nlevels] = {0, 3, 1, 5, 2, 4}; // N = 40 lmax=3
	constexpr int ls[10] = {0, 1, 2, 0, 3, 1};
	// Number of energies for each l
	constexpr int Enums[lmax + 1] = {1, 1}; // N = 8
	// constexpr int Enums[lmax + 1] = {2, 1, 1}; // N = 20
	// constexpr int Enums[lmax + 1] = {2, 2, 1, 1}; // N = 40

	// Background properties
	constexpr double rhoB = 1 / (4. / 3 * M_PI * rs * rs * rs);
	const double Rc = std::pow(N * 3.0 / (rhoB * 4. * M_PI), 1. / 3);

	// Spatial mesh
	constexpr uint64_t M = 1e4;
	constexpr double a = 1e-5;
	constexpr double b = 20;
	constexpr double h = (b - a) / M;

	// Potential
	double V[M];

	// Eigenvectors, eigenvalues and density
	double Y[M*Nlevels];
	double Yxx[M*Nlevels];
	double rho[M];
	double Es[Nlevels];
	double y0,y1;
	
	// Energy mesh
	double Eh,Ea;
	constexpr double Eb = 0;
	constexpr uint64_t ME = 1e3;

	std::vector<double> ymax(ME);

	std::cout << "Problem constants:"
			  << "\nN: " << N
			  << "\nM: " << M
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  //   << "\nEa: " << Ea
			  << "\nEb: " << Eb
			//   << "\nEh: " << Eh
			  << "\nME: " << ME
			  << "\nrs: " << rs
			  << "\nrhoB: " << rhoB
			  << "\nRc: " << Rc
			//   << "\nmin(V): " << - 2 * M_PI * rhoB * Rc * Rc
			  << "\n"
			  << std::endl;

	int itot = 0,i;
	double E;
	for (int l = 0; l <= lmax; l++)
	{
		// Effective potential
		for (uint64_t m = 0; m < M; m++)
		{
			V[m] = V_eff_0(a + m * h, rhoB, Rc, l);
		}

		E=Ea = min(V, M);
		Eh = (Eb - Ea) / ME;
		y0 = std::pow(a, l + 1);
		y1 = std::pow(a + h, l + 1);

		// Solve the differential eq with numerov for each state
		for (uint8_t g = 0; g < Enums[l]; g++)
		{
			// Actual occupation order
			i=is[itot];

			E=Es[i] = numerov_find_energy(y0, y1, h, M, V, E+1e-10, Eb, Eh);
			std::cout << l << " " << Es[i] << std::endl;

			// for (uint64_t m = 0; m < ME; m++)
			// {
			// 	ymax[m] = numerov_integrate_yxmax(y0, y1, h, M, V, Ea + m * Eh);
			// }

			// Compute the WF
			(Y + i * M)[0] = y0;
			(Y + i * M)[1] = y1;
			numerov_integrate(Y + i * M, h, M, V, Es[i]);
			// Normalize
			double norm=0;
			for (uint64_t m = 0; m < M; m++)
			{
				norm += h * (Y + i * M)[m]*(Y + i * M)[m];
			}
			norm=std::sqrt(norm);
			for (uint64_t m = 0; m < M; m++)
			{
				(Y + i * M)[m]/=norm;
			}
			// Compute derivative
			diff_2_5points_allmesh((Y + i * M), M, (Yxx + i * M), h);

			itot++;
		}
	}

	// Compute rho
	for (uint64_t m = 0; m < M; m++)
	{
		// Construct rho
		double tmp = 0;
		itot=0;
		for (int l = 0; l <= lmax; l++)
		{
			for (int g = 0; g < Enums[l]; g++)
			{
				i=is[itot];
				tmp += 2 * (2 * l + 1) * Y[i * M + m] * Y[i * M + m];
				itot++;
			}
		}
		rho[m] = tmp;
	}
	
	// Compute the two differents values of energy to check the convergence
	double Eeigen = DFT_Eeigen(
		rho, M, a, h, Nlevels, ls, Es,
		[](double) { return 0; }, [](double) { return 0; }, [](double) { return 0; }, [](double) { return 0; }, [](double, double, uint64_t, const double *) { return 0; });


	double EMF = DFT_Efunc(
		rho, Y, Yxx, M, a, h, Nlevels, ls, rhoB, Rc,
		[](double) { return 0; }, [](double) { return 0; }, [](double, double, uint64_t, const double *) { return 0; }, V_ext);

	std::cout << Eeigen << " " << EMF << std::endl;

	{
	std::ofstream file{"data/numerov.dat"};
	file << "r phi rho \n";
	for (uint64_t m = 0; m < M; m++)
	{
		file << a + m * h;
		// for (int i = 0; i <= 3; i++)
		// {
		// 	rho[i*M+m] = (i==0?0:rho[(i-1)*M+m])+ Y[i*M+m]*Y[i*M+m];
		// 	file << ' ' << Y[i*M+m] << ' '<< rho[i*M+m]; // degeneracy
		// }
		file<<' '<<rho[m];
		file << '\n';
	}
	file.close();
	}

	// {
	// std::ofstream file{"data/y_max.dat"};
	// file << "E y_max \n";
	// for (uint64_t m = 0; m < ME; m++)
	// {
	// 	file << vE[m] << ' ' << ymax[m] << '\n';
	// }
	// file.close();
	// }

		
	
	
	return 0;
	
}