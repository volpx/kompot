#include "uniconst.hpp"
#include "vector_help.hpp"
#include "HF.hpp"
#include "MF.hpp"
#include "gsl_help.hpp"
#include "numerov.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <static_math/cmath.h>

// #include <eigen3/Eigen/Dense>
// #define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

// Correlation energy parameters
constexpr double p = 1.0;
constexpr double A = 0.031091;
constexpr double alpha_1 = 0.21370;
constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};


// Correlation energy
double E_c(const double r_s)
{
	const double DEN = 2 * A * (beta[0] * std::pow(r_s, 0.5) + beta[1] * r_s + beta[2] * std::pow(r_s, 3.0 / 2) + beta[3] * std::pow(r_s, p + 1));
	return -2 * A * (1 + alpha_1 * r_s) * std::log(1 + 1.0 / DEN);
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

	const double r_s = std::pow(3 / (4 * M_PI * rho), 1. / 3);
	const double dr_s = -1. / 3 * std::pow(3. / 4 / M_PI, 1. / 3) * std::pow(rho, -4. / 3);
	const double arg1 = 2 * A * (beta[0] * std::pow(r_s, 0.5) + beta[1] * r_s + beta[2] * std::pow(r_s, 3. / 2) + beta[3] * std::pow(r_s, p + 1));
	const double arg2 = 1 + 1. / arg1;
	return dr_s *(
		-2*A*alpha_1*r_s*std::log(arg2)
		-2*A*(1+alpha_1*r_s)*1./(arg2)*(-(2*A*
			(beta[0]*0.5/std::sqrt(r_s)+beta[1]+beta[2]*3./2*std::sqrt(r_s)+
			beta[3]*(p+1)*std::pow(r_s,p)))/(arg1*arg1))
	);
}

// External potential
double V_ext( const double rhoB, const double r, const double R_c)
{
	double v_ext = 2 * M_PI * rhoB * 
	((r > R_c) ? 
		(-2. / 3 * std::pow(R_c, 3.0) / r + R_c * R_c): 
		1. / 3 * r * r);
	
	return v_ext;
}

double V_eff(const double rhoB, const double R_c, const double r, const int l)
{
	return V_ext(rhoB, R_c, r) + 0.5 * l * (l + 1) / (r * r);
}

// double V_eff(const double rhoB, const double R_c, const double r, const int l)
// {
// 	return 0.5 * r* r + 0.5 * l * (l + 1) / (r * r);
// }

int main()
{
	std::cout << "DFT: Independent Electron Model" << std::endl;
	
	// Problem constants
	constexpr double r_s_Na = 3.93;
	constexpr double r_s_K = 4.86;
	constexpr double r_s = r_s_K;
	constexpr uint32_t Ne = 19;
	constexpr uint8_t Nlevels = 4;

	// Initial guess
	constexpr double rho_0 = 1 / (4. / 3 * M_PI * r_s * r_s * r_s);
	
	// const uint8_t lmax=2;
	const double R_c = std::pow(Ne * 3.0 / (rho_0 * 4. / 3 * M_PI), 1. / 3);

	// Spatial mesh
	constexpr uint64_t M = 1e5;
	constexpr double a = 1e-5;
	constexpr double b = 3000;
	constexpr double h = (b - a) / M;
	std::vector<double> r(M);
	arange(r, a, h);

	double Y[M*Nlevels];
	double rho[M*Nlevels];
	std::vector<double> E(Nlevels);
	double y0,y1;
	constexpr int ls[]={0,1,2,0};
	constexpr double Eb=2;
	constexpr double Eh=1e-3;
	constexpr double Eas[]={0.1,0.1,0.1,0.224+Eh};
	constexpr uint64_t ME=(Eb)/Eh;
	std::vector<double> vE(ME);
	std::vector<double> ymax(ME);
	std::vector<double> y(M);
	arange(vE,0,Eh);
	int l;
	std::function<double(double)> V;

	for (int i = 0; i < Nlevels; i++)
	// for (int i = 0; i < Nlevels; i=100)
	{
		// Solve the differential eq with numerov for each state
		l = ls[i];
		V = [&](double r) -> double {
			return V_eff(rho_0, R_c, r, l);
			// return 0;
		};
		
		y0 = 1e-2*std::pow(a, l + 1);
		y1 = 1e-2*std::pow(a + h, l + 1);
		E[i] = numerov_find_energy(y0, y1, a, h, M, V, Eas[i], Eb, Eh); 
		std::cout << E[i] << std::endl;
		for (uint64_t m = 0; m < ME; m++)
		{
			ymax[m] = numerov_integrate_yxmax(y0, y1, a, h, M, V, vE[m]);
		}

		// Compute the WF
		(Y + i * M)[0] = y0;
		(Y + i * M)[1] = y1;
		numerov_integrate(Y + i * M, a, h, M, V, E[i]);//degenerazione
		// numerov_integrate(Y + i * M, a, h, M, V, 0.5/100);
	}

	{
	std::ofstream file{"data/numerov.dat"};
	file << "r phi rho \n";
	for (uint64_t m = 0; m < ME; m++)
	{
		file << a + m * h;
		for (int i = 0; i <= 3; i++)
		{
			rho[i*M+m] = (i==0?0:rho[(i-1)*M+m])+ Y[i*M+m]*Y[i*M+m];
			file << ' ' << Y[i*M+m] << ' '<< rho[i*M+m];
		}
		file << '\n';
	}
	file.close();
	}
	
	{
	std::ofstream file{"data/y_max.dat"};
	file << "E y_max \n";
	for (uint64_t m = 0; m < ME; m++)
	{
		file << 0 + m * Eh << ' ' << ymax[m] << '\n';
	}
	file.close();
	}

		
	
	// Compute the two differents values of energy to check the convergence 
	return 0;
	
}