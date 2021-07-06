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
double V_ext( const double rhoB, const double r, const double Rc)
{
	double v_ext = 2 * M_PI * rhoB * 
	((r > Rc) ? 
		(-2. / 3 * std::pow(Rc, 3.0) / r): 
		1. / 3 * r * r - Rc * Rc);
	
	return v_ext;
}

double V_eff(const double rhoB, const double Rc, const double r, const int l)
{
	return V_ext(rhoB, Rc, r) + 0.5 * l * (l + 1) / (r * r);
}

// double V_eff(const double rhoB, const double Rc, const double r, const int l)
// {
// 	return 0.5 * r* r + 0.5 * l * (l + 1) / (r * r);
// }

int main()
{
	std::cout << "DFT: Independent Electron Model" << std::endl;
	
	// Problem constants
	constexpr double r_s_Na = 3.93;
	constexpr double r_s_K = 4.86;
	constexpr double rs = r_s_K;
	constexpr uint32_t Ne = 8;
	constexpr uint8_t Nlevels = 4;

	// Initial guess
	constexpr double rhoB = 1 / (4. / 3 * M_PI * rs * rs * rs);
	// const uint8_t lmax=2;
	const double Rc = std::pow(Ne * 3.0 / (rhoB * 4. * M_PI), 1. / 3);

	// Harmonic comparison
	// constexpr double k = 4 * M_PI * rhoB / 3;
	// constexpr double ke = smath::sqrt(k);
	

	// Spatial mesh
	constexpr uint64_t M = 1e4;
	constexpr double a = 1e-5;
	constexpr double b = 20;
	constexpr double h = (b - a) / M;
	// std::vector<double> r(M);
	// arange(r, a, h);

	double Y[M*Nlevels];
	double y0,y1;
	double rho[M*Nlevels];
	std::vector<double> E(Nlevels);
	constexpr int ls[]={0,1,2,0};
	// Starting E
	double Eas[Nlevels]={-1.2,-1.2,-0.8,-0.5};
	constexpr double Eb=3.5;
	constexpr double Eh=1e-2;
	constexpr uint64_t ME=(Eb+1.2)/Eh;
	std::vector<double> vE(ME);
	std::vector<double> ymax(ME);
	std::vector<double> y(M);
	arange(vE,-1.2,Eh);
	int l;
	std::function<double(double)> V;

	std::cout << "Problem constants:"
			  << "\n Ne: " << Ne
			  << "\nM: " << M
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  //   << "\nEa: " << Ea
			  << "\nEb: " << Eb
			  << "\nEh: " << Eh
			  << "\nME: " << ME
			  << "\nrs: " << rs
			  << "\nrhoB: " << rhoB
			  << "\nRc: " << Rc
			  << "\nkh: " << 4 * M_PI * rhoB / 3
			  << "\n"
			  << std::endl;

	for (int i = 0; i < Nlevels; i++)
	// for (int i = 3; i < Nlevels; i=100)
	{
		// Solve the differential eq with numerov for each state
		l = ls[i];
		V = [&](double r) -> double {
			return V_eff(rhoB, Rc, r, l);
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
		numerov_integrate(Y + i * M, a, h, M, V, E[i]);
		// numerov_integrate(Y + i * M, a, h, M, V, 0.5/100);
	}

	// {
	// std::ofstream file{"data/numerov.dat"};
	// file << "r phi rho \n";
	// for (uint64_t m = 0; m < M; m++)
	// {
	// 	file << a + m * h;
	// 	for (int i = 0; i <= 3; i++)
	// 	{
	// 		rho[i*M+m] = (i==0?0:rho[(i-1)*M+m])+ Y[i*M+m]*Y[i*M+m];
	// 		file << ' ' << Y[i*M+m] << ' '<< rho[i*M+m]; // degeneracy
	// 	}
	// 	file << '\n';
	// }
	// file.close();
	// }
	
	{
	std::ofstream file{"data/y_max.dat"};
	file << "E y_max \n";
	for (uint64_t m = 0; m < ME; m++)
	{
		file << vE[m] << ' ' << ymax[m] << '\n';
	}
	file.close();
	}

		
	
	// Compute the two differents values of energy to check the convergence 
	return 0;
	
}