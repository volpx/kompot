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

// #include <static_math/cmath.h>

// External potential
double V_ext(const double r, const double rhoB, const double Rc)
{
	double v_ext = 2 * M_PI * rhoB * 
		((r <= Rc) ? 
			1. / 3 * r * r - Rc * Rc:
			(-2. / 3 * std::pow(Rc, 3.0) / r));
	
	return v_ext;
}

double V_eff(const double r, const double rhoB, const double Rc, const int l)
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

	#if 1
	constexpr uint32_t N = 8;
	constexpr uint32_t lmax = 1;
	constexpr int Nlevels = 2; // N = 8
	constexpr int is[Nlevels] = {0, 1}; // N = 8 lmax=1
	constexpr int Enums[lmax + 1] = {1, 1}; // N = 8
	#endif

	#if 0 
	constexpr uint32_t N = 20;
	constexpr int lmax = 2;
	constexpr int Nlevels = 4; // N = 20
	constexpr int is[Nlevels] = {0, 3, 1, 2}; // N = 20 lmax=2
	constexpr int Enums[lmax + 1] = {2, 1, 1}; // N = 20
	#endif

	#if 0
	constexpr uint32_t N = 40;
	constexpr int lmax = 3;
	constexpr int Nlevels = 6; // N = 40
	constexpr int is[Nlevels] = {0, 3, 1, 5, 2, 4}; // N = 40 lmax=3
	constexpr int Enums[lmax + 1] = {2, 2, 1, 1}; // N = 40
	#endif

	constexpr int ls[10] = {0, 1, 2, 0, 3, 1};

	// constexpr uint32_t N = 8;
	// constexpr int Nlevels = 2; // N = 8
	// // constexpr int Nlevels = 4; // N = 20
	// // constexpr int Nlevels = 6; // N = 40

	// // Angular momentum
	// constexpr int lmax = 1;
	// // Occupation order
	// constexpr int is[Nlevels] = {0, 1}; // N = 8 lmax=1
	// // constexpr int is[Nlevels] = {0, 3, 1, 2}; // N = 20 lmax=2
	// // constexpr int is[Nlevels] = {0, 3, 1, 5, 2, 4}; // N = 40 lmax=3
	// // Number of energies for each l
	// constexpr int Enums[lmax + 1] = {1, 1}; // N = 8
	// // constexpr int Enums[lmax + 1] = {2, 1, 1}; // N = 20
	// // constexpr int Enums[lmax + 1] = {2, 2, 1, 1}; // N = 40

	// Background properties
	constexpr double rhoB = 1 / (4. / 3 * M_PI * rs * rs * rs);
	const double Rc = std::pow(N * 3.0 / (rhoB * 4. * M_PI), 1. / 3);

	// Spatial mesh
	constexpr double a = 1e-4;
	constexpr double b = 30;
	constexpr double h = 1e-3;
	constexpr size_t M = static_cast<size_t>((b-a)/h);
	size_t m_explode=M;

	// Potential
	double V[M];
	double rhor2[M];
	double rhor[M];

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
			  << std::endl;

	// Variables 
	int itot = 0;
	int i;
	double E;
	double Eeigen = 1.0;
	double Efunc = 0;
	double Ymin, norm;
	fill(rho, M, 0);

	for (uint32_t l = 0; l <= lmax; l++)
	{
		// Effective potential
		for (uint64_t m = 0; m < M; m++)
		{
			V[m] = V_eff(a + m * h, rhoB, Rc, l);
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

			// Move up a bit from the previous value/minimum
			Ea=E+1e-10;
			std::cout 
				<< l << g 
				<< " Ea: "<<Ea<<std::endl;
			E = Es[i] = numerov_find_energy(y0, y1, h, M, V, Ea, Eb, Eh);

			// Compute the WF
			(Y + i * M)[0] = y0;
			(Y + i * M)[1] = y1;
			numerov_integrate(Y + i * M, h, M, V, Es[i]);

			// Search where it explodes
			Ymin=std::abs((Y + i * M)[M-1]);
			for (uint64_t m=0; m < M; m++)
			{
				if(Ymin >= std::abs((Y + i * M)[M-1-m])){
					Ymin=std::abs((Y + i * M)[M-1-m]);
				}
				else
				{
					m_explode=M-m+1;
					break;
				}
			}

			// Normalize
			norm=0;
			for (uint64_t m = 0; m < M; m++)
			{
				if(m>=m_explode)
				{
					(Y + i * M)[m] = (Y + i * M)[m_explode-1] 
						* std::exp(-(h/2*(m-m_explode)));
				}
				norm += h * (Y + i * M)[m] * (Y + i * M)[m];
			}
			norm = std::sqrt(norm);
			for (uint64_t m = 0; m < M; m++)
			{
				(Y + i * M)[m] /= norm*std::sqrt(4*M_PI);
			}

			// Normalize
			// norm=0;
			// for (uint64_t m = 0; m < M; m++)
			// {
			// 	norm += h * (Y + i * M)[m]*(Y + i * M)[m];
			// }
			// norm=std::sqrt(norm);

			// for (uint64_t m = 0; m < M; m++)
			// {
			// 	(Y + i * M)[m]/=norm;
			// }

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
		for (uint32_t l = 0; l <= lmax; l++)
		{
			for (int g = 0; g < Enums[l]; g++)
			{
				i=is[itot];
				tmp += 2 * (2 * l + 1) * Y[i * M + m] * Y[i * M + m] / ((a + m * h) * (a + m * h)); ;
				itot++;
			}
		}
		rho[m] = tmp;
		rhor[m] = 0;//rho[m]*(a+m*h);
		rhor2[m] = 0;//rhor[m]*(a+m*h);
	}
	
	// Compute the two differents values of energy to check the convergence

	Eeigen = DFT_Eeigen(rho,rhor,rhor2, M, a, h, Nlevels, ls, Es,
			 [](double) { return 0; }, [](double) { return 0; }, 
			[](double) { return 0; }, [](double,double,double,size_t,const double *,const double *) { return 0; });

	Efunc = DFT_Efunc(rho,rhor,rhor2, Y, Yxx, M, a, h, Nlevels, ls, rhoB, Rc, 
		[](double) { return 0; }, [](double) { return 0; }, [](double,double,double,size_t,const double *,const double *) { return 0; }, V_ext);

	// double Eeigen = DFT_Eeigen(
	// 	rho, M, a, h, Nlevels, ls, Es,
	// 	[](double) { return 0; }, [](double) { return 0; }, [](double) { return 0; }, [](double) { return 0; }, [](double, double, uint64_t, const double *) { return 0; });


	// double EMF = DFT_Efunc(
	// 	rho, Y, Yxx, M, a, h, Nlevels, ls, rhoB, Rc,
	// 	[](double) { return 0; }, [](double) { return 0; }, [](double, double, uint64_t, const double *) { return 0; }, V_ext);

	std::cout << Eeigen << " " << Efunc << std::endl;

	{
	std::ofstream file{"data/numerov.dat"};
	file << Eeigen << ' ' << Efunc << '\n';
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