#include "uniconst.hpp"
#include "vector_help.hpp"
#include "HF.hpp"
#include "MF.hpp"
#include "gsl_help.hpp"
#include "numerov.hpp"
#include "DFT.hpp"
#include "differentiate.hpp"
#include "integrate.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <static_math/cmath.h>

// Exchange energy
// e_x(rho(r)) (scalar r)
double e_x(const double rho)
{
	return -3. / 4 * std::pow(3 * rho / M_PI, 1. / 3);
}

// de_x/drho Derivative of exchange
double de_x(const double rho)
{
	return -1. / 4 * std::pow(3. / M_PI, 1. / 3) * std::pow(rho, -2. / 3);
}

// Exchange potential 
double V_e_x(const double rho, std::function<double(double)> e_x){
	return 4. / 3 * e_x(rho);
}

// Correlation energy
// e_c(rho(r)) (scalar r)
double e_c(const double rho)
{
	// Correlation energy parameters
	constexpr double A = 0.031091;
	constexpr double alpha_1 = 0.21370;
	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};

	const double rs = std::pow(3. / (4 * M_PI * rho), 1. / 3);
	const double DEN = 2 * A * (beta[0] * std::pow(rs, 0.5) + beta[1] * rs + beta[2] * std::pow(rs, 3.0 / 2) + beta[3] * rs * rs);
	return -2 * A * (1 + alpha_1 * rs) * std::log(1 + 1.0 / DEN);
}

// de_c/drho Derivative of correlation
double de_c(const double rho)
{

	// Correlation energy parameters
	constexpr double A = 0.031091;
	constexpr double alpha_1 = 0.21370;
	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};

	// 
	const double c1 = std::pow(3. / (4 * M_PI), 1. / 3);
	const double rs = c1 * std::pow(rho,-1./3);
	const double c2 = 2 * A * alpha_1 * c1;
	const double den = 2 * A * (beta[0] * std::sqrt(rs) + beta[1] * rs + beta[2] * std::pow(rs, 3. / 2) + beta[3] * rs * rs);
	const double c3 = 1 + 1./den;

	return + c2 / (3 * std::pow(rho, 4./3)) * std::log(c3)
		   - (2 * A + c2 / std::pow(rho, 1./3)) / (den * den * c3) * 2 * A
		   * (+ beta[0] * std::sqrt(c1) / (6 * std::pow(rho, 7./6))
			  + beta[1] * c1 / (3 * std::pow(rho, 4./3))
			  + 0.5 * beta[2] * std::pow(c1 / rho,3./2)		
			  + 2. / 3 * beta[3] * c1 * c1 / std::pow(rho,5./3)
			);
}

// Correlation potential
double V_e_c(const double rho, std::function<double(double)> e_c, std::function<double(double)> de_c)
{
	return e_c(rho) + de_c(rho) * rho;
}

// External potential
double V_ext(const double r, const double rhoB, const double Rc)
{
	double v_ext = 2 * M_PI * rhoB * 
	((r < Rc) ? 
		1. / 3 * r * r - Rc * Rc:
		(-2. / 3 * std::pow(Rc, 3.0) / r));
	
	return v_ext;
}

// Direct term
double U_r(const double r, const double a, const double h,const uint64_t M, const double rho[])
{
	double u = 0;
	double rp = 0;

	// Direct term integral
	uint64_t m=0;
	while(m < M && rp<=r)
	{
		rp = a + m * h;
		u += h * (rp * rp * rho[m]);
		m++;
	}
	u=u/r;
	while(m < M)
	{
		rp = a + m * h;
		u += h * rp * rho[m];
		m++;
	}
	return 4 * M_PI * u;
}

double V_rot(
	const double r, const int l)
{
	return 0.5 * l * (l + 1) / (r * r);
}

// V_eff(r) compute V_eff at distance r (with index i)
double V_eff(
	const double r, const uint64_t i,
	const double rho[], const uint64_t M, 
	const double h, const double a,
	const double Rc,const double rhoB)
{
	// Vext + Vrot + Vdirect + Exc + dExc/drho * rho
	return V_ext(r ,rhoB, Rc)
		// + U_r(r, a, h, M, rho) 
		+ 4. / 3 * e_x(rho[i]) + (e_c(rho[i]) + de_c(rho[i]) * rho[i]);
}






// // Potential definition

// // External potential
// double V_ext(const double r, const double rhoB, const double Rc)
// {
// 	double v_ext = 2 * M_PI * rhoB * 
// 	((r < Rc) ? 
// 		(1. / 3 * r * r - Rc * Rc):
// 		(-2. / 3 * std::pow(Rc, 3.0) / r));
	
// 	return v_ext;
// }

// // Centrifugal potential
// double V_rot(const double r, const int l)
// {
// 	return 0.5 * l * (l + 1) / (r * r);
// }

// // Direct potential
// double U_r(const double r, const double a, const double h,const uint64_t M, const double rho[])
// {
// 	double u = 0;
// 	double rp = a;

// 	// Integral's computation
// 	uint64_t m=0;
// 	while(m < M && rp<=r)
// 	{
// 		rp = a + m * h;
// 		u += h * (rp * rp * rho[m]);
// 		m++;
// 	}
// 	u=u/r;
// 	while(m < M)
// 	{
// 		rp = a + m * h;
// 		u += h * rp * rho[m];
// 		m++;
// 	}
// 	return 4 * M_PI * u;
// }

// // Exchange energy
// double e_x(const double rho)
// {
// 	return -3. / 4 * std::pow(3 * rho / M_PI, 1. / 3);
// }

// // Correlation energy
// double e_c(const double rho)
// {
// 	// Correlation energy parameters
// 	constexpr double A = 0.031091;
// 	constexpr double alpha_1 = 0.21370;
// 	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};

// 	const double rs = std::pow(3. / (4 * M_PI * rho), 1. / 3);
// 	const double DEN = 2 * A * (beta[0] * std::pow(rs, 0.5) + beta[1] * rs + beta[2] * std::pow(rs, 3.0 / 2) + beta[3] * rs * rs); // beta_3 r_s^(p+1) but p = 1
// 	return -2 * A * (1 + alpha_1 * rs) * std::log(1 + 1.0 / DEN);
// }

// // Derivative of the correlation energy
// double de_c(const double rho)
// {

// 	// Correlation energy parameters
// 	constexpr double A = 0.031091;
// 	constexpr double alpha_1 = 0.21370;
// 	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};

// 	// Definition of the correlation energy
// 	const double c1 = std::pow(3. / (4 * M_PI), 1. / 3);
// 	const double rs = c1 * std::pow(rho, -1. / 3);
// 	const double c2 = 2 * A * alpha_1 * c1;
// 	const double den = 2 * A * (beta[0] * std::sqrt(rs) + beta[1] * rs + beta[2] * std::pow(rs, 3. / 2) + beta[3] * rs * rs); // beta_3 r_s^(p+1) but p = 1
// 	const double c3 = 1 + 1./den;

// 	return + c2 / (3 * std::pow(rho, 4. / 3)) * std::log(c3)
// 		   - (2 * A + c2 / std::pow(rho, 1. / 3)) / (den * den * c3) * 2 * A
// 		   * (+ beta[0] * std::sqrt(c1) / (6 * std::pow(rho, 7. / 6))
// 			  + beta[1] * c1 / (3 * std::pow(rho, 4. / 3))
// 			  + 0.5 * beta[2] * std::pow(c1 / rho, 3. / 2)		
// 			  + 2 * beta[3] * c1 * c1 / (3 * std::pow(rho, 5. / 3))
// 			);
// }

// // Exchange correlation potential
// double V_xc(const double rho)
// {
// 	return e_x(rho) / 0.75 + e_c(rho) + de_c(rho) * rho;
// }

// // KS potential
// double V_KS(const double r, const uint64_t i,
// 	const double rho[], const uint64_t M, 
// 	const double h, const double a,
// 	const double Rc,const double rhoB)
// {
// 	return V_xc(rho[i]) + V_ext(r, rhoB, Rc) + U_r(r, a, h, M, rho) ;
// }



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
	constexpr double rhoB = 1. / (4. / 3 * M_PI * rs * rs * rs);
	const double Rc = std::pow(N, 1. / 3) * rs;

	// Spatial mesh
	constexpr uint64_t M = 3e4;
	constexpr double a = 1e-4;
	constexpr double b = 20;
	constexpr double h = (b - a) / M;
	uint64_t m_explode=M;
	uint64_t m_explode_min=M;
	uint64_t m_explode_min_before=M; 

	// Potential
	double V[M];
	double V_0[M];
	double V_c[M];
	double V_x[M];
	double U_R[M];
	double V_Ext[M];

	// Eigenvectors, eigenvalues and density
	double Y[M*Nlevels];
	double Yxx[M*Nlevels];
	double rho[M];
	double rho_new[M];
	double Es[Nlevels];
	double y0,y1;
	
	// Energy mesh
	double Eh,Ea;
	constexpr double Eb = 10;
	constexpr uint64_t ME = 1e3;

	// Mixing parameter
	constexpr double alpha = 1e-2;
	// Accuracy
	constexpr double epsilon = 1e-6;

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
			  << "\n"
			  << std::endl;


	int itot = 0;
	int i;
	double E;
	double Eeigen = 1.0;
	double Efunc = 0;
	for (uint64_t m = 0; m < M; m++)
	{
		rho[m]=0;
	}
	uint64_t step = 0;
	while(std::abs(Eeigen - Efunc) > epsilon && step <=30)
	{
		std::cout<<"Step: "<<step<<std::endl;
		std::cout<<"m_explode_min_before: "<<m_explode_min_before<<std::endl;
		itot=0;

		for (uint64_t m = 0; m < M; m++)
		{
			rho_new[m]=0;
		}
		m_explode_min=M;

		// Compute potential without the l dependancy
		for (uint64_t m = 0; m < M; m++)
		{
			if(step == 0)
				V_0[m] = V_ext(a + m * h, rhoB, Rc);
			else
			{
				// V_c[m] = V_e_c(rho[m], e_c, de_c);
				// V_x[m] = V_e_x(rho[m], e_x);
				// U_R[m] = U_r(a + h * m, a, h, M, rho);
				// V_Ext[m] = V_ext(a + m * h, rhoB, Rc);
				if (m < m_explode_min_before)
				{
					V_0[m] = V_eff(a + m * h, m, rho, m_explode_min_before, h, a, Rc, rhoB);
					// V_0[m] = V_KS(a + m * h, m, rho, m_explode_min_before, h, a, Rc, rhoB);
				}
				else
				{
					V_0[m] = V_0[m_explode_min_before - 1] * std::exp(h * (1.0 * m_explode_min_before - m - 1) / 5);
					// V_l[m] = 0; //V_eff_0(a + m * h, m, rho, m_explode_min_before, h, Rc, rhoB, l);
				}
			}
		}
		
		for (int l = 0; l <= lmax; l++)
		{
			// Solve for value l
			for (uint64_t m = 0; m < M; m++)
			{
				V[m] = V_0[m] + V_rot(a + m * h, l);
			}

			E=Ea = min(V, m_explode_min_before);
			Eh = (Eb - Ea) / ME;
			y0 = std::pow(a, l + 1);
			y1 = std::pow(a + h, l + 1);

			// Solve the differential eq with numerov for each state
			for (int g = 0; g < Enums[l]; g++)
			{
				// Actual occupation order
				i=is[itot];

				Ea=E+1e-10;
				E = Es[i] = numerov_find_energy(y0, y1, h, M, V, Ea, Eb, Eh);
				std::cout << '\t' << l << g << " Ea: "<<Ea<< " E: " << Es[i] << std::endl;
				if(std::isnan(E))
				{
					std::cout<<"is nan"<<'\n';
					for (uint64_t m = 0; m < ME; m++)
					{
						ymax[m] = numerov_integrate_yxmax(y0, y1, h, M, V, Ea + m * Eh);
					}

					{
						std::ofstream file{"data/y_max.dat"};
						file << "E y_max \n";
						for (uint64_t m = 0; m < ME; m++)
						{
							file << Ea+m*Eh << ' ' << ymax[m] << '\n';
						}
						file.close();
					}

					numerov_integrate(Y + i * M, h, M, V, Ea);

					{
						std::ofstream file{"data/numerov_sc_y.dat"};
						file << "r phi rho V_l\n";
						for (uint64_t m = 0; m < M; m++)
						{
							file << a + m * h;
							// for (int i = 0; i <= 3; i++)
							// {
							// 	rho[i*M+m] = (i==0?0:rho[(i-1)*M+m])+ Y[i*M+m]*Y[i*M+m];
							// 	file << ' ' << Y[i*M+m] << ' '<< rho[i*M+m]; // degeneracy
							// }
							file<<' '<<(Y + i * M)[m]<<' '<<V_0[m];
							file << '\n';
						}
						file.close();
					}

					return 0;
				}

				// Compute the WF
				(Y + i * M)[0] = y0;
				(Y + i * M)[1] = y1;
				numerov_integrate(Y + i * M, h, M, V, Es[i]);

				// Search where it explodes
				double val=std::abs((Y + i * M)[M-1]);
				for (uint64_t m=0; m < M; m++)
				{
					if(val >= std::abs((Y + i * M)[M-1-m])){
						val=std::abs((Y + i * M)[M-1-m]);
					}
					else
					{
						m_explode=M-m+1;
						break;
					}
				}
				// std::cout<<m_explode<<std::endl;
				m_explode_min = std::min(m_explode, m_explode_min);

				// Normalize
				double norm=0;
				for (uint64_t m = 0; m < M; m++)
				{
					if(m<m_explode)
						norm += h * (Y + i * M)[m] * (Y + i * M)[m];
					else
						(Y + i * M)[m]=0;
				}
				norm = std::sqrt(norm);
				// std::cout<<norm<<'\n';
				for (uint64_t m = 0; m < M; m++)
				{
					(Y + i * M)[m] /= norm;
				}
				// Compute derivative
				diff_2_5points_allmesh((Y + i * M), M, (Yxx + i * M), h);

				for (uint64_t m = 0; m < M; m++)
				{
					rho_new[m]+=2 * (2 * l + 1) * 
							// Psi=R_nl*Y_lm
							// R_nl=u_nl/r   u is the Y here
							Y[i * M + m] * Y[i * M + m] / ((a + m * h) * (a + m * h));
				}
				itot++;
			}
		}

		// Compute rho
		for (uint64_t m = 0; m < M; m++)
		{
			rho[m] = (1 - alpha) * rho[m] + alpha * rho_new[m];
		}
		
		// Compute the two differents values of energy to check the convergence
		Eeigen = DFT_Eeigen(rho, m_explode_min, a, h, Nlevels, ls, Es, 
			e_c, de_c, e_x, de_x, U_r);
		Efunc = DFT_Efunc(rho, Y, Yxx, m_explode_min, a, h, Nlevels, ls, rhoB, Rc, 
			e_c, e_x, U_r, V_ext);
		
		std::cout << "\tE_eigen: "<< Eeigen << " " <<"E_func: "<< Efunc << std::endl;
		
		m_explode_min_before=m_explode_min;
		
		step++;
	}

	{
	std::ofstream file{"data/numerov_sc.dat"};
	file << "r rho V_l V_c V_x U_R V_Ext\n";
	for (uint64_t m = 0; m < M; m++)
	{
		file << a + m * h;
		// for (int i = 0; i <= 3; i++)
		// {
		// 	rho[i*M+m] = (i==0?0:rho[(i-1)*M+m])+ Y[i*M+m]*Y[i*M+m];
		// 	file << ' ' << Y[i*M+m] << ' '<< rho[i*M+m]; // degeneracy
		// }
		file<<' '<<rho[m]<<' '<<V_0[m] << ' '<< V_c[m] << ' ' << V_x[m] << ' ' << U_R[m] << ' ' << V_Ext[m];
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