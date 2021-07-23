#include "uniconst.hpp"
#include "vector_help.hpp"
#include "numerov.hpp"
#include "DFT.hpp"
#include "differentiate.hpp"
#include "integrate.hpp"
#include "functions.hpp"

#include <cstdint>
#include <cstddef>
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

// Rotational potential
double V_rot(
	const double r, const int l)
{
	return 0.5 * l * (l + 1) / (r * r);
}


// Direct term
double U_r(
	const double r, const double a, const double h,
	const size_t M, const double rhor[], const double rhor2[])
{
	double u = 0;
	const size_t m_split = static_cast<size_t>((r-a)/h);
	
	u += 1./r*integrator_simpson_cubic(rhor2,m_split+1,h);
	u += integrator_simpson_cubic(rhor+m_split,M-m_split,h);

	return 4 * M_PI * u;
}

// Exchange energy
// e_x(rho(r)) (scalar r)
double e_x(const double rho)
{
	return -3. / 4 * std::pow(3 * rho / M_PI, 1. / 3);
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
	// const double DEN = 2 * A * (beta[0] * std::pow(rs, 0.5) + beta[1] * rs + beta[2] * std::pow(rs, 3.0 / 2) + beta[3] * rs * rs);
	// return -2 * A * (1 + alpha_1 * rs) * std::log(1 + 1.0 / DEN);
	return -2.*A*(1.+alpha_1*rs) * 
		std::log(1.+1./(2. * A * rs * (
			beta[0]/std::sqrt(rs) + beta[1] +
			beta[2] * std::sqrt(rs) +
			beta[3] * std::pow(rs,1.)
		)));
}

// de_c/drho Derivative of correlation
double de_c(const double rho)
{

	// Correlation energy parameters
	constexpr double A = 0.031091;
	constexpr double alpha_1 = 0.21370;
	constexpr double beta[4] = {7.5957, 3.5876, 1.6382, 0.49294};

	// 
	const double rho13 = std::pow(rho,1./3);
	const double c1 = std::pow(3. / (4 * M_PI), 1. / 3);
	const double c2 = 2. * A * alpha_1 * c1;
	const double rs = c1 / rho13;
	// const double den = 2 * A * (beta[0] * std::sqrt(rs) + beta[1] * rs + beta[2] * std::pow(rs, 3. / 2) + beta[3] * rs * rs);
	const double den = 2. * A * rs *(
		beta[0]/std::sqrt(rs) + beta[1] +
		beta[2]* std::sqrt(rs) + beta[3] * std::pow(rs,1.)
	);
	const double c3 = 1. + 1./den;

	return c2/(3. * rho13*rho)*std::log(c3)
		- (2.*A+c2/(rho13))/(c3*den*den)*2.*A*
	(
		beta[0]*std::sqrt(c1) / (6.*std::pow(rho13,7./2))+
		beta[1]*c1/(3.*rho*rho13) +
		beta[2]*0.5*std::pow(c1/rho,3./2) +
		beta[3]*2.*c1*c1 / (3. *std::pow(rho13,5.))
	);
}

// Exchange correlation potential
double V_xc(const double rho)
{
	return 4./3 * e_x(rho) + e_c(rho) + de_c(rho) * rho;
}

// KS potential
double V_KS(const double r, const size_t i,
	const double rho[], const size_t M,
	const double rhor[], const double rhor2[],
	const double h, const double a,
	const double Rc,const double rhoB)
{
	return + U_r(r, a, h, M, rhor, rhor2) 
		   + (rho[i] == 0 ? 0 : V_xc(rho[i])) 
		   + V_ext(r, rhoB, Rc);
}

int main()
{
	std::cout << "DFT: Independent Electron Model" << std::endl;

	// Problem constants
	constexpr double r_s_Na = 3.93;
	constexpr double r_s_K = 4.86;
	constexpr double rs = r_s_Na;
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

	// Background properties
	constexpr double rhoB = 1. / (4. / 3 * M_PI * rs * rs * rs);
	const double Rc = std::pow(N, 1. / 3) * rs;

	// Spatial mesh
	constexpr double a = 1e-4;
	constexpr double b = 30;
	constexpr double h = 1e-3;
	constexpr size_t M = static_cast<size_t>((b-a)/h);
	size_t m_explode=M;

	// Potential
	double V[M];
	double V_0[M];
	double V_xcs[M];
	double U_rs[M];
	double V_exts[M];
	double rhor2[M];
	double rhor[M];

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
	constexpr double epsilon = 1e-4;

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


	// Variables
	int itot = 0;
	int i;
	double E;
	double Eeigen = 1.0;
	double Efunc = 0;
	double Ymin, norm;
	fill(rho, M, 0);
	uint64_t step = 0;

	// First potential
	for (uint64_t m = 0; m < M; m++)
	{
		V_0[m] = V_ext(a + m * h, rhoB, Rc);
	}

	// Start self consistent procedure
	while (std::abs(Eeigen - Efunc) > epsilon)
	{
		std::cout << "Step: " << step << std::endl;

		// Initialization
		itot = 0;
		fill(rho_new, M, 0);

		for (int l = 0; l <= lmax; l++)
		{
			// Solve for value l

			// Add rotational term
			for (uint64_t m = 0; m < M; m++)
			{
				V[m] = V_0[m] + V_rot(a + m * h, l);
			}

			E = Ea = min(V, M);
			Eh = (Eb - Ea) / ME;
			y0 = std::pow(a, l + 1);
			y1 = std::pow(a + h, l + 1);

			// Solve the differential eq with numerov for each state
			for (int g = 0; g < Enums[l]; g++)
			{
				// Actual occupation order
				i=is[itot];

				// Move up a bit from the previous value/minimum
				Ea=E+1e-10;
				std::cout 
					<< '\t' << l << g 
					<< " Ea: "<<Ea;
				E = Es[i] = numerov_find_energy(y0, y1, h, M, V, Ea, Eb, Eh);
				
				#if 0
				if(std::isnan(E))
				{
					std::cout<<"is nan\n";
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
				#endif

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

				// Compute derivative
				diff_2_5points_allmesh((Y + i * M), M, (Yxx + i * M), h);

				for (uint64_t m = 0; m < M; m++)
				{
					rho_new[m] += 2 * (2 * l + 1) * 
							// Psi=R_nl*Y_lm
							// R_nl=u_nl/r   u is the Y in this code
							Y[i * M + m] * Y[i * M + m] / ((a + m * h) * (a + m * h));
				}

				std::cout
					<< " E: " << Es[i] 
					<< " (" << Ymin << ")"
					<< std::endl;

				itot++;
			}
		}

		// Compute rho
		for (uint64_t m = 0; m < M; m++)
		{
			rho[m] = (1 - alpha) * rho[m] + alpha * rho_new[m];
			rhor[m] = rho[m]*(a+m*h);
			rhor2[m] = rhor[m]*(a+m*h);
		}

		// Compute new potential without the l dependancy
		for (uint64_t m = 0; m < M; m++)
		{

			V_0[m] = V_KS(a + m * h, m, rho, M, rhor, rhor2, h, a, Rc, rhoB);
		}
		
		// Compute the two differents values of energy to check the convergence
		Eeigen = DFT_Eeigen(rho,rhor,rhor2, M, a, h, Nlevels, ls, Es, e_c, e_x, 
			V_xc, U_r);
		Efunc = DFT_Efunc(rho,rhor,rhor2, Y, Yxx, M, a, h, Nlevels, ls, rhoB, Rc, 
			e_c, e_x, U_r, V_ext);
		
		std::cout << "\tE_eigen: "<< Eeigen << " " <<"E_func: "<< Efunc << std::endl;

		#if 0
		{
			std::ofstream file{"data/numerov_sc.dat"};
			file << "r rho V_l V_c V_x U_R V_Ext\n";
			for (uint64_t m = 0; m < M; m++)
			{
				V_xcs[m] = V_xc(rho[m]);
				U_rs[m] = U_r(a + h * m, a, h, M,rhor,rhor2);
				V_exts[m] = V_ext(a + m * h, rhoB, Rc);

				file
					<< a + m * h
					<< ' ' << rho[m] 
					<< ' ' << V_0[m] 
					<< ' ' << V_xcs[m] 
					<< ' ' << U_rs[m] 
					<< ' ' << V_exts[m]
					<< '\n';
			}
			file.close();
		}
		PressEnterToContinue();
		#endif

		step++;
	}

	{
		std::ofstream file{"data/numerov_sc.dat"};
		file << "r rho V_l V_c V_x U_R V_Ext\n";
		for (uint64_t m = 0; m < M; m++)
		{
			V_xcs[m] = V_xc(rho[m]);
			U_rs[m] = U_r(a + h * m, a, h, M,rhor,rhor2);
			V_exts[m] = V_ext(a + m * h, rhoB, Rc);

			file
				<< a + m * h
				<< ' ' << rho[m] 
				<< ' ' << V_0[m] 
				<< ' ' << V_xcs[m] 
				<< ' ' << U_rs[m] 
				<< ' ' << V_exts[m]
				<< '\n';
		}
		file.close();
	}

	size_t m_Rc = static_cast<size_t>((Rc - a) / h);
	for (size_t m = 0; m < M; m++)
	{
		rhor2[m]*=4*M_PI;
		if (m < m_Rc)
		{
			rhor2[m] = 0;
		}
	}

	double deltaN = integrator_simpson_cubic(rhor2, M, h);
	double alphaN = std::pow(Rc, 3.) * (1 + deltaN / N);
	std::cout << "DeltaN: " << deltaN << ' ' << "Alpha: " << alphaN << std::endl;

	return 0;
	
}