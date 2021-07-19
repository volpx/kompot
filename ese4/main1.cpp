#include "uniconst.hpp"
#include "vec3d.hpp"
#include "simulation.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <static_math/cmath.h>

// LJ potential
double V_LJ(const double r)
{
	return 4 * (std::pow(r, -12.0) - std::pow(r, -6.0));
}

// WF
double WF(
	const Vec3D pos[], const size_t N,
	const double b,const double L)
{
	double sum = 0;
	Vec3D other{0,0,0};
	double r_ij;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			other = pos[j];

			if ((r_ij = compute_alias(pos[i], other, L)) > 0)
			{
				sum += -0.5 * std::pow(b / r_ij, 5.0);
			}
		}
	}
	return std::exp(sum);
}

double Density(
	const Vec3D pos[], const size_t N,
	const double b, const double L)
{
	// std::cout<<WF(pos, N, b)<<std::endl;
	double wf = WF(pos, N, b, L);
	return wf * wf;
}

// Local kinetic energy
double get_local_T(
	const Vec3D pos[], const uint32_t N,
	const double b, const double L)
{
	double T = 0;
	double r_li_s = 0;
	const double c = 30 * std::pow(b, 5.0);
	const double c1 = -c / 6; // -5b^5
	Vec3D other{0, 0, 0};
	Vec3D sum{0, 0, 0};
	double tmp;

	for (uint32_t l = 0; l < N; l++)
	{
		sum.clear();
		for (uint32_t i = 0; i < N; i++)
		{
			if (l != i)
			{
				other = pos[i];
				// r_li_s = (pos[l]-pos[i]).norm();
				r_li_s = compute_alias(pos[l], other, L);
				if (r_li_s > 0)
				{
					tmp = std::pow(r_li_s, -7.0);
					// Laplacian term
					T += -0.5 * (c + 2 * c1) * tmp;
					// Gradient term
					sum += tmp * (pos[l] - other);
				}
			}
		}
		T += 0.25 * c1 * c1 * sum.norm2();
	}

	return -T;
}

//Jackson-Feenberg kinetic energy
double get_local_TJF(
	const Vec3D pos[], const uint32_t N,
	const double b, const double L)
{
	double T = 0;
	double r_li_s = 0;
	double c = 30 * std::pow(b, 5.0);
	double c1 = -c / 6; // 5b^5
	Vec3D other;

	for (uint32_t l = 0; l < N; l++)
	{
		for (uint32_t i = 0; i < N; i++)
		{
			if (l != i)
			{
				other = pos[i];

				r_li_s = compute_alias(pos[l], other, L);
				if (r_li_s > 0)
				{
					T += (c + 2 * c1) * std::pow(r_li_s, -7.0);
				}
			}
		}
	}

	return 0.25 * T;
}

double get_potential(
	const Vec3D pos[], const uint32_t N,
	const double L)
{
	double d;
	Vec3D alias;
	double V = 0;

	// Sum on the pairs
	for (size_t i{0}; i < N; i++)
	{
		for (size_t j{0}; j < i; j++)
		{
			// Copy the particle j
			alias = pos[j];

			// Check if there is interaction
			// and also compute the correct alias
			if ((d = compute_alias(pos[i], alias, L)) > 0)
			{
				V += V_LJ(d);
			}
			// Else: no interaction, un-perfect packing of spheres
		}
	}
	return V;
}

int main()
{
	std::cout << "VMC: 4He superfluideseese" << std::endl;

	// Problem parameters
	constexpr double epsilon = uni::kB_r * 10.22;  // J
	constexpr double sigma = 0.2556e-9;			   // m
	constexpr double rho_experimental = 21.86e27;  // #/m^3
	constexpr double m = 4.002602 * uni::Da_to_kg; // Kg

	// Number of VMC steps
	constexpr uint32_t S = 1e6;
	// Thermalization steps
	constexpr uint32_t Sth = S / 10;
	constexpr uint32_t Sdata = S - Sth;
	// Number of particles
	constexpr uint32_t N = 32;
	// Averaging length
	constexpr uint32_t K = 1e3;

	// Derived parameters
	constexpr double hbar2_2m = uni::hbar_r * uni::hbar_r /
								(2 * m * epsilon * sigma * sigma);
	const double L = std::cbrt(N / rho_experimental) / sigma;
	const double b0 = std::pow(16. / (25 * hbar2_2m), 0.1);

	// Number of variational parameters
	constexpr uint32_t M = 20;
	// Parameter mesh
	const double ba = b0;
	const double bb = b0*1.1;
	const double bh = (bb-ba)/M;
	double bs[M];
	for (uint32_t i = 0; i < M; i++)
	{
		bs[i] = ba + i * bh;
	}
	
	const double deltamax=L/2;
	const double deltamin=1e-4;
	double delta = 0.4;
	

	std::cout
		<< "Problem constants:"
		<< "\nS: " << S
		<< "\nN: " << N
		<< "\nepsilon: " << epsilon
		<< "\nsigma: " << sigma
		<< "\nm: " << m
		<< "\nrho_exp: " << rho_experimental
		<< "\nL: " << L
		<< "\nhbar2_2m: " << hbar2_2m
		<< "\nb0: " << b0
		<< "\ndelta0: " << delta
		<< "\n"
		<< std::endl;

	// Final quantities (sums)
	double T = 0;
	double TJF = 0;
	double V = 0;
	WArray Ts(Sdata, Sth);
	WArray TJFs(Sdata, Sth);
	WArray Vs(Sdata, Sth);

	// Variables
	Vec3D *posa = new Vec3D[N];
	Vec3D *posb = new Vec3D[N];
	Vec3D *pos1 = posa;
	Vec3D *pos2 = posb;
	Vec3D *tmp = posa;
	double new_prob, b, P1, P2;
	CArray probs(K);

	for (uint32_t m = 0; m < M; m+=M)
	{
		b = bs[m];
		std::cout << "b: " << b << std::endl;
		// Particles initialization
		init_lattice(pos1, N, L, 2, 4);
		apply_periodic_bounds(pos1, N, L);
		P1 = Density(pos1, N, b, L);

		Ts.fill(0);
		TJFs.fill(0);
		Vs.fill(0);

		// Compute new quantities
		T = hbar2_2m * get_local_T(pos1, N, b, L);
		TJF = hbar2_2m * get_local_TJF(pos1, N, b, L);
		V = get_potential(pos1, N, L);
		std::cout 
			<< "\nInitial K_L: " << T
			<< "\nInitial K_JF: " << TJF
			<< "\nInitial potential: " << V
			<< std::endl;

		// Parameter m
		for (uint32_t s = 0; s < S; s++)
		{
			// Step s
			// Print progress
			if(s % (S/100)==0)
			{
				std::cout << s << ' ' << average(probs.data, K) << ' ' << delta << std::endl;
			}

			// Previus positions in pos1
			// New positions in pos2
			for (uint32_t i = 0; i < N; i++)
			{
				pos2[i].x = pos1[i].x + delta * (randu() - 0.5);
				pos2[i].y = pos1[i].y + delta * (randu() - 0.5);
				pos2[i].z = pos1[i].z + delta * (randu() - 0.5);
			}
			// Periodic boundary codition
			apply_periodic_bounds(pos2, N, L);

			// Acception-rejection method
			P2 = Density(pos2, N, b, L);
			new_prob = P2 / P1;
			probs[s] = new_prob;
			// Touchup the delta
			if (s < Sth && s % K == 0)
			{
				delta *= (1 + 1e-2 * (average(probs.data, K) / 0.5 - 1));
				delta = (delta>deltamax)? deltamax:(delta<deltamin)?deltamin:delta;
				// std::cout << s<< ' '<< average(probs.data, K) << ' ' << delta << std::endl;
			}

			if (randu() < new_prob)
			{
				// Accept the new state

				// Compute new quantities
				T = hbar2_2m * get_local_T(pos2, N, b, L);
				TJF = hbar2_2m * get_local_TJF(pos2, N, b, L);
				V = get_potential(pos2, N, L);

				tmp = pos2;
				pos2 = pos1;
				pos1 = tmp;
				P1 = P2;
			}

			// Sample the quantities (after "thermalization")
			if (Ts.inWindow(s))
			{
				Ts[s] = T;
				TJFs[s] = TJF;
				Vs[s] = V;
			}
		}

		// TODO: compute variance knowing that points are not perfectly independent
		std::cout
			<< "Delta: " << delta << " Probs:" << average(probs.data,K)
			<< "\nLocal kinetic energy: " << average(Ts.data, Sdata) << " +/- " << stddev(Ts.data, Sdata,1)
			<< "\nJF    kinetic energy: " << average(TJFs.data, Sdata) << " +/- " << stddev(TJFs.data, Sdata,1)
			<< "\nPotential     energy: " << average(Vs.data, Sdata) << " +/- " << stddev(Vs.data, Sdata,1)
			<< "\n"
			<< std::endl;
	}

	{
		std::ofstream file{"data/T.dat"};
		file << "step T TJf prob" << '\n';
		for (uint32_t s = 0; s < S; s++)
		{
			file << s << ' ' << Ts[s] << ' ' << TJFs[s] << ' ' << probs.get(s) << '\n';
		}
		file.close();
	}

	// interaction cutoff
	// autocorrelation

	return 0;
}