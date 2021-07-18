#pragma once

#include "vec3d.hpp"
#include "functions.hpp"
#include "differentiate.hpp"
#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <functional>

/* Init the cubic lattice positions with:
 *  L box width
 *  n cells along the side
 *  q particles for every cell
 *
 * Return:
 *  0 on succesful init
 *
 * NOTE:
 *  q determines the configuration used
 */
#define INIT_LATTICE_N_TOO_SMALL -3
#define INIT_LATTICE_Q_NOT_SUPPORTED -4
int init_lattice(Vec3D pos[],const size_t N, double L = 1,
				 uint32_t n = 1, uint32_t q = 1);

// Periodic conditions for box of length L
void apply_periodic_bounds(Vec3D pos[], const size_t N, double L);

// Compute the alias of other for interaction on me
// in the condition of periodic boundaries of lenght L
double compute_alias(const Vec3D &me, Vec3D &other, const double L);

/*
Does the interaction between the particles in my system
*/
template <typename IntCall,
		  typename IntSetup>
int do_interactions(
	const Vec3D pos[], size_t N,
	IntCall *interaction,
	IntSetup *start_interaction = nullptr,
	const double L = 1.)
{
	// TODO: maybe check for return codes
	double d;
	Vec3D alias;

	for (size_t i{0}; i < N; i++)
	{
		// House keeping at the begin
		if (start_interaction)
			(*start_interaction)(i);

		// Interaction on i
		// caused by all other particles
		for (size_t j{0}; j < N; j++)
		{
			// Skip the self case
			if (i != j)
			{
				// Copy the particle j
				alias = pos[j];

				// Check if there is interaction
				// and also compute the correct alias
				if ((d = compute_alias(pos[i], alias, L)) > 0)
				{
					(*interaction)(i, j, alias, d);
				}
				// Else: no interaction, un-perfect packing of spheres
			}
		}
	}
	return 0;
}

#if 0

// Give nomal distributed components to vector vec
void init_distribute_maxwell_boltzmann(Vec3D &vec, const double std = 1);
void init_distribute_maxwell_boltzmann(
	std::vector<Vec3D> &vel, const double std = 1);
// Velocity verlet
// namespace VV
// {
void new_position(
	std::vector<double> &pos1,
	const std::vector<double> &pos0,
	const std::vector<double> &vel0,
	const std::vector<double> &acc0,
	double dt, double L);
double new_velocity(
	std::vector<double> &vel1,
	const std::vector<double> &vel0,
	const std::vector<double> &acc0,
	const std::vector<double> &acc1,
	double dt);

// Compute the new accelerations and return the potential energy
double new_acceleration(
	// Where to store the accelerations
	std::vector<Vec3D> &acc,
	// Positions of the particles
	const std::vector<Vec3D> &pos,
	// Potential of interaction
	const std::function<double(double)> &Vr,
	// Mass of single particle (all equal)
	const double m,
	// Box length
	const double L,
	// Optional derivative of the potential
	const std::function<double(double)> &dVr_arg = nullptr);
// } // namespace VV


// Useful instantiation
extern const std::vector<Vec3D> aliaser;

#define AUTCO_EARLY_DATA 1
#define AUTCO_LATE_DATA 2
// AutoCorr manages the calculation of autocorrelations
class AutoCorr
{
public:
	AutoCorr() = delete;
	AutoCorr(const AutoCorr &) = delete;

	// Instantiate a new auto-correlation:
	// c_istart: starting index for correlation, after initial fluctuation
	// c_istop: total length of the simulation
	// c_length: correlation length
	// Nq: number of particles
	AutoCorr(
		uint32_t c_istart,
		uint32_t c_istop,
		uint32_t c_length,
		uint32_t Nq);

	// Single addition of data to the correlation
	// q: vector of new points
	// i: "time" index of the simulation
	int add_data(
		const std::vector<Vec3D> &q,
		uint32_t i);

	// Overloading for final output of the data
	double operator[](size_t i) const;

	uint32_t get_means_number() const;

private:
	// Starting index for correlation, after initial fluctuation
	const uint32_t c_istart;
	// Total length of the simulation
	const uint32_t c_istop;
	// Correlation length
	const uint32_t c_length;
	// Number of particles
	const uint32_t Nq;
	// Final number of means that will be done on the correlation
	// taking in account the
	const uint32_t c_means_number;

	// Actual correlation vector (the output)
	std::vector<double> c;
	// Storage for all the initial conditions step by step
	std::vector<Vec3D> c_0;
	// Quick save of also the RMS value fo the initial conditions
	std::vector<double> c_02m;
};

/*
Compute the autocorrelation function
x: series to computhe the autocorrelation of
corr: resulting autocorrelation function,
	its length is taken as the correlation length
*/
template <typename T>
void autocorrelation(std::vector<T> &corr,
					 const std::vector<T> &x)
{
	// Mean of x
	T m{mean(x)};
	// deviation from mean, numerator, denominator
	T xim, n, d;

	for (size_t t = 0; t < corr.size(); t++)
	{
		// Autocorrelation at offset t

		n = 0; // Numerator
		d = 0; // Denominator

		// Sum all the contributions
		for (size_t i = 0; i < x.size() - t; i++)
		{
			xim = x[i] - m;
			n += xim * (x[i + t] - m);
			d += xim * xim;
		}

		// Save
		corr[t] = n / d;
	}
}

#endif
