#include "simulation.hpp"

int init_lattice(Vec3D pos[], const size_t N,
				 const double L, const uint32_t n, const uint32_t q)
{
	if (q * n * n * n > N)
	{
		// too few particles spaces
		return INIT_LATTICE_N_TOO_SMALL;
	}
	if (!((q == 1) || (q == 2) || (q == 4)))
	{
		// Only q=1,2,4 is supported
		return INIT_LATTICE_Q_NOT_SUPPORTED;
	}
	// Length of a single cell
	const double dl{L / n};

	// Place the particles
	for (uint32_t i{0}; i < n; i++)
	{
		for (uint32_t j{0}; j < n; j++)
		{
			for (uint32_t k{0}; k < n; k++)
			{
				if (q == 1)
				{
					pos[q * (i * n * n + j * n + k)] =
						Vec3D{i * dl, j * dl, k * dl};
				}
				else if (q == 2)
				{
					pos[q * (i * n * n + j * n + k) + 0] =
						Vec3D{i * dl, j * dl, k * dl};
					pos[q * (i * n * n + j * n + k) + 1] =
						Vec3D{i * dl + dl / 2,
							  j * dl + dl / 2,
							  k * dl + dl / 2};
				}
				else if (q == 4)
				{
					pos[q * (i * n * n + j * n + k) + 0] =
						Vec3D{i * dl, j * dl, k * dl};
					pos[q * (i * n * n + j * n + k) + 1] =
						Vec3D{i * dl,
							  j * dl + dl / 2,
							  k * dl + dl / 2};
					pos[q * (i * n * n + j * n + k) + 2] =
						Vec3D{i * dl + dl / 2,
							  j * dl,
							  k * dl + dl / 2};
					pos[q * (i * n * n + j * n + k) + 3] =
						Vec3D{i * dl + dl / 2,
							  j * dl + dl / 2,
							  k * dl};
				}
			}
		}
	}
	return N - q * n * n * n;
}


inline int sign(double x)
{
	return x > 0 ? 1 : -1;
}
double compute_alias(const Vec3D &me, Vec3D &other, const double L)
{
	// I'd like to thank @zecco for this immensly valuable contribution
	if (std::fabs(me.x - other.x) > L / 2)
	{
		other.x += sign(me.x - other.x) * L;
	}
	if (std::fabs(me.y - other.y) > L / 2)
	{
		other.y += sign(me.y - other.y) * L;
	}
	if (std::fabs(me.z - other.z) > L / 2)
	{
		other.z += sign(me.z - other.z) * L;
	}

	double d{(me - other).norm()};
	// Return the distance
	// but a negative value if is >= L/2
	return d < L / 2 ? d : -d;
}

void apply_periodic_bounds(Vec3D pos[], const size_t N, double L)
{
	for (size_t i{0}; i < N; i++)
	{
		pos[i].x -= L * std::floor(pos[i].x / L);
		pos[i].y -= L * std::floor(pos[i].y / L);
		pos[i].z -= L * std::floor(pos[i].z / L);
	}
}

#if 0
void init_distribute_maxwell_boltzmann(Vec3D &vec, const double std)
{
	// Give a normal distribution to the coponents
	vec.x = std * randn();
	vec.y = std * randn();
	vec.z = std * randn();
}

// Distribute the velocities as MB then adjusts for CM motion and total std
void init_distribute_maxwell_boltzmann(
	std::vector<Vec3D> &vel, const double std)
{
	const size_t N{vel.size()};
	Vec3D cm{0, 0, 0};
	double exp_std{0};

	for (size_t i{0}; i < N; ++i)
	{
		// Give a normal distribution to the components
		vel[i].x = std * randn();
		vel[i].y = std * randn();
		vel[i].z = std * randn();

		// Sum to cm
		cm += vel[i];
	}

	// Corrections
	exp_std = 0;
	for (size_t i{0}; i < N; ++i)
	{
		vel[i] -= (1.0 / N) * cm;

		// overall energy
		exp_std += vel[i].norm2();
	}

	// Obtained std
	exp_std = std::sqrt(exp_std / (3 * N));

	// scale to match std
	for (size_t i{0}; i < N; ++i)
	{
		vel[i] /= (exp_std / std);
	}
}
double new_acceleration(
	std::vector<Vec3D> &acc,
	const std::vector<Vec3D> &pos,
	const std::function<double(double)> &Vr,
	const double m,
	const double L,
	const std::function<double(double)> &dVr_arg)
{
	// Needed for evenual 5 points differentiation
	constexpr double diff_eps{1e-10};

	// Number of particles
	size_t N{pos.size()};

	double d;
	Vec3D tmp_acc;
	Vec3D alias;
	std::function<double(double)> dVr;
	// Cumulative value to return
	double Vreturn{0};

	// Check if derivative has been supplied
	if (dVr_arg)
	{
		dVr = dVr_arg;
	}
	else
	{
		// Computational derivative
		dVr = [&](double d) -> double {
			return derive_5points(Vr, d, diff_eps);
		};
	}

	for (size_t i{0}; i < N; i++)
	{
		// Clear the acceleration before start adding contributions
		acc[i].clear();

		// Acceleration on i
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
					// Compute the derivative on function v around d
					// Compute the acc contribute for the pair ij
					// Add the contribute to the acceleration
					// using the difference versor
					acc[i] += -dVr(d) / m / d * (pos[i] - alias);
					// add potential on half of the cases
					if (i < j)
						Vreturn += Vr(d);
				}
				// Else: no interaction, un-perfect packing of spheres
			}
		}
	}
	return Vreturn;
}


const std::vector<Vec3D> aliaser{
	// in order of feasability
	Vec3D{0, 0, 0},

	Vec3D{0, 0, -1},
	Vec3D{0, 0, 1},
	Vec3D{0, -1, 0},
	Vec3D{0, 1, 0},
	Vec3D{-1, 0, 0},
	Vec3D{1, 0, 0},

	Vec3D{-1, -1, 0},
	Vec3D{1, 1, 0},
	Vec3D{-1, 0, -1},
	Vec3D{1, 0, 1},
	Vec3D{-1, 0, 1},
	Vec3D{-1, 1, 0},
	Vec3D{0, -1, -1},
	Vec3D{0, -1, 1},
	Vec3D{0, 1, -1},
	Vec3D{0, 1, 1},
	Vec3D{1, -1, 0},
	Vec3D{1, 0, -1},

	Vec3D{-1, -1, -1},
	Vec3D{-1, -1, 1},
	Vec3D{-1, 1, -1},
	Vec3D{-1, 1, 1},
	Vec3D{1, -1, -1},
	Vec3D{1, -1, 1},
	Vec3D{1, 1, -1},
	Vec3D{1, 1, 1},
};

// AutoCorr constructor
AutoCorr::AutoCorr(
	uint32_t c_istart,
	uint32_t c_istop,
	uint32_t c_length,
	uint32_t Nq)
	: c_istart{c_istart},
	  c_istop{c_istop},
	  c_length{c_length},
	  Nq{Nq},
	  c_means_number{c_istop - c_istart - c_length + 1},

	  c(c_length),
	  c_0(Nq * c_length),
	  c_02m(c_length)
{
	// Clear all the initial conditions
	for (uint32_t i{0}; i < Nq * c_length; i++)
	{
		c_0[i].clear();
	}
	// And also the results vector
	fill(c, 0);
}

// Handles the single addition of new data
int AutoCorr::add_data(
	const std::vector<Vec3D> &q,
	uint32_t i)
{
	// Check temporal position in simulation
	if (i >= c_istart && i < c_istop)
	{
		// wrapping index needed for cvv_v0 and cvv_v02m
		// it shouldn't really be necessary to subtract c_istart
		uint32_t c_i = (i - c_istart) % c_length;

		// Setup initial conditions of correlations
		// until the incomplete sequences at the end
		if (i < c_istop - c_length + 1)
		{
			// Normal operation
			// Calculate already the RMS
			c_02m[c_i] = 0;
			for (uint32_t j = 0; j < Nq; j++)
			{
				// Save the new c_0 offsetting location
				// using c_i and calculate the norm2
				c_0[c_i * Nq + j] = q[j];
				c_02m[c_i] += q[j].norm2();
			}

			// Mean
			// Can be avoided because cancels out
			// with the other mean at the numerator
			// c_02m[i] /= Nq;
		}
		else
		{
			// Last cvv_length steps
			// setting to zero will prevent
			// adding the corresponding incomplete contributes
			c_02m[c_i] = 0;
		}

		// Contribution on c calculations
		double c_tmp;
		uint32_t jjj;
		for (uint32_t jj = 0; jj < c_length; jj++)
		{
			// Contribution on c[jj]

			// The index of c_02m and c_0 corresponding to jj
			jjj = (c_i - jj + c_length) % c_length;

			// Skip the steps for which c_02m==0
			// they happens if I'm in the last steps for which
			// it has been hard set above
			// or can happen if the mean velocity of *ALL* the particles
			// is null either way it would skyrocket the correlation
			if (c_02m[jjj] != 0)
			{
				// Incremental mean
				c_tmp = 0;
				for (uint32_t j = 0; j < Nq; j++)
				{
					// Scalar product
					c_tmp += q[j] * c_0[jjj * Nq + j];
				}
				// Mean on particles
				// Avoided as stated before on computing c_02m,
				// decommenting *BOTH* lines would not change the result
				// c_tmp /= Nq;
				// Normalize on the <c0^2>
				c_tmp /= c_02m[jjj];
				// Add contribution
				c[jj] += c_tmp;
			}
		}
	}
	else if (i < c_istart)
	{
		return AUTCO_EARLY_DATA;
	}
	else if (i >= c_istop)
	{
		return AUTCO_LATE_DATA;
	}
	return 0;
}

// Getter for the data
double AutoCorr::operator[](size_t i) const
{
	// Means of the final correlation result
	return c[i] / c_means_number;
}

uint32_t AutoCorr::get_means_number() const
{
	return c_means_number;
}
#endif
