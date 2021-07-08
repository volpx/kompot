#include <numerov.hpp>
#include <iostream>

// Normal numerov integration using y[0],y[1] as initial conditions
void numerov_integrate(
	double y[],
	const double a, const double h,
	const uint64_t N, std::function<double(double)> V,
	const double E)
{
	const double h2 = h * h;
	double k2, k_i2, k_ii2;
	k_ii2 = 2. * (E - V(a));
	k_i2 = 2. * (E - V(a + h));

	for (uint64_t i = 2; i < N; i++)
	{
		// New k
		k2 = 2. * (E - V(a + i * h));

		y[i] = (2. * y[i - 1] * (1. - 5. / 12 * h2 * k_i2) -
				y[i - 2] * (1. + 1. / 12 * h2 * k_ii2)) /
			   (1. + h2 / 12.0 * k2);

		k_ii2 = k_i2;
		k_i2 = k2;
	}
}
void numerov_integrate(
	double y[], const double h,
	const uint64_t N, const double V[],
	const double E)
{
	const double h2 = h * h;
	double k2, k_i2, k_ii2;
	k_ii2 = 2. * (E - V[0]);
	k_i2 = 2. * (E - V[1]);

	for (uint64_t i = 2; i < N; i++)
	{
		// New k
		k2 = 2. * (E - V[i]);

		y[i] = (2. * y[i - 1] * (1. - 5. / 12 * h2 * k_i2) -
				y[i - 2] * (1. + 1. / 12 * h2 * k_ii2)) /
			   (1. + h2 / 12.0 * k2);

		k_ii2 = k_i2;
		k_i2 = k2;
	}
}

// Function that returns the value of y(xmax) for given E
double numerov_integrate_yxmax(
	double y_ii, double y_i,
	const double a, const double h,
	const uint64_t N, const std::function<double(double)> V,
	const double E)
{
	const double h2 = h * h;

	double y = 0;
	double k2, k_i2, k_ii2;

	k_ii2 = 2. * (E - V(a));
	k_i2 = 2. * (E - V(a + h));

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V(a + i * h));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);

		if (std::abs(y) > 1e100)
		{
			return y;
		}

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
	return y;
}

double numerov_integrate_yxmax(
	double y_ii, double y_i,
	const double h,
	const uint64_t N, const double V[],
	const double E)
{
	const double h2 = h * h;

	double y = 0;
	double k2, k_i2, k_ii2;

	k_ii2 = 2. * (E - V[0]);
	k_i2 = 2. * (E - V[1]);

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V[i]);

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);

		if (std::abs(y) > 1e100)
		{
			return y;
		}

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
	return y;
}

// Print best estimate of Energy eigenvalues using secants method
double numerov_find_energy(
	double y_ii, double y_i,
	const double a, const double h,
	const uint64_t N, std::function<double(double)> V,
	const double Ea, const double Eb, const double Eh)
{
	double E_after_secants = NAN;
	double y_xmax_ii = numerov_integrate_yxmax(y_ii, y_i, a, h, N, V, Ea);
	double y_xmax_i = 0;

	double E = Ea + Eh;
	while (E < Eb)
	{
		y_xmax_i = numerov_integrate_yxmax(y_ii, y_i, a, h, N, V, E);

		// Search for a change of sign
		if (y_xmax_ii * y_xmax_i < 0)
		{
			// Get a better value by secants method
			double E_after_secants = findzero_secants_xeps(
				[y_ii, y_i, a, h, N, V](double E) -> double { return numerov_integrate_yxmax(y_ii, y_i, a, h, N, V, E); },
				E, E - Eh, 1e-7, E - Eh, E);
			return E_after_secants;
		}

		y_xmax_ii = y_xmax_i;
		E += Eh;
	}

	return E_after_secants;
}

double numerov_find_energy(
	double y_ii, double y_i,
	const double h,
	const uint64_t N, const double V[],
	const double Ea, const double Eb, const double Eh)
{
	double E_after_secants = NAN;
	double y_xmax_ii = numerov_integrate_yxmax(y_ii, y_i, h, N, V, Ea);
	double y_xmax_i = 0;

	uint64_t i = 1;
	double E = Ea + i * Eh;
	while (E < Eb)
	{
		// if(!(i%100))
		// {
		// 	std::cout<<i<<"\tTrying E="<<E<<std::endl;
		// }
		y_xmax_i = numerov_integrate_yxmax(y_ii, y_i, h, N, V, E);

		// Search for a change of sign
		if (y_xmax_ii * y_xmax_i < 0)
		{
			// Get a better value by secants method
			double E_after_secants = findzero_secants_xeps(
				[y_ii, y_i, h, N, V](double E) -> double { return numerov_integrate_yxmax(y_ii, y_i, h, N, V, E); },
				E, E - Eh, 1e-7, E - Eh, E);
			return E_after_secants;
		}

		y_xmax_ii = y_xmax_i;
		i++;
		E = Ea + i * Eh;
	}

	return E_after_secants;
}