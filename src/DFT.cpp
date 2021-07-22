#include <DFT.hpp>

double DFT_Eeigen(
	const double rho[], const double rhor[], const double rhor2[], const uint64_t M, const double a, const double h,
	const uint8_t Nlevels, const int ls[], const double es[],
	std::function<double(double)> e_c, std::function<double(double)> e_x,
	std::function<double(double)> V_xc,
	std::function<double(double, double, double, size_t, const double[], const double[])> U)
{
	double E = 0;

	for (uint8_t i = 0; i < Nlevels; i++)
	{
		E += 2 * (2 * ls[i] + 1) * es[i];
	}

	std::function<double(size_t)> integrand = [&](size_t m) -> double
	{
		return (a + m * h) * (a + m * h) * rho[m] * (+e_x(rho[m]) + e_c(rho[m]) - V_xc(rho[m]) - 0.5 * U(a + m * h, a, h, M, rhor, rhor2));
	};
	E += 4 * M_PI * integrator_simpson_cubic(integrand, M, h);

	return E;
}


double DFT_Efunc(
	const double rho[], const double rhor[], const double rhor2[], const double y[], const double yxx[], const uint64_t M,
	const double a, const double h, const uint8_t Nlevels, const int ls[],
	const double rhoB, const double Rc,
	std::function<double(double)> e_c, std::function<double(double)> e_x,
	std::function<double(double, double, double, size_t, const double[], const double[])> U,
	std::function<double(const double, const double, const double)> V_ext)
{
	double E = 0;
	
	for (uint8_t i = 0; i < Nlevels; i++)
	{
		std::function<double(size_t)> integrand = [&](size_t m) -> double
		{
			return 0.5 * (-y[i * M + m] * yxx[i * M + m] + ls[i] * (ls[i] + 1) / ((a + h * m) * (a + h * m)) * y[i * M + m] * y[i * M + m]);
		};
		E += 2 * (2 * ls[i] + 1) * integrator_simpson_cubic(integrand, M, h);
	}


	std::function<double(size_t)> integrand = [&](size_t m) -> double
	{
		return (a + m * h) * (a + m * h) * rho[m] * (
			+ V_ext((a + m * h), rhoB, Rc) 
			+ 0.5 * U(a + m * h, a, h, M, rhor, rhor2) 
			+ e_x(rho[m]) + e_c(rho[m]));
	};
	E += integrator_simpson_cubic(integrand, M, h);
	return 4 * M_PI * E;
}