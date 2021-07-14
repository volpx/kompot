#include <DFT.hpp>

double DFT_Eeigen(
	const double rho[], const uint64_t M, const double a, const double h,
	const uint8_t Nlevels, const int ls[], const double es[],
	std::function<double(double)> e_c, std::function<double(double)> de_c,
	std::function<double(double)> e_x, std::function<double(double)> de_x,
	std::function<double(double, double, double, uint64_t, const double[])> U)
{
	double E = 0;

	for (uint8_t i = 0; i < Nlevels; i++)
	{
		E += 2 * (2 * ls[i] + 1) * es[i];
	}
	
	// E_x E_c
	// Integral
	double I = 0;
	for(uint64_t m=0; m<M; m++)
	{
		I += h * ((a + m * h) * (a + m * h) * rho[m] * (

				+ e_x(rho[m]) + e_c(rho[m]) 
				+ (
					+ e_x(rho[m]) + e_c(rho[m])
					+ rho[m] * ( de_c(rho[m]) + de_x(rho[m]) )
				  )
				- 0.5 * U(a + m * h, a, h, M, rho))
				);
	}
	E += 4 * M_PI * I;

	return E;
}

double DFT_Efunc(
	const double rho[],const double y[], const double yxx[], const uint64_t M, 
	const double a, const double h, const uint8_t Nlevels,const int ls[], 
	const double rhoB, const double Rc,
	std::function<double(double)> e_c,	std::function<double(double)> e_x,
	std::function<double(double, double,double, uint64_t,const double[])> U,
	std::function<double(const double,const  double,const  double)> V_ext)
{
	double E = 0;
	double I = 0;
	
	for (uint8_t i = 0; i < Nlevels; i++)
	{
		// integral
		I=0;
		for (uint64_t m = 0; m < M; m++)
		{
			I += h * 0.5 * (
				-y[i * M + m] * yxx[i * M + m]
				+ ls[i] * (ls[i] + 1) / ((a + h * m) * (a + h * m)) * y[i * M + m] * y[i * M + m]
				);
		}
		E += 2 * (2 * ls[i] + 1) * I;
	}
	
	I = 0;
	for(uint64_t m=0; m<M; m++)
	{
		I = h * ((a + m * h) * (a + m * h) * rho[m] * (
			+ V_ext((a + m * h), rhoB, Rc) 
			+ 0.5 * U(a + m * h, a, h, M, rho) 
			+ e_x(rho[m]) + e_c(rho[m]))
			);
	}
	E += I;
	return 4 * M_PI * E;
}