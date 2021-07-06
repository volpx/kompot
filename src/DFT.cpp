#include <DFT.hpp>

double DFT_Eeigen(
	const double rho[],const uint64_t M, const double a, const double h,
	const uint8_t Nlevels,const int ls[],const double es[],
	std::function<double(double)> e_c,std::function<double(double)> de_c,
	std::function<double(double)> e_x,std::function<double(double)> de_x,
	std::function<double(double, double, uint64_t,const double[])> U)
{
	double E = 0;

	for (uint8_t i = 0; i < Nlevels; i++)
	{
		E += 2 * (2 * ls[i] + 1) * es[i];
	}
	
	for(uint64_t m=0; m<M; m++)
	{
		// E_ex
		E += (e_c(rho[m]) + e_x(rho[m])) * rho[m] * h;
		
		// Direct term
		E += -(de_c(rho[m])+de_x(rho[m]) + 0.5 * U(a + h *m,h,M,rho))*rho[m] * h;
	}

	// E_x E_c

	return E;
}

double DFT_EMF(
	const double rho[],const double y[], const double yxx[], const uint64_t M, 
	const double a, const double h, const uint8_t Nlevels,const int ls[], 
	const double rhoB, const double Rc,
	std::function<double(double)> e_c,	std::function<double(double)> e_x,
	std::function<double(double, double, uint64_t,const double[])> U,
	std::function<double(const double,const  double,const  double)> V_ext)
{
	double E = 0;
	double tmp = 0;
	
	for (uint8_t i = 0; i < Nlevels; i++)
	{
		// integral
		tmp=0;
		for (uint64_t m = 0; m < M; m++)
		{
			tmp += 0.5 * (-y[i * M + m] * yxx[i * M + m] + ls[i] * (ls[i] + 1) / ((a + h * m) * (a + h * m)) * y[i * M + m] * y[i * M + m]) * h;
		}
		E += 2 * (2 * ls[i] + 1) * tmp;
	}
	
	for(uint64_t m=0; m<M; m++)
	{
		// E_ex
		E += (e_c(rho[m]) + e_x(rho[m])) * rho[m] * h;

		// Direct term
		E += (V_ext(a + h * m, rhoB, Rc) + 0.5 * U(a + h * m, h, M, rho)) * rho[m] * h;
	}

	return E;
}