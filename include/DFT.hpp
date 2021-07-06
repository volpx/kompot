#pragma once

#include <cstdint>
#include <functional>


double DFT_Eeigen(
	const double rho[],const uint64_t M, const double a, const double h,
	const uint8_t Nlevels,const int ls[],const double es[],
	std::function<double(double)> e_c,std::function<double(double)> de_c,
	std::function<double(double)> e_x,std::function<double(double)> de_x,
	std::function<double(double, double, uint64_t,const double[])> U);

double DFT_EMF(
	const double rho[],const double y[], const double yxx[], const uint64_t M, 
	const double a, const double h, const uint8_t Nlevels,const int ls[], 
	const double rhoB, const double Rc,
	std::function<double(double)> e_c,	std::function<double(double)> e_x,
	std::function<double(double, double, uint64_t,const double[])> U,
	std::function<double(const double,const  double,const  double)> V_ext);