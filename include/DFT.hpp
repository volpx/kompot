#pragma once

#include "integrate.hpp"

#include <cstdint>
#include <functional>
#include <cmath>


double DFT_Eeigen(
	const double rho[],const double rhor[],const double rhor2[], const uint64_t M, const double a, const double h,
	const uint8_t Nlevels, const int ls[], const double es[],
	std::function<double(double)> e_c, std::function<double(double)> e_x,
	std::function<double(double)> V_xc,
	std::function<double(double, double,double,size_t, const double[], const double[])> U);

double DFT_Efunc(
	const double rho[],const double rhor[],const double rhor2[], const double y[], const double yxx[], const uint64_t M, 
	const double a, const double h, const uint8_t Nlevels,const int ls[], 
	const double rhoB, const double Rc,
	std::function<double(double)> e_c,	std::function<double(double)> e_x,
	std::function<double(double, double,double,size_t, const double[], const double[])> U,
	std::function<double(const double,const  double,const  double)> V_ext);
	