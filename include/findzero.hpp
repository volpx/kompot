#pragma once

#include <cstdint>
#include <cmath>
#include <functional>

double findzero_secants_xeps(
	std::function<double(double)> F,
     double x0,double x1,const double xeps,
	const double xmin,const double xmax);

double findzero_secants_yeps(
	std::function<double(double)> F,
	 double x0, double x1, const double yeps,
	const double xmin, const double xmax);
