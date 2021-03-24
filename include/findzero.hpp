#ifndef __FINDZERO_HPP__
#define __FINDZERO_HPP__

#include <cstdint>
#include <cmath>
#include <functional>

double findzero_secants_xeps(
	std::function<double(double)> F,
     double x0,double x1,const double epsilon,
	const double xmin,const double xmax);


#endif // __FINDZERO_HPP__