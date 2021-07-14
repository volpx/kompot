#include "integrate.hpp"

double integrator_simpson_cubic(
	const double y[],const uint64_t M,const double h)
{
	double I=0;

	// for (uint64_t i{0}; i<M; i++)
	// {
	// 	I+=h*y[i];
	// }
	// return I;

	for (uint32_t i{1}; i<M-1; i+=2){
		I+=(y[i-1]+4*y[i]+y[i+1]);
	}
	if (M%2==1){
		I+=y[M-1]*(5./4);
		I+=y[M-2]*(5./2);
		I-=y[M-3]/4;
	}

	return I*h/3;
}