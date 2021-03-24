#include "findzero.hpp"

double findzero_secants_xeps(
	std::function<double(double)> F, double x0,double x1,const double epsilon,
	const double xmin,const double xmax){

	double corr{-F(x1)/( (F(x1)-F(x0))/(x1-x0) )}; 

	while(x1!=x0 && std::abs(corr) > epsilon){
		corr=-F(x1)/( (F(x1)-F(x0))/(x1-x0) );
		x0=x1; 
		x1=x1+corr; 
		x1=x1>xmax?xmax:x1<xmin?xmin:x1;
	}
	return x1;
}