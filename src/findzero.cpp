#include "findzero.hpp"

#include <iostream>

double findzero_secants_xeps(
	std::function<double(double)> F, double x0,double x1,const double xeps,
	const double xmin,const double xmax){

	double fx0;
	double fx1=F(x0);
	double corr=xeps+1; 

	while(x1!=x0 && std::abs(corr) > xeps){
		fx0=fx1;
		fx1=F(x1);
		corr = -(x1 - x0) / (fx1 - fx0) * fx1;

		x0=x1; 
		x1=x1+corr; 
		x1=(x1>xmax)?xmax:((x1<xmin)?xmin:x1);
	}
	return x1;
}

double findzero_secants_yeps(
	std::function<double(double)> F, double x0, double x1, const double yeps,
	const double xmin, const double xmax)
{

	double fx1 = F(x0);
	double fx0=fx1+yeps+1;
	double corr;

	while (x1 != x0 && std::abs(fx1) > yeps)
	{
		fx0 = fx1;
		fx1 = F(x1);
		corr = -(x1 - x0) / (fx1 - fx0) * fx1;

		x0 = x1;
		x1 = x1 + corr;
		x1 = (x1 > xmax) ? xmax : ((x1 < xmin) ? xmin : x1);
		
	}
	return x1;
}