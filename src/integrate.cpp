#include "integrate.hpp"

double integrator_simpson_cubic(
	const double y[],const size_t M,const double h)
{
	double res = 0;
	switch (M){
	case 0:
		break;
	case 1:
		res=h*y[0];
		break;
	default:
		res = y[0] + y[M - 1];
		for (size_t i = 1; i <= M / 2 - 1; i++)
		{
			res += 4 * y[2 * i - 1] + 2 * y[2 * i];
		}
		res += 4 * y[M - 2];
		res *= h / 3;
		break;
	}
	return res;
}

double integrator_simpson_cubic(
	std::function<double(size_t)> y, const size_t M, const double h)
{
	double res = 0;
	switch (M){
	case 0:
		break;
	case 1:
		res=h*y(0);
		break;
	default:
		res = y(0) + y(M - 1);
		for (size_t i = 1; i <= M / 2 - 1; i++)
		{
			res += 4 * y(2 * i - 1) + 2 * y(2 * i);
		}
		res += 4 * y(M - 2);
		res *= h / 3;
		break;
	}
	return res;
}
