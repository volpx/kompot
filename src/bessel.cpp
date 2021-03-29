#include "bessel.hpp"

double sf_bessel_jl(int l, double x)
{
	if (l = -1)
	{
		return std::cos(x) / x;
	}
	else if (l = 0)
	{
		return std::sin(x) / x;
	}
	else if (l > 0)
	{
		return (2 * l + 1) / x * sf_bessel_jl(l - 1, x) - sf_bessel_jl(l - 2, x);
	}
	else
	{
		return 1.0 / 0;
	}
}

double sf_bessel_yl(int l, double x)
{
	if (l = -1)
	{
		return std::sin(x) / x;
	}
	else if (l = 0)
	{
		return -std::cos(x) / x;
	}
	else if (l > 0)
	{
		return (2 * l + 1) / x * sf_bessel_yl(l - 1, x) - sf_bessel_yl(l - 2, x);
	}
	else
	{
		return 1. / 0;
	}
}