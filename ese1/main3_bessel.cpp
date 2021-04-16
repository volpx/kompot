#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <fstream>

double j_m1(double x)
{
	return std::cos(x) / x;
}

double j_0(double x)
{
	return std::sin(x) / x;
}

double n_m1(double x)
{
	return std::sin(x) / x;
}

double n_0(double x)
{
	return -std::cos(x) / x;
}


double jl_true(int l,double x){
	double res=0;
	switch (l)
	{
	case -1:
		res = std::cos(x) / x;
		break;
	case 0:
		res = std::sin(x) / x;
		break;
	case 1:
		res = std::sin(x) / (x*x) - std::cos(x)/x;
		break;
	case 2:
		res = (3/(x*x)-1)*std::sin(x)/x-3*std::cos(x) / (x*x);
		break;
	case 3:
		res = (15/(x*x*x)-6/x)*std::sin(x) / x-(15/(x*x)-1)*std::cos(x)/x;
		break;
	default:
		break;
	}
	return res;
}
double yl_true(int l,double x){
	double res=0;
	switch (l)
	{
	case -1:
		res = std::sin(x) / x;
		break;
	case 0:
		res = -std::cos(x) / x;
		break;
	case 1:
		res = -std::cos(x) / (x*x) - std::sin(x)/x;
		break;
	case 2:
		res = -(3/(x*x)-1)*std::cos(x)/x-3*std::sin(x) / (x*x);
		break;
	case 3:
		res = -(15/(x*x*x)-6/x)*std::cos(x) / x-(15/(x*x)-1)*std::sin(x)/x;
		break;
	default:
		break;
	}
	return res;
}


int main()
{
	constexpr uint64_t N = 1e2;
	constexpr double a = 0. + 1. / N;
	constexpr double b = 10.;
	constexpr double h = (b - a) / N;
	constexpr uint64_t lmax = 3;

	std::vector<double> x(N);
	map(x, [a, h](uint64_t i) -> double { return a + i * h; });

	std::vector<double> j((lmax + 2) * N);
	std::vector<double> n((lmax + 2) * N);
	std::vector<double> jerr((lmax + 2) * N);
	std::vector<double> nerr((lmax + 2) * N);
	// Initialize l=-1
	for (uint64_t i = 0; i < N; i++)
	{
		j[0 * N + i] = j_m1(x[i]);
		// jerr[0 * N + i] = j[0 * N + i] - std::cos(x) / x;
		n[0 * N + i] = n_m1(x[i]);
		// nerr[0 * N + i] = n[0 * N + i] - std::sin(x) / x;
	}
	// Initialize l=0
	for (uint64_t i = 0; i < N; i++)
	{
		j[1 * N + i] = j_0(x[i]);
		
		n[1 * N + i] = n_0(x[i]);
	}

	for (uint64_t k = 2; k <= lmax + 1; k++)
	{
		// Computing this l
		uint64_t l = k - 1;
		for (uint64_t i = 0; i < N; i++)
		{
			// Compute j-functions
			j[k * N + i] = (2 * (l - 1) + 1) / x[i] * j[(k - 1) * N + i] - j[(k - 2) * N + i];
			// j_gsl[k * N + i] = gsl_sf_bessel_jl(l, x[i]);
			jerr[k * N + i] = j[k * N + i] - jl_true(l, x[i]);

			// Compute n-functions
			n[k * N + i] = (2 * (l - 1) + 1) / x[i] * n[(k - 1) * N + i] - n[(k - 2) * N + i];
			// n_gsl[k * N + i] = gsl_sf_bessel_yl(l, x[i]);
			nerr[k * N + i] = n[k * N + i] - yl_true(l, x[i]);
		}
	}

	std::ofstream file{"data/bessel_err.dat"};
	file << "#x l=1 jerr nerr l=2 jerr nerr\n";
	// Write to file
	for (uint64_t i = 0; i < N; i++)
	{
		file << x[i];
		for (uint64_t k = 2; k <= lmax + 1; k++)
		{
			file << ' ' << jerr[k * N + i]
				 << ' ' << nerr[k * N + i];
		}
		file << '\n';
	}
	file.close();

	return 0;
}