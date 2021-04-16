#include "vector_help.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_hermite.h>

double V(double x){
	return 0.5*x*x;
}

double y_true(double x, int n){
	return gsl_sf_hermite(n ,x)*std::exp(-0.5*x*x);
}

void save_function(const double E,const int even, const double a, const double b, const uint64_t N, std::ofstream &file)
{
	const double h = (b - a) / N;
	const double h2 = h * h;

	double y = 0, y_i, y_ii;
	double k2, k_i2, k_ii2;

	if (even==0){
		y_ii = 1 * std::exp(-0.5 * a *a);
		y_i = 1 * std::exp(-0.5 * (a + h) * (a + h));
	}
	else{
		y_ii = a;
		y_i = a+h;
	}
	k_ii2 = 2. * (E - V(a));
	k_i2 = 2. * (E - V(a+h));

	file << a << ' ' << y_ii << '\n';
	file << a + h << ' ' << y_i << '\n';

	for (uint64_t i = 2; i < N; i++)
	{

		// New k
		k2 = 2. * (E - V(a + i * h));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);

		file << a + i * h << ' ' << y << '\n';

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
}

double explode_at(const double E, const int even, const double a,const double b, const double h,
	const double ymax,const double xth)
{
	const double h2 = h * h;

	double y = 0, y_i, y_ii;
	double k2, k_i2, k_ii2;

	if (even == 0)
	{
		y_ii = 1 * std::exp(-0.5 * a * a);
		y_i = 1 * std::exp(-0.5 * (a + h) * (a + h));
	}
	else
	{
		y_ii = a;
		y_i = a + h;
	}
	k_ii2 = 2. * (E - V(a));
	k_i2 = 2. * (E - V(a + h));

	for (uint64_t i = 2; a+h*i<b; i++)
	{

		// New k
		k2 = 2. * (E - V(a + i * h));

		y = (2. * y_i * (1. - 5. / 12 * h2 * k_i2) -
			 y_ii * (1. + 1. / 12 * h2 * k_ii2)) /
			(1. + h2 / 12.0 * k2);
		
		if (a+i*h > xth && std::abs(y)>ymax){
			return a+i*h;
		}

		k_ii2 = k_i2;
		k_i2 = k2;
		y_ii = y_i;
		y_i = y;
	}
	return b;
}

int main1()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t M=1e3;      // Number of energies
	constexpr double Ea=0;
	constexpr double Eb=0.5+7*1+0.25;

	constexpr uint64_t N=1e8; // number of mesh points
	constexpr double a=0.+1./1e6;   // start included
	constexpr double b=10;   // stop excluded

	constexpr double h=(b-a)/N; // Mesh spacing
	constexpr double h2=h*h;

	std::cout << "Problem constants:"
			  << "\nM: " << M
			  << "\nEa: " << Ea
			  << "\nEb: " << Eb
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << std::endl;

	// Create an output of the wavefunction
	char filename[50];
	sprintf(filename, "data/eigen1d_mesh%ld.dat",N);
	std::ofstream file{filename};
	save_function(0.5,0, a, b, N, file);
	file.close();

	return 0;
}

int main()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t N = 1e6;			// number of mesh points
	constexpr double a = 0. + 1. / 1e6; // start included
	constexpr double b = 10;			// stop excluded

	constexpr double h = (b - a) / N; // Mesh spacing
	constexpr double h2 = h * h;

	std::cout << "Problem constants:"
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << std::endl;

	const double ha=1./1e7;
	const double hb=1./1e1;
	const uint64_t M=1000;
	std::vector<double> xexpl(M);

	std::ofstream file{"data/mesh_var.dat"};
	file<<"# h xexpl\n";
	for (uint64_t i=0;i<M;i++){
		double h=ha*std::exp((double)i/M*std::log(hb/ha));
		xexpl[i]=explode_at(0.5,0,a,b,h,0.5,2);
		file<<h<<' '<<xexpl[i]<<'\n';
	}
	file.close();

	return 0;
}