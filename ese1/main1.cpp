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

int main()
{
	std::cout << "Ciao ese1!" << std::endl;

	constexpr uint64_t N=1e5; // number of mesh points
	constexpr double a=5.;   // start included
	constexpr double b=-5.;   // stop excluded
	constexpr int n=1; // Energy quantum number

	constexpr double h=(b-a)/N; // Mesh spacing
	constexpr double E0=0.5+1.0*n;
	constexpr double E=E0;
	constexpr double h2=h*h;

	std::cout << "Problem constants:"
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << "\nE0: " << E0
			  << std::endl;

	std::vector<double> x(N);
	for (uint64_t i=0;i<N;i++){
		x[i]=a+i*h;
	}
	std::vector<double> y(N);
	// std::vector<double> k(N);

	// Initialization
	y[0]=y_true(x[0], n);
	y[1]=y_true(x[1], n);

	// Numerov algorithm
	double k_ii2 = 2.*(E - V(x[0]));
	double k_i2 = 2.*(E - V(x[1]));
	double k2;
	for (uint64_t i = 2; i < N; i++)
	{
		// if ((E - V(x[i]))<0){
		//     std::cout<<(V(x[i]))<<std::endl;
		// }
		// New k
		k2 = 2.*(E - V(x[i]));

		y[i] = (
			2.*y[i - 1] * (1. - 5. / 12 * h2 * k_i2) 
			- y[i - 2] * (1. + 1. / 12 * h2 * k_ii2 )) 
			/ (1. + h2 / 12.0 * k2 );

		// Move down the k for next iteration
		k_ii2 = k_i2;
		k_i2 = k2;
	}
	

	// Open main file of output
	std::ofstream file{"data/numerov1.dat"};
	file << "#x y_true y_numerov\n";
	for (uint64_t i = 0; i < N; i++){
		// Save the result
		file
			<< x[i] << ' '
			<< y_true(x[i],n) << ' '
			<< y[i] << ' '
			<< '\n';
	}


	return 0;
}