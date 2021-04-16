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

	constexpr uint64_t M=1e3;      // Number of energies
	constexpr double Ea=0;
	constexpr double Eb=0.5+7*1+0.25;

	constexpr uint64_t N=1e5; // number of mesh points
	constexpr double a=0.;   // start included
	constexpr double b=5;   // stop excluded


	constexpr double h=(b-a)/N; // Mesh spacing
	constexpr double h2=h*h;

	std::cout << "Problem constants:"
			  << "\nM: " << M
			  << "\nEa: " << Ea
			  << "\nEb: " << Eb
			  << "\nEh: " << (Eb-Ea)/(M-1)
			  << "\na: " << a
			  << "\nb: " << b
			  << "\nh: " << h
			  << "\nN: " << N
			  << std::endl;

	std::vector<double> vE(M);
	linspace(vE,Ea,Eb);

	std::vector<double> vy_xmax_odd(M);
	std::vector<double> vy_xmax_even(M);

	std::vector<double> x(N);
	for (uint64_t i=0;i<N;i++){
		x[i]=a+i*h;
	}
	std::vector<double> y_even(N);
	std::vector<double> y_odd(N);
	// std::vector<double> k(N);


	// Numerov algorithm
	for(uint64_t j = 0; j<M ; j++){

		if (j%100==0){
			std::cout<<"P: " << 1.0*j/M <<std::endl; 
		}

		// Initialization
		// y[0]=y_true(x[0], 0);
		// y[1]=y_true(x[1], 0);
		y_even[0]=1;
		y_even[1]=1*std::exp(-0.5*x[1]*x[1]);
		
		y_odd[0]=0;
		y_odd[1]=x[1];

		double E=vE[j];
		double k_ii2 = 2.*(E - V(x[0]));
		double k_i2 = 2.*(E - V(x[1]));
		double k2;

		for (uint64_t i = 2; i < N; i++)
		{
		
			// New k
			k2 = 2.*(E - V(x[i]));

			y_even[i] = (
				2.*y_even[i - 1] * (1. - 5. / 12 * h2 * k_i2) 
				- y_even[i - 2] * (1. + 1. / 12 * h2 * k_ii2 )) 
				/ (1. + h2 / 12.0 * k2 );

			y_odd[i] = (
				2.*y_odd[i - 1] * (1. - 5. / 12 * h2 * k_i2) 
				- y_odd[i - 2] * (1. + 1. / 12 * h2 * k_ii2 )) 
				/ (1. + h2 / 12.0 * k2 );
				
			// Move down the k for next iteration
			k_ii2 = k_i2;
			k_i2 = k2;
		}

		vy_xmax_even[j]=y_even[N-1];
		vy_xmax_odd[j]=y_odd[N-1];
	}

	std::ofstream file{"data/energy1.dat"};
	file << "#E y_xmax_even y_xmax_odd\n";
	for (uint64_t j = 0; j < M; j++){
		// Save the result
		file
			<< vE[j] << ' '
			<< vy_xmax_even[j] << ' '
			<< vy_xmax_odd[j]<<' '
			<< '\n';
	}

	// Print best estimate of Energy eigenvalues
	uint64_t ii=0;
	// uint64_t j=0;
	for (uint64_t i=0;i<M-1;i++){
		if (vy_xmax_even[i]*vy_xmax_even[i+1]<0){
			std::cout << "Eigen Even " << ii++ << ": " 
				<< vE[i + 1] - (vE[i + 1] - vE[i]) / (vy_xmax_even[i + 1] - vy_xmax_even[i]) * vy_xmax_even[i + 1] 
				<< '\n';
		}
		if(vy_xmax_odd[i]*vy_xmax_odd[i+1]<0){
			std::cout<<"Eigen Odd  "<< ii++ << ": "
				<<vE[i+1]-(vE[i+1]-vE[i])/(vy_xmax_odd[i+1]-vy_xmax_odd[i])*vy_xmax_odd[i+1]
				<<'\n';
		}
	}


	return 0;
}