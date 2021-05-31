#include "uniconst.hpp"
#include "vector_help.hpp"
#include "HF.hpp"
#include "MF.hpp"
#include "gsl_help.hpp"

#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

// #include <eigen3/Eigen/Dense>
// #define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

double chi_p(const double a, const double r)
{
	return std::exp(-a * r*r);
}
// Overlap integral with gaussians
double S_pq(const double a_p, const double a_q)
{
	return std::pow(M_PI / (a_p + a_q), 3.0 / 2);
}
// Laplpacian integral with gaussians 
double T_pq(const double a_p, const double a_q)
{
	return 3.0 * a_p * a_q * std::pow(M_PI, 3.0 / 2) / std::pow(a_p + a_q, 5.0 / 2);
}
// 1/r integral with gaussians
double Coulomb_pq(const double a_p, const double a_q)
{
	return 2.0 * M_PI /(a_p + a_q);
}
// <pq|1/(r12)|rs>
double Coulomb_pqrs(const double a_p, const double a_q,const double a_r,const double a_s)
{
	return Coulomb_pq(a_p, a_r) * Coulomb_pq(a_q, a_s);
}
// <pq|0|rs>
double Zero_pqrs(const double, const double,const double,const double)
{
	return 0;
}

int main()
{
	std::cout << "Hartree-Fock Hydrogen atom" << std::endl;

	// Steps
	// constexpr uint64_t N=1e4;
	// Mixing
	constexpr double alpha=1e-3;
	// Basis
	constexpr uint32_t Nbase = 4;
	constexpr double a[Nbase] = {13.00773, 1.962079, 0.444529, 0.1219492};

	// Number of orbitals in the Slater determinant
	constexpr uint32_t Norbitals=1;
	constexpr uint32_t Z=1;

	constexpr double epsilon = 1e-7;
	double E_eigen = 1.0;
	double E_expected = 0.0;

	// Allocate workspace for computation of eigenvalues problems
	gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(Nbase);
	// Overlap matrix
	gsl_matrix *S = gsl_matrix_calloc(Nbase, Nbase);
	gsl_vector *SVal = gsl_vector_calloc(Nbase);
	gsl_matrix *SVec = gsl_matrix_calloc(Nbase, Nbase);

	// X matrix
	gsl_matrix *X = gsl_matrix_calloc(Nbase, Nbase);

	// Hcore matrix
	gsl_matrix *Hcore = gsl_matrix_calloc(Nbase, Nbase);
	
	// Fock matrix
	gsl_matrix *F = gsl_matrix_calloc(Nbase, Nbase);
	
	// C = V*C' coefficient matrixes and handlers
	gsl_matrix *Ca = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Cb = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Cold = Ca;
	gsl_matrix *Cnew = Cb;
	gsl_matrix *Ctmp;

	// Eigenvalues of HF
	gsl_vector *E = gsl_vector_calloc(Nbase);
	
	// Calculation support matrixes
	gsl_matrix *A = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Fp = gsl_matrix_calloc(Nbase, Nbase);

	// Compute constant matrixes
	for (uint32_t p = 0; p < Nbase; p++)
	{
		for (uint32_t q = 0; q <Nbase; q++)
		{
			gsl_matrix_set(S, p, q, S_pq(a[p], a[q]));
			gsl_matrix_set(Hcore, p, q, T_pq(a[p], a[q]) - Z * Coulomb_pq(a[p], a[q]));
		}
	}

	// Diagonalize the S matrix and build the X matrix
	gsl_eigen_symmv(S,SVal,SVec,ws);
	gsl_eigen_symmv_sort(SVal,SVec,GSL_EIGEN_SORT_VAL_ASC);

	for (uint64_t i = 0; i < Nbase; i++)
	{
		for (uint64_t j = 0; j < Nbase; j++)
		{
			double U_ij = gsl_matrix_get(SVec,i,j);
			gsl_matrix_set(X, i, j, U_ij / std::sqrt(gsl_vector_get(SVal, j)));
		}
	}

	// Main loop
	uint64_t i=0;
	while(std::abs(E_eigen-E_expected)>epsilon)
	{
		// Compute F matrix
		gsl_matrix_memcpy(F, Hcore);
		// HF_Fock_add_interaction(F, a, Cold, Coulomb_pqrs);

		// Solve Roothaan and obtain Cnew
		HF_solve_Roothaan(
			Cnew, E,
			F, X, ws,
			A, Fp);

		// Mixing Cnew = (1-alpha)*Cold + alpha*Cnew
		gsl_matrix_scale(Cold, 1 - alpha);
		gsl_matrix_scale(Cnew, alpha);
		gsl_matrix_add(Cnew, Cold);

		// Compute the energy eigen
		E_eigen=HF_Eeigen(
			Cnew,E,Nbase,
			a,Norbitals,
			Zero_pqrs);

		// Compute the energy expected
		E_expected=MF_Eexpected(Cnew,a,Norbitals,T_pq,Coulomb_pq,Zero_pqrs);

		if (i%1000)
			std::cout<<i<<' '<< E_eigen<< ' '<< E_expected<<'\n';

		// Exchange C matrixes
		Ctmp = Cold;
		Cold = Cnew;
		Cnew = Ctmp;

		i++;
	}
	std::cout << "Final E_eigen: "<<E_eigen << std::endl;
	std::cout<<"C: "<<std::endl;
	print_matrix(std::cout, Cnew);
	std::cout<<"eigens:"<<std::endl;
	print_vector(std::cout, E);

	// Free the workspace
	gsl_eigen_symmv_free(ws);
	// Free matrixes space
	gsl_matrix_free(S);
	gsl_matrix_free(SVec);
	gsl_vector_free(SVal);

	gsl_matrix_free(X);
	gsl_matrix_free(A);
	gsl_matrix_free(Fp);
	gsl_vector_free(E);

	gsl_matrix_free(Hcore);
	gsl_matrix_free(F);

	gsl_matrix_free(Ca);
	gsl_matrix_free(Cb);

	// // Open main file of output
	// std::ofstream file{"data/numerov1.dat"};
	// file << "#x y_true y_numerov\n";
	// for (uint64_t i = 0; i < N; i++){
	// 	// Save the result
	// 	file
	// 		<< x[i] << ' '
	// 		<< y_true(x[i],n) << ' '
	// 		<< y[i] << ' '
	// 		<< '\n';
	// }

	return 0;
}
