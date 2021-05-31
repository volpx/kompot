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
// #include <static_math/static_math.h>

double chi_p(const double a, const double r)
{
	return std::exp(-a * r * r);
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
	return 2 * std::pow(M_PI, 5.0 / 2) / ((a_p + a_r) * (a_q + a_s) * std::sqrt(a_p + a_q + a_r + a_s));
}
// <pq|0|rs>
double Zero_pqrs(const double, const double,const double,const double)
{
	return 0;
}

int main()
{
	std::cout << "Hartree-Fock He atom" << std::endl;

	// Steps
	constexpr uint64_t M=1e5;
	// Mixing
	constexpr double alpha=1e-2;
	// Basis
	constexpr uint32_t Nbase = 4;
	constexpr double a[Nbase] = {14.899983, 2.726485, 0.757447, 0.251390};

	// Number of orbitals in the Slater determinant
	constexpr uint32_t Z=2;
	constexpr uint32_t N=2;

	// Case for Nup=Ndown=N/2
	constexpr uint32_t Norbitals=N/2;

	constexpr double epsilon = 1e-7;
	double E_eigen = 1.0;
	double E_expected = 0.0;

	// Allocate workspace for computation of eigenvalues problems
	gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(Nbase);
	// Overlap matrix
	gsl_matrix *S = gsl_matrix_calloc(Nbase, Nbase);
	gsl_vector *SVal = gsl_vector_calloc(Nbase);
	gsl_matrix *SVec = gsl_matrix_calloc(Nbase, Nbase);

	// X matrix of coordinate change
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

	gsl_matrix *Pa = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Pb = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Pold = Pa;
	gsl_matrix *Pnew = Pb;
	gsl_matrix *Ptmp;

	// Eigenvalues of HF
	gsl_vector *E = gsl_vector_calloc(Nbase);
	
	// Calculation support matrixes
	gsl_matrix *A = gsl_matrix_calloc(Nbase, Nbase);
	gsl_matrix *Fp = gsl_matrix_calloc(Nbase, Nbase);

	// Compute constant matrixes
	for (uint32_t p = 0; p < Nbase; p++)
	{
		for (uint32_t q = 0; q < Nbase; q++)
		{
			gsl_matrix_set(S, p, q, S_pq(a[p], a[q]));
			gsl_matrix_set(Hcore, p, q, T_pq(a[p], a[q]) - Z * Coulomb_pq(a[p], a[q]));
		}
	}

	print_matrix(std::cout,Hcore);

	// Diagonalize the S matrix and build the X matrix
	gsl_eigen_symmv(S,SVal,SVec,ws);
	gsl_eigen_symmv_sort(SVal,SVec,GSL_EIGEN_SORT_VAL_ASC);
	
	for (uint64_t i = 0; i < Nbase; i++)
	{
		for (uint64_t j = 0; j < Nbase; j++)
		{
			// Scale the eigenvectors by their eigenvalue
			gsl_matrix_set(X, i, j, gsl_matrix_get(SVec, i, j) / std::sqrt(gsl_vector_get(SVal, j)));
		}
	}

	// Main loop
	uint64_t i=0;
	while(std::abs(E_eigen-E_expected)>epsilon)
	{
		// Compute F matrix
		gsl_matrix_memcpy(F, Hcore);
		HF_Fock_add_interaction(F, Norbitals, a, Pold, Coulomb_pqrs);

		// Solve Roothaan and obtain Cnew
		HF_solve_Roothaan(
			Cnew, E,
			F, X, ws,
			A, Fp);

		// Normalize the functions
		// for (uint32_t k = 0; k < Nbase; k++)
		// {
		// 	// compute the norm of |phi_k>
		// 	double norm2 = 0;
		// 	for (uint32_t p = 0; p < Nbase; p++)
		// 	{
		// 		for (uint32_t q = 0; q < Nbase; q++)
		// 		{
		// 			norm2+=gsl_matrix_get(Cnew,p,k)*gsl_matrix_get(Cnew,q,k)*gsl_matrix_get(S,p,q);
		// 		}
		// 	}
		// 	// scale the column k of Cnew for a normalized |phi_k>
		// 	for (uint32_t q = 0; q < Nbase; q++){
		// 		double *Cnew_qk=gsl_matrix_ptr(Cnew,q,k);
		// 		(*Cnew_qk)/=std::sqrt(norm2);
		// 	}
		// }

		// Compute the density matrix
		for (uint32_t p = 0; p < Nbase; p++)
		{
			for (uint32_t q = 0; q < Nbase; q++)
			{
				double *P_pq = gsl_matrix_ptr(Pnew, p, q);
				(*P_pq)=0;
				for (uint32_t k = 0; k < Norbitals; k++)
				{
					(*P_pq) += gsl_matrix_get(Cnew, p, k) * gsl_matrix_get(Cnew, q, k);
				}
			}
		}

		// Compute the energy eigen
		E_eigen=HF_Eeigen(
			Pnew,E,Nbase,
			a,Norbitals,
			Coulomb_pqrs);
		// E_eigen=gsl_vector_get(E,0);

		// Compute the energy expected
		E_expected = MF_Eexpected(Z, Pnew, a, Norbitals, Hcore, Coulomb_pqrs);


		// Mixing Cnew = (1-alpha)*Cold + alpha*Cnew
		// gsl_matrix_scale(Cold, 1 - alpha);
		// gsl_matrix_scale(Cnew, alpha);
		// gsl_matrix_add(Cnew, Cold);

		// Mixing Pnew = (1-alpha)*Pold + alpha*Pnew
		for (uint32_t p = 0; p < Nbase; p++)
		{
			for (uint32_t q = 0; q < Nbase; q++)
			{
				double *P_pq = gsl_matrix_ptr(Pnew, p, q);
				(*P_pq) = (1-alpha)* gsl_matrix_get(Pold, p, q) + alpha*(*P_pq);
			}
		}

		// Exchange C matrixes
		Ctmp = Cold;
		Cold = Cnew;
		Cnew = Ctmp;
		// Exchange P matrixes
		Ptmp = Pold;
		Pold = Pnew;
		Pnew = Ptmp;

		if (i%1==0){
			std::cout<<i<<' '<< E_eigen<< ' '<< E_expected<<'\n';
			// print_vector(std::cout,E);
		}
		i++;
	}
	std::cout << "Final E_eigen: "<<E_eigen << std::endl;
	std::cout<<"C: "<<std::endl;
	print_matrix(std::cout, Cold);
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
	gsl_matrix_free(Pa);
	gsl_matrix_free(Pb);

	return 0;
}