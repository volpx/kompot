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
double Coulomb_pqts(const double a_p, const double a_q,const double a_r,const double a_s)
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
	std::cout << "Hartree-Fock Be atom 8 gaussians contracted to 2 functions" << std::endl;

	// Steps
	// constexpr uint64_t M=1e5;
	// Mixing
	constexpr double alpha=1e-2;
	// Basis
	constexpr uint32_t Nbase = 2;
	constexpr uint32_t Nsubbase = 4;
	constexpr double a[2 * Nsubbase] = {0.7064859542e2, 0.1292782254e2, 0.3591490662e1, 0.1191983464e1, 0.3072833610e1, 0.6652025433, 0.2162825386, 0.8306680972e-1};
	constexpr double b[2 * Nsubbase] = {0.5675242080e-1, 0.2601413550, 0.5328461143, 0.2916254405, -0.6220714565e-1, 0.2976804596e-4, 0.5588549221, 0.4977673218};

	// Number of orbitals in the Slater determinant
	constexpr uint32_t Z=4;
	constexpr uint32_t N=4;

	// Case for Nup=Ndown=N/2
	constexpr uint32_t Norbitals=N/2;

	constexpr double epsilon = 1e-8;
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

	// V2 2 body "matrix"
	std::vector<double> V2(Nbase*Nbase*Nbase*Nbase);
	
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
			double tmp1,tmp2;
			tmp1=tmp2=0;
			for (uint32_t i = 0; i < Nsubbase; i++)
			{
				for (uint32_t j = 0; j < Nsubbase; j++)
				{
					// Sum for S
					tmp1+=b[p*Nsubbase+i]*b[q*Nsubbase+j]*S_pq(a[p*Nsubbase+i],a[q*Nsubbase+j]);

					// Sum for Hcore
					tmp2+=b[p*Nsubbase+i]*b[q*Nsubbase+j]*
					(T_pq(a[p*Nsubbase+i], a[q*Nsubbase+j]) - Z * Coulomb_pq(a[p*Nsubbase+i], a[q*Nsubbase+j]));
				}
			}
			// Set matrix element <p|q>
			gsl_matrix_set(S, p, q, tmp1);
			// Set matrix element <p|Hcore|q>
			gsl_matrix_set(Hcore, p, q, tmp2);

			for (uint32_t t = 0; t < Nbase; t++)
			{
				for (uint32_t s = 0; s < Nbase; s++)
				{
					// Compute matrix element V2[p,q,t,s]=<pq|V2|ts>
					double tmp=0;
					for (uint32_t i = 0; i < Nsubbase; i++)
					{
						for (uint32_t j = 0; j < Nsubbase; j++)
						{
							for (uint32_t n = 0; n < Nsubbase; n++)
							{
								for (uint32_t m = 0; m < Nsubbase; m++)
								{
									// b_pi*b_qj
									tmp+=b[p*Nsubbase+i]*b[q*Nsubbase+j]*b[t*Nsubbase+n]*b[s*Nsubbase+m]
										*Coulomb_pqts(a[p*Nsubbase+i],a[q*Nsubbase+j],a[t*Nsubbase+n],a[s*Nsubbase+m]);
								}
							}
						}
					}
					V2[p*(Nbase*Nbase*Nbase)+q*(Nbase*Nbase)+t*(Nbase)+s]=tmp;
				}
			}
		}
	}

	// print_matrix(std::cout,Hcore);

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
		HF_Fock_add_interaction(F, Pold, V2);

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
			Norbitals, V2);
		// E_eigen=gsl_vector_get(E,0);

		// Compute the energy expected
		E_expected = MF_Eexpected(Z, Pnew, a, Norbitals, Hcore, V2);


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

		// if (i%1==0){
		// 	std::cout<<i<<' '<< E_eigen<< ' '<< E_expected<<'\n';
		// 	// print_vector(std::cout,E);
		// }
		i++;
	}
	std::cout<<"Iterations: "<<i<<std::endl;
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