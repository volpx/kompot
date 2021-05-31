#include "HF.hpp"

int HF_Fock_add_interaction(
	gsl_matrix *F,
	const gsl_matrix *P,
	const std::vector<double> V2)
{
	const uint32_t Nbase=F->size1;

	for (uint32_t p = 0; p < Nbase; p++)
	{
		// row p
		for (uint32_t q = 0; q < Nbase; q++)
		{
			// column q

			// get pointer to matrix element p,q
			double *element=gsl_matrix_ptr(F, p, q);

			double tmp1=0;
			for (uint32_t t = 0; t < Nbase; t++)
			{
				for (uint32_t s = 0; s < Nbase; s++)
				{
					// tmp1+=gsl_matrix_get(P,s,t)*
					// 	(2*V_pqrs(a[p],a[s],a[q],a[t])
					// 	-V_pqrs(a[p],a[s],a[t],a[q]));
					tmp1+=gsl_matrix_get(P,s,t)*
						(2*V2[p*(Nbase*Nbase*Nbase)+s*(Nbase*Nbase)+q*(Nbase)+t]
						  -V2[p*(Nbase*Nbase*Nbase)+s*(Nbase*Nbase)+t*(Nbase)+q]);
				}
			}

			// Add to the matrix element
			(*element)+=tmp1;
		}
	}

	return 0;
}

int HF_solve_Roothaan(
	gsl_matrix *C,
	gsl_vector *E,
	const gsl_matrix *F,
	const gsl_matrix *V,
	gsl_eigen_symmv_workspace *ws,
	gsl_matrix *A,
	gsl_matrix *Fp)
{
	gsl_matrix_set_zero(A);
	gsl_matrix_set_zero(Fp);

	// Compute Fp=F'=V^T F V
	gsl_blas_dgemm(
		CblasTrans, CblasNoTrans,
		1.0, V, F,
		0.0, A);
	gsl_blas_dgemm(
		CblasNoTrans, CblasNoTrans,
		1.0, A, V,
		0.0, Fp);

	// Solve eigenvalue problem for Fp
	gsl_eigen_symmv(Fp, E, A, ws);
	gsl_eigen_symmv_sort(E,A,GSL_EIGEN_SORT_VAL_ASC);

	// C=V C'
	gsl_blas_dgemm(
		CblasNoTrans, CblasNoTrans,
		1.0, V, A,
		0.0, C);
	
	return 0;
}

double HF_Eeigen(
	const gsl_matrix *P,
	const gsl_vector *E,
	const uint32_t Nbase,
	const uint32_t Norbitals,
	const std::vector<double> V2)
{
	double Etot = 0;

	for (uint32_t k = 0; k < Norbitals; k++)
	{
		Etot += 2 * gsl_vector_get(E, k);
	}

	for (uint32_t p = 0; p < Nbase; p++)
	{
		for (uint32_t q = 0; q < Nbase; q++)
		{
			for (uint32_t t = 0; t < Nbase; t++)
			{
				for (uint32_t s = 0; s < Nbase; s++)
				{
					// Etot += -gsl_matrix_get(P,p,q) *gsl_matrix_get(P,t,s)* 
					// (2 * V_pqrs(a[p], a[t], a[q], a[s]) - V_pqrs(a[p], a[t], a[s], a[q]));

					Etot += -gsl_matrix_get(P,p,q) *gsl_matrix_get(P,t,s)* 
					(2 * V2[p*(Nbase*Nbase*Nbase)+t*(Nbase*Nbase)+q*(Nbase)+s] 
					    -V2[p*(Nbase*Nbase*Nbase)+t*(Nbase*Nbase)+s*(Nbase)+q]);
				}
			}
		}
	}

	return Etot;
}
