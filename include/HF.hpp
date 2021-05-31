#pragma once

// Hartree Fock routines

#include <cstdint>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <functional>
#include <vector>

#include <iostream>

int HF_Fock_add_interaction(
	gsl_matrix *F,
	const gsl_matrix *C,
	std::vector<double> V2);

int HF_solve_Roothaan(
	// Outputs
	// C : Coefficients matrix
	// E : Eigenvalue matrix
	gsl_matrix *C, 
	gsl_vector *E, 
	// Inputs
	// F : Fock operator matrix
	// V : Eigenvectors/sqrt(lambda) of the overlap matrix
	const gsl_matrix *F, 
	const gsl_matrix *V,
	// Already allocated workspace etc...
	gsl_eigen_symmv_workspace *ws,
	gsl_matrix *A,
	gsl_matrix *Fp);

// 
double HF_Eeigen(
	const gsl_matrix *C,
	const gsl_vector *E,
	const uint32_t Nbase,
	const uint32_t Norbitals,
	std::vector<double> V2);