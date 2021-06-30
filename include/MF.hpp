#pragma once

// Mean Field routines

#include <cstdint>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <functional>
#include <vector>

#include <iostream>

#define FILLED_NO 0
#define FILLED_HALFUP_HALFDOWN 1

double MF_Eexpected(const int Z,const gsl_matrix *P, const double a[], 
					const uint32_t Norbitals, const gsl_matrix *Hcore,
					const std::vector<double> V_pqts);