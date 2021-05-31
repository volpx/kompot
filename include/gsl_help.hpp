#pragma once

#include <cstdint>
#include <gsl/gsl_matrix.h>
#include <ostream>


void print_matrix(std::ostream& out, const gsl_matrix* m);
void print_vector(std::ostream& out, const gsl_vector* v);