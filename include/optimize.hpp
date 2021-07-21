#pragma once

#include <cstdint>
#include <cstddef>
#include <iostream>

extern "C"
{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
}

/* Fit the ydata as a function of xdata to an exponential:
 * y = l0 * exp(+l1 * x) + l2 
 *
 */
int fit_to_pexp(double lambda[3],
			   const double datay[],
			   const size_t N,
			   const double datax[] = nullptr);

/* Fit the ydata as a function of xdata to an exponential:
 * y = l0 * exp(-l1 * x) + l2 
 *
 */
int fit_to_nexp(double lambda[3],
			   const double datay[],
			   const size_t N,
			   const double datax[] = nullptr);

