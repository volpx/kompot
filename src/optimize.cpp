#include "optimize.hpp"

int fit_to_nexp(double lambda[3],
			   const double datay[],
			   const size_t N,
			   const double datax[])
{
	// number of parameters
	const size_t p = 3;
	double lambda_init[p] = {lambda[0], lambda[1], lambda[2]};

	struct data
	{
		size_t N;
		const double *x;
		const double *y;
	};

	auto exp_f{
		[](const gsl_vector *lambda, void *data,
		   gsl_vector *f) -> int {
			size_t n = static_cast<struct data *>(data)->N;
			const double *x = static_cast<struct data *>(data)->x;
			const double *y = static_cast<struct data *>(data)->y;
			double l0 = gsl_vector_get(lambda, 0);
			double l1 = gsl_vector_get(lambda, 1);
			double l2 = gsl_vector_get(lambda, 2);
			if(x)
			{
				for (size_t i = 0; i < n; i++)
				{
					/* Model Yi = l0 * exp(-l1 * xi) + l2 */
					double Yi = l0 * std::exp(-l1 * x[i]) + l2;
					gsl_vector_set(f, i, Yi - y[i]);
				}
			}
			else
			{
				for (size_t i = 0; i < n; i++)
				{
					/* Model Yi = l0 * exp(-l1 * t_i) + l2 */
					double Yi = l0 * std::exp(-l1 * i) + l2;
					gsl_vector_set(f, i, Yi - y[i]);
				}
			}
			return GSL_SUCCESS;
		}};

	auto exp_df{
		[](const gsl_vector *lambda, void *data,
		   gsl_matrix *J) -> int {
			size_t n = static_cast<struct data *>(data)->N;
			const double *x = static_cast<struct data *>(data)->x;
			double l0 = gsl_vector_get(lambda, 0);
			double l1 = gsl_vector_get(lambda, 1);
			if(x)
			{
				for (size_t i = 0; i < n; i++)
				{
					/* Jacobian matrix J(i,j) = dfi / dxj, */
					/* where fi = (Yi - yi)/sigma[i],*/
					/* Yi = A * exp(-lambda * xi) + b */
					/* and the xj are the parameters (A,lambda,b) */
					double e = std::exp(-l1 * x[i]);
					gsl_matrix_set(J, i, 0, e);
					gsl_matrix_set(J, i, 1, -x[i] * l0 * e);
					gsl_matrix_set(J, i, 2, 1.0);
				}

			}
			else
			{
				for (size_t i = 0; i < n; i++)
				{
					/* Jacobian matrix J(i,j) = dfi / dxj, */
					/* where fi = (Yi - yi)/sigma[i],*/
					/* Yi = A * exp(-lambda * t_i) + b */
					/* and the xj are the parameters (A,lambda,b) */
					double e = std::exp(-l1 * i);
					gsl_matrix_set(J, i, 0, e);
					gsl_matrix_set(J, i, 1, -1.0*i * l0 * e);
					gsl_matrix_set(J, i, 2, 1.0);
				}
			}

			return GSL_SUCCESS;
		}};

	const gsl_multifit_nlinear_type *multifittype = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters fdf_params =
		gsl_multifit_nlinear_default_parameters();
	fdf_params.trs = gsl_multifit_nlinear_trs_lm;
	/* allocate workspace with default parameters */
	gsl_multifit_nlinear_workspace *w =
		gsl_multifit_nlinear_alloc(multifittype, &fdf_params, N, p);
	if (w == NULL)
	{
		std::cerr << "gsl_multifit_nlinear_alloc: " << w << std::endl;
		return -1;
	}
	gsl_vector_view lambda_init_wiew = gsl_vector_view_array(lambda_init, p);
	gsl_multifit_nlinear_fdf fdf;
	// define the function to be minimized
	fdf.f = exp_f;
	// set to NULL for finite-difference Jacobian
	fdf.df = exp_df;
	// not using geodesic acceleration
	fdf.fvv = NULL;
	fdf.n = N;
	fdf.p = p;
	struct data d = {N, datax, datay};
	fdf.params = &d;
	gsl_multifit_nlinear_init(
		&lambda_init_wiew.vector,
		&fdf,
		w);

	int status, info;

	/* solve the system with a maximum of 100 iterations */
	const double xtol = 1e-8;
	const double gtol = 1e-8;
	const double ftol = 0.0;
	status = gsl_multifit_nlinear_driver(
		100, xtol, gtol, ftol,
		NULL, NULL, &info, w);

	// std::cout
	// 	<< "status: " << status << '\n'
	// 	<< "niter:" << gsl_multifit_nlinear_niter(w) << '\n'
	// 	<< "reason stop: "
	// 	<< ((info == 1) ? "small step size" : "small gradient") << '\n'
	// 	<< "l0: " << gsl_vector_get(w->x, 0) << '\n'
	// 	<< "l1: " << gsl_vector_get(w->x, 1) << '\n'
	// 	<< "l2: " << gsl_vector_get(w->x, 2) << '\n'
	// 	<< std::endl;
	lambda[0] = gsl_vector_get(w->x, 0);
	lambda[1] = gsl_vector_get(w->x, 1);
	lambda[2] = gsl_vector_get(w->x, 2);

	gsl_multifit_nlinear_free(w);

	return 0;
}
