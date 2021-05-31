#include "gsl_help.hpp"

void print_matrix(std::ostream& out, const gsl_matrix* m)
{
	for (uint64_t i = 0; i < m->size1; i++)
	{
		for (uint64_t j = 0; j< m->size2; j++)
		{
			out << ' ' << gsl_matrix_get(m,i,j);
		}
		out<< '\n';
	}
}

void print_vector(std::ostream& out, const gsl_vector* v)
{
	for (uint64_t i = 0; i < v->size; i++)
	{
		out << ' ' << gsl_vector_get(v,i);
	}
	out<< '\n';
}