#include <MF.hpp>

double MF_Eexpected(const uint8_t filled, const gsl_matrix *P, 
					const double a[], const uint32_t Norbitals,
					const gsl_matrix *Hcore,
					const std::vector<double> V2)
{
	uint32_t Nbase = P->size1;
	double Etot = 0;
	double coeff=0;

	for (uint32_t p = 0; p < Nbase; p++)
	{
		for (uint32_t q = 0; q < Nbase; q++)
		{
			// One body part
			Etot += (2) * gsl_matrix_get(P, p, q) * gsl_matrix_get(Hcore, p, q);

			// Two body part
			double tmp = 0;
			for (uint32_t t = 0; t < Nbase; t++)
			{
				for (uint32_t s = 0; s < Nbase; s++)
				{
					// tmp += gsl_matrix_get(P,p,q) * gsl_matrix_get(P, t, s) * 
					// 	(2 * V_pqts(a[p], a[t], a[q], a[s]) 
					// 	- V_pqts(a[p], a[t], a[s], a[q]));
					
					tmp += gsl_matrix_get(P,p,q) * gsl_matrix_get(P, t, s) * 
						(2 * V2[p*(Nbase*Nbase*Nbase)+t*(Nbase*Nbase)+q*(Nbase)+s] 
						   - V2[p*(Nbase*Nbase*Nbase)+t*(Nbase*Nbase)+s*(Nbase)+q]);
				}
			}
			Etot += tmp;
		}
	}
	return Etot;
}