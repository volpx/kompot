#include "differentiate.hpp"

// https://web.media.mit.edu/~crtaylor/calculator.html
int diff_2_5points_allmesh(
	const double y[], const uint64_t M, 
	double y_xx[], const double h)
{
	if (M < 4)
	{
		return -1;
	}
	const double h2 = h * h;
	y_xx[0] = (1 * y[0 + 0] - 2 * y[0 + 1] + 1 * y[0 + 2]) / (h2);
	y_xx[1] = (1 * y[1 - 1] - 2 * y[1 + 0] + 1 * y[1 + 1] + 0 * y[1 + 2]) / (h2);
	for (uint64_t i = 2; i < M - 2; i++)
	{
		y_xx[i] = (-1 * y[i - 2] + 16 * y[i - 1] - 30 * y[i + 0] + 16 * y[i + 1] - 1 * y[i + 2]) / (12 * h2);
	}
	y_xx[M - 2] = (0 * y[(M - 2) - 2] + 1 * y[(M - 2) - 1] - 2 * y[(M - 2) + 0] + 1 * y[(M - 2) + 1]) / (h2);
	y_xx[M - 1] = (1 * y[(M - 1) - 2] - 2 * y[(M - 1) - 1] + 1 * y[(M - 1) + 0]) / (h2);

	return 0;
}