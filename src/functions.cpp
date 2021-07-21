#include "functions.hpp"

double randn()
{
	// Statically save the two values
	constexpr double two_pi{M_PI * 2};
	constexpr double scale{1.0 / 32767};

	// If generated==true it already contains the next value
	thread_local double next;
	// Track if already generated
	thread_local bool generated;
	generated = !generated;
	if (generated)
	{
		// Generate new
		double q, w;
		// Don't want a zero, bad for the logarithm
		do
		{
			q = (rand() % 32767) * scale;
		} while (q == 0);
		// Scale the value between [0,1)
		w = (rand() % 32767) * scale;
		// Save for next call
		next = std::sqrt(-2.0 * std::log(q)) * std::cos(two_pi * w);
		// Return value
		return std::sqrt(-2.0 * std::log(q)) * std::sin(two_pi * w);
	}
	else
	{
		// Already good from previous call
		return next;
	}
}

double randu()
{
	constexpr double scale{1.0 / 32767};
	return (rand() % 32767) * scale;
}

int number_of_significant_digits(double x, double corr)
{
	return (int)std::log10(x / corr);
}

void PressEnterToContinue(){
	int c;
	std::cout << "Press ENTER to continue... " << std::endl;
	do 
		c = getchar(); 
	while ((c != '\n') && (c != EOF));
}