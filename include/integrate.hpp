#pragma once

#include <cstdint>

double integrator_simpson_cubic(
	const double y[],
	const uint64_t M,
	const double h);