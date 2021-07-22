#pragma once

#include <cstdint>
#include <cstddef>
#include <functional>

double integrator_simpson_cubic(
	const double y[],
	const size_t M,
	const double h);

double integrator_simpson_cubic(
	std::function<double(size_t)> y,
	const size_t M,
	const double h);

