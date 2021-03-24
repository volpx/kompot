#include "vector_help.hpp"

void map(std::vector<double> &vec, const std::function<double(uint64_t)> f)
{
	const uint64_t M{static_cast<uint64_t>(vec.size())};
	for (uint64_t i{0}; i < M; ++i)
	{
		vec[i] = f(i);
	}
}
void arange(std::vector<double> &vec, const double start, const double step)
{
	std::function<double(uint64_t)> f{
		[&](uint64_t i) -> double { return start + i * step; }};
	map(vec, f);
}
void linspace(std::vector<double> &vec, const double xmin, const double xmax)
{
	const double h{(xmax - xmin) / (vec.size() - 1)};
	arange(vec, xmin, h);
}
void fill(std::vector<double> &vec, const double val)
{
	std::function<double(uint64_t)> f{
		[&](uint64_t) -> double { return val; }};
	map(vec, f);
}

uint64_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map,
				 const uint64_t start,
				 uint64_t stop)
{

	uint64_t ind{start};
	double val{map ? map(vec[0]) : vec[0]};

	if (map != nullptr)
	{
		for (uint64_t i{0}; i < static_cast<uint64_t>(vec.size()); ++i)
		{
			if (map(vec[i]) < val)
			{
				val = map(vec[i]);
				ind = i;
			}
		}
	}
	else
	{
		for (uint64_t i{0}; i < static_cast<uint64_t>(vec.size()); ++i)
		{
			if (vec[i] < val)
			{
				val = vec[i];
				ind = i;
			}
		}
	}

	return ind;
}
