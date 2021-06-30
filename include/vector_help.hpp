#ifndef __VECTOR_HELP_HPP__
#define __VECTOR_HELP_HPP__

#include <cstdint>
#include <vector>
#include <functional>


// vector manipulation
void map(std::vector<double> &vec, std::function<double(uint64_t)> f);
void map(std::vector<double> &y, const std::vector<double> &x, const std::function<double(double)> f);
void arange(std::vector<double> &vec, const double start, const double step);
void linspace(std::vector<double> &vec, const double xmin, const double xmax);
void fill(std::vector<double> &vec, double val);

uint64_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr,
				 const uint64_t start=0,
				 uint64_t stop=1);

#endif // __VECTOR_HELP_HPP__