#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <functional>


// vector manipulation
void map(double vec[],const size_t N, std::function<double(uint64_t)> f);
void map(double vec[],const size_t N, const std::vector<double> &x, const std::function<double(double)> f);
void arange(double vec[],const size_t N, const double start, const double step);
void linspace(double vec[],const size_t N, const double xmin, const double xmax);
void fill(double vec[],const size_t N, double val);

uint64_t ind_min(const std::vector<double> &vec,
				 std::function<double(double)> map = nullptr,
				 const uint64_t start=0,
				 uint64_t stop=1);


double min(
	const double vec[],
	const uint64_t N);
double average(
	const double vec[],
	const uint64_t N);

// Windowed array
class WArray{
public:
	WArray(){};
	WArray(const size_t N, const size_t a=0);
	WArray(double ptr[], const size_t N,const size_t a=0);
	// WArray(const WArray &other);
	// WArray(WArray &&other);

	bool owner_of_ptr=false;
	double *data=nullptr;
	size_t N=0;
	size_t a=0; // Start of window

	bool inWindow(size_t i);
	void set(const size_t i,const double val=0);
	double& get(const size_t i);
	double& operator[](const size_t i);
	void fill(const double val=0);

	~WArray();
};

// Circular array
class CArray{
public:
	CArray(){};
	CArray(const size_t N);
	CArray(double ptr[], const size_t N);

	bool owner_of_ptr=false;
	double *data=nullptr;
	size_t N=0;

	void set(const size_t i,const double val=0);
	double& get(const size_t i);
	double& operator[](const size_t i);
	void fill(const double val=0);

	~CArray();
};
