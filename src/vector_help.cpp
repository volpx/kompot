#include "vector_help.hpp"

size_t ind_min(
	const double vec[], const size_t N,
	const size_t start,
	std::function<double(double)> map)
{
	size_t idx_min = 0;
	double val_min=map(vec[idx_min]);
	for (size_t i = start; i < N; i++)
	{
		if (val_min > map(vec[i]))
		{
			idx_min = i;
		}
	}
	return idx_min;
}

double min(
	const double vec[],
	const size_t N)
{
	double min = vec[0];
	for (size_t i = 0; i < N; i++)
	{
		if (min > vec[i])
		{
			min = vec[i];
		}
	}
	return min;
}

double average(
	const double vec[],
	const size_t N)
{
	double mean = 0;
	for (size_t i = 0; i < N; i++)
	{
		mean += vec[i];
	}
	return mean / N;
}

double variance(
	const double vec[],
	const size_t N,
	const int ddof)
{
	double sum = 0,sum2=0;

	for (size_t i = 0; i < N; i++)
	{
		sum  += vec[i];
		sum2 += vec[i]*vec[i];
	}
	return (sum2-sum*sum/N) / (N-ddof);
}

double stddev(
	const double vec[],
	const size_t N,
	const int ddof)
{
	return std::sqrt(variance(vec,N,ddof));
}


// Windowed array
WArray::WArray(const size_t N, const size_t a)
	: owner_of_ptr{true},
	  data{new double[N]},
	  N{N},
	  a{a}
{
}

WArray::WArray(double ptr[], const size_t N,const size_t a)
	: owner_of_ptr{false},
	  data{ptr},
	  N{N},
	  a{a}
{
}

// WArray::WArray(const WArray &other)
// 	: owner_of_ptr{true},
// 	  data{new double[other.N]},
// 	  N{other.N},
// 	  a{other.a}
// {
// 	for(size_t i=0; i<N; i++)
// 	{
// 		this->data[i]=other.data[i];
// 	}
// }

// WArray::WArray(const WArray &other)
// 	: owner_of_ptr{true},
// 	  data{new double[other.N]},
// 	  N{other.N},
// 	  a{other.a}
// {
// 	for(size_t i=0; i<N; i++)
// 	{
// 		this->data[i]=other.data[i];
// 	}
// }

bool WArray::inWindow(const size_t i)
{
	if(i>=this->a && i<(this->a+this->N))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void WArray::set(const size_t i,const double val)
{
	if(this->inWindow(i))
	{
		this->data[i-a]=val;
	}
}

double &WArray::get(const size_t i)
{
	return this->data[i-a];
}
double& WArray::operator[](const size_t i)
{
	return this->get(i);
}

void WArray::fill(const double val)
{
	for (size_t i{0}; i < this -> N; i++)
		this->data[ i ]=val;
}

WArray::~WArray()
{
	if(this->owner_of_ptr)
	{
		delete[] this->data;
	}
}

// Circular array 
CArray::CArray(const size_t N)
	: owner_of_ptr{true},
	  data{new double[N]},
	  N{N}
{
}

CArray::CArray(double ptr[], const size_t N)
	: owner_of_ptr{false},
	  data{ptr},
	  N{N}
{
}

void CArray::set(const size_t i,const double val)
{
	this->data[ i % this->N]=val;
}

double &CArray::get(const size_t i)
{
	return this->data[i % this->N];
}
double& CArray::operator[](const size_t i)
{
	return this->get(i);
}

void CArray::fill(const double val)
{
	for (size_t i{0}; i < this -> N; i++)
		this->data[ i ]=val;
}

CArray::~CArray()
{
	if(this -> owner_of_ptr)
	{
		delete[] this->data;
	}
}
