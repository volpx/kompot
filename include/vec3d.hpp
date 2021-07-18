#pragma once

#include <cstdint>

#include <cmath>
#include <cstdlib>
#include <ostream>

class Vec3D
{
public:
	// Empty constructor is -well- empty
	Vec3D(){};
	// Normal constructor (could overload the above)
	Vec3D(double x, double y, double z);
	// Copying constructor
	Vec3D(const Vec3D &other);
	// Moving constructor
	Vec3D(Vec3D &&other);

	// Actual class parameters
	double x;
	double y;
	double z;

	// Getters for the norm and norm^2
	double norm() const;
	double norm2() const;

	// Make a versor by either modify myself or create a copy
	Vec3D &normalize();
	Vec3D getNormalized() const;

	// Set to null vector
	Vec3D &clear();

	// Overload of operators to work on some abstraction
	Vec3D &operator=(const Vec3D &a);

	Vec3D &operator+=(const Vec3D &a);
	Vec3D &operator-=(const Vec3D &a);
	Vec3D &operator*=(const double k);
	Vec3D &operator/=(const double k);

	friend Vec3D operator-(const Vec3D &a, const Vec3D &b);
	friend Vec3D operator+(const Vec3D &a, const Vec3D &b);
	friend Vec3D operator*(const double k, const Vec3D &a);

	friend Vec3D &operator-(Vec3D &&a, const Vec3D &b);
	friend Vec3D &operator+(Vec3D &&a, const Vec3D &b);
	friend Vec3D &operator*(const double k, Vec3D &&a);

	// Overload of the scalar product operator
	// TODO: remove me
	friend double operator*(const Vec3D &a, const Vec3D &b);

	// Overload handy output operator
	friend std::ostream &operator<<(std::ostream &out, const Vec3D &vec);
};

