#include "vec3d.hpp"

Vec3D::Vec3D(double x, double y, double z)
	: x{x},
	  y{y},
	  z{z}
{
}

Vec3D::Vec3D(const Vec3D &other)
	: x{other.x},
	  y{other.y},
	  z{other.z}
{
}

Vec3D::Vec3D(Vec3D &&other)
	: x{other.x},
	  y{other.y},
	  z{other.z}
{
}

double Vec3D::norm() const
{
	return std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

double Vec3D::norm2() const
{
	return this->x * this->x + this->y * this->y + this->z * this->z;
}

Vec3D &Vec3D::normalize()
{
	double n{this->norm()};
	this->x /= n;
	this->y /= n;
	this->z /= n;

	return *(this);
}

Vec3D Vec3D::getNormalized() const
{
	// Use copy constructor to create a copy instance
	Vec3D res{*this};
	res.normalize();
	return res;
}

Vec3D &Vec3D::clear()
{
	this->x = this->y = this->z = 0;
	return *this;
}

Vec3D &Vec3D::operator=(const Vec3D &a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
	return *(this);
}

Vec3D &Vec3D::operator+=(const Vec3D &a)
{
	this->x += a.x;
	this->y += a.y;
	this->z += a.z;
	return *(this);
}

Vec3D &Vec3D::operator-=(const Vec3D &a)
{
	this->x -= a.x;
	this->y -= a.y;
	this->z -= a.z;
	return *(this);
}

Vec3D &Vec3D::operator*=(const double k)
{
	this->x *= k;
	this->y *= k;
	this->z *= k;
	return *(this);
}

Vec3D &Vec3D::operator/=(const double k)
{
	this->x /= k;
	this->y /= k;
	this->z /= k;
	return *(this);
}

Vec3D operator-(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3D operator+(const Vec3D &a, const Vec3D &b)
{
	return Vec3D{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3D operator*(const double k, const Vec3D &a)
{
	return Vec3D{a.x * k, a.y * k, a.z * k};
}

Vec3D &operator-(Vec3D &&a, const Vec3D &b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}
Vec3D &operator+(Vec3D &&a, const Vec3D &b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
Vec3D &operator*(const double k, Vec3D &&a)
{
	a.x *= k;
	a.y *= k;
	a.z *= k;
	return a;
}

double operator*(const Vec3D &a, const Vec3D &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

std::ostream &operator<<(std::ostream &out, const Vec3D &vec)
{
	out << "[" << vec.x << "," << vec.y << "," << vec.z << "]";
	return out;
}
