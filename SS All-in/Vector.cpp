#include "Vector.h"
#include <iostream>

Vector::Vector(int q)
{
	n = q;
	ar = new double[n]; 
}

Vector::~Vector()
{
	delete[] ar;
}

Vector::Vector(Vector &other)
{
	n = other.n;
	ar = new double[n];
	for (int i = 0; i < n; i++)
	{
		ar[i] = other.ar[i];
	}
}

Vector Vector::operator+(Vector &other)
{
	Vector res = Vector(n);
	if (n != other.n)
	{
		std::cout << "Incorrect vector sizes!" << std::endl;
		for (int i = 0; i < n; i++)
			res.ar[i] = 0;
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			res.ar[i] = ar[i] = other.ar[i];
		}
	}
	return Vector(res);
}

double Vector::operator*(Vector &other)
{
	if (n != other.n)
	{
		std::cout << "Incorrect vector sizes!" << std::endl;
		return 0;
	}
	else
	{
		double result = 0;
		for (int i = 0; i < n; i++)
		{
			result += ar[i] * other.ar[i];
		}
		return result;
	}
}

Vector Vector::operator*(double k)
{
	Vector res = Vector(n);
	for (int i = 0; i < n; i++)
	{
		res.ar[i] = ar[i] * k;
	}
	return Vector(res);
}

Vector Vector::operator-(Vector &other)
{
	Vector res = Vector(n);
	if (n != other.n)
	{
		std::cout << "Incorrect vector sizes!" << std::endl;
		for (int i = 0; i < n; i++)
			res.ar[i] = 0;
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			res.ar[i] = ar[i] - other.ar[i];
		}
	}
	return Vector(res);
}

void Vector::print()
{
	std::cout << "Vector: " << std::endl;
	for (int i = 0; i < n; i++)
	{
		std::cout << ar[i] << std::endl;
	}
	std::cout << std::endl;
}

void Vector::init()
{
	int q;
	std::cout << "Vector length: ";
	std::cin >> q; 
	std::cout << std::endl;
	n = q;
	if (ar != nullptr) delete[] ar;
	ar = new double[n];
	for (int i = 0; i < n; i++)
	{
		std::cin >> ar[i];
	}
}

double Vector::norm1()
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += abs(ar[i]);
	}
	return res;
}

double Vector::norm2()
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += ar[i]*ar[i];
	}
	return sqrt(res);
}

double Vector::norminf()
{
	double res = abs(ar[0]);
	for (int i = 1; i < n; i++)
	{
		double q = abs(ar[i]);
		if (q>res) res = q;
	}
	return res;
}

Vector &Vector::operator=(Vector &other)
{
	n = other.n;
	if (ar!=nullptr) delete[] ar;
	ar = nullptr;
	ar = new double[n];
	for (int i = 0; i < n; i++)
	{
		ar[i] = other.ar[i];
	}
	return *this;
}