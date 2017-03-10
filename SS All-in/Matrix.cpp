#include "Matrix.h"
#include <iostream>
#include <iomanip>


Matrix::Matrix(int q, int w)
{
	m = q; n = w;
	ar = new double *[m];
	for (int i = 0; i < m; i++)
	{
		ar[i] = new double[n];
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < m; i++)
	{
		delete[] ar[i];
	}
	delete[] ar;
}

Matrix::Matrix(Matrix& other)
{
	m = other.m;
	n = other.n;
	ar = new double *[m];
	for (int i = 0; i < m; i++)
	{
		ar[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			ar[i][j] = other.ar[i][j];
		}
	}
}

Matrix Matrix::operator+(Matrix &other)
{
	Matrix res = Matrix(m, n);
	if (m != other.m || n != other.n)
	{
		std::cout << "Incorrect matrix sizes!" << std::endl;
		for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			res.ar[i][j] = 0;
		}
	}
	else
	{
		for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			res.ar[i][j] = ar[i][j] + other.ar[i][j];
		}
	}
	return Matrix(res);
}

Matrix Matrix::operator*(Matrix &other)
{
	Matrix res = Matrix(m, other.n);
	if (n != other.m)
	{
		std::cout << "Incorrect matrix sizes!" << std::endl;
		for (int i = 0; i < m; i++)
		for (int j = 0; j < other.n; j++)
		{
			res.ar[i][j] = 0;
		}
	}
	else
	{
		for (int i = 0; i < m; i++)
		for (int j = 0; j < other.n; j++)
		{
			res.ar[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				res.ar[i][j] += ar[i][k] * other.ar[k][j];
			}
		}
	}
	return Matrix(res);
}

Matrix Matrix::operator*(double &k)
{
	Matrix res = Matrix(m, n);
	for (int i = 0; i < m; i++)
	for (int j = 0; j < n; j++)
		res.ar[i][j] = ar[i][j] * k;
	return Matrix(res);
}

Vector Matrix::operator*(Vector &vec)
{
	Vector res = Vector(m);
	if (n != vec.n)
	{
		std::cout << "Matrix and vector sizes do not match!" << std::endl;
		for (int j = 0; j < n; j++)
		{
			res.ar[j] = 0;
		}
	}
	else
	{
		for (int i = 0; i < m; i++)
		{
			res.ar[i] = 0;
			for (int j = 0; j < n; j++)
			{
				res.ar[i] += ar[i][j] * vec.ar[j];
			}
		}
	}
	return Vector(res);
}

Matrix Matrix::operator-(Matrix &other)
{
	Matrix res = Matrix(m, n);
	if (m != other.m || n != other.n)
	{
		std::cout << "Incorrect matrix sizes!" << std::endl;
		for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			res.ar[i][j] = 0;
		}
	}
	else
	{
		for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			res.ar[i][j] = ar[i][j] - other.ar[i][j];
		}
	}
	return Matrix(res);
}

Matrix &Matrix::operator=(Matrix &other)
{
	if (ar != nullptr)
	{
		for (int i = 0; i < m; i++)
			delete[] ar[i];
		delete[] ar;
	}
	m = other.m; 
	n = other.n;
	ar = nullptr;
	ar = new double *[m];
	for (int i = 0; i < m; i++)
	{
		ar[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			ar[i][j] = other.ar[i][j];
		}
	}
	return *this;
}

Matrix Matrix::transpose()
{
	Matrix res = Matrix(n, m);
	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
	{
		res.ar[i][j] = ar[j][i];
	}
	return Matrix(res);
};

double Matrix::norm1()
{
	double res = 0;
	for (int j = 0; j < n; j++)
	{
		double tres = 0;
		for (int i = 0; i < m; i++)
		{
			tres += abs(ar[i][j]);
		}
		if (tres > res) res = tres;
	}
	return res;
}

double Matrix::norminf()
{
	double res = 0;
	for (int i = 0; i < m; i++)
	{
		double tres = 0;
		for (int j = 0; j < n; j++)
		{
			tres += abs(ar[i][j]);
		}
		if (tres > res) res = tres;
	}
	return res;
}

void Matrix::init()
{
	if (ar != nullptr)
	{
		for (int i = 0; i < m; i++)
			delete[] ar[i];
		delete[] ar;
	}
	int q, w;
	std::cout << "Matrix size: ";
	std::cin >> q;
	std::cout << " x ";
	std::cin >> w;
	std::cout << std::endl;
	m = q; n = w;
	ar = new double *[m];
	for (int i = 0; i < m; i++)
	{
		ar[i] = new double[m];
		for (int j = 0; j < n; j++)
			std::cin >> ar[i][j];
		std::cout << std::endl;
	}
}

void Matrix::print()
{
	std::cout << std::setprecision(8) << std::endl;
	std::cout << "Matrix: " << std::endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << std::scientific << ar[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Matrix::swapC(int a, int b) 
{
	double temp;
	for (int i = 0; i < m; i++)
	{
		temp = ar[i][a];
		ar[i][a] = ar[i][b];
		ar[i][b] = temp;
	}
}