#include "Tasks.h"
#include <iostream>
#include <iomanip>

double eps = 0.0000001;

Vector si(Matrix a, Vector b, Vector x0)
{
	std::cout << "Simple Iteration" << std::endl;

	Matrix A = Matrix(a.n, a.n + 1);

	for (int i = 0; i < A.m; i++)
	{
		for (int j = 0; j < A.m; j++)
			A.ar[i][j] = a.ar[i][j];
		A.ar[i][A.m] = b.ar[i];
	}

	for (int i = 0; i < A.m; i++)
	{
		double lead = A.ar[i][i];
		for (int j = 0; j < A.n; j++)
			A.ar[i][j] /= lead;
	}

	Matrix B = Matrix(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			B.ar[i][j] = A.ar[i][j] * (i != j);

	B.print();

	std::cout << "Norm: " << B.norm1() << std::endl;
	std::cout << std::endl;

	bool fin = false;
	int k = 0;
	Vector x = x0;
	std::cout << 0 << std::scientific << std::setprecision(8) << ' ' << x.ar[0] << ' ' << x.ar[1] << ' ' << x.ar[2] << std::endl;
	while (!fin)
	{
		Vector old = x;
		for (int i = 0; i < x.n; i++)
		{
			x.ar[i] = A.ar[i][A.m];
			for (int j = 0; j < i; j++)
				x.ar[i] -= A.ar[i][j] * old.ar[j];
			for (int j = i + 1; j < x.n; j++)
				x.ar[i] -= A.ar[i][j] * old.ar[j];
		}
		k++;
		std::cout << k << std::scientific << std::setprecision(8) << ' ' << x.ar[0] << ' ' << x.ar[1] << ' ' << x.ar[2] << std::endl;
		fin = (abs(x.ar[0] - old.ar[0]) < eps) && (abs(x.ar[1] - old.ar[1]) < eps) && (abs(x.ar[2] - old.ar[2]) < eps);
	}

	(a*x - b).print();
	return Vector(x);
}

Vector nm(Matrix a, Vector b, Vector x0)
{
	std::cout << "Nekrasov's method" << std::endl;

	Matrix A = Matrix(a.n, a.n + 1);

	for (int i = 0; i < A.m; i++)
	{
		for (int j = 0; j < A.m; j++)
			A.ar[i][j] = a.ar[i][j];
		A.ar[i][A.m] = b.ar[i];
	}

	for (int i = 0; i < A.m; i++)
	{
		double lead = A.ar[i][i];
		for (int j = 0; j < A.n; j++)
			A.ar[i][j] /= lead;
	}
	
	bool fin = false;
	int k = 0;
	Vector x = x0;
	std::cout << 0 << std::scientific << std::setprecision(8) << ' ' << x.ar[0] << ' ' << x.ar[1] << ' ' << x.ar[2] << std::endl;
	while (!fin)
	{
		Vector old = x;
		for (int i = 0; i < x.n; i++)
		{
			x.ar[i] = A.ar[i][A.m];
			for (int j = 0; j < i; j++)
				x.ar[i] -= A.ar[i][j] * x.ar[j];
			for (int j = i+1; j < x.n; j++)
				x.ar[i] -= A.ar[i][j] * old.ar[j];
		}
		k++;
		std::cout << k << std::scientific << std::setprecision(8) << ' ' << x.ar[0] << ' ' << x.ar[1] << ' ' << x.ar[2] << std::endl;
		fin = (abs(x.ar[0] - old.ar[0]) < eps) && (abs(x.ar[1] - old.ar[1]) < eps) && (abs(x.ar[2] - old.ar[2]) < eps);
	}

	(a*x - b).print();
	return Vector(x);
}

void Task2()
{
	Matrix A = Matrix(3, 3);
	A.ar[0][0] = 6.6872;
	A.ar[0][1] = 0.80267;
	A.ar[0][2] = -2.0646;
	A.ar[1][0] = 0.80267;
	A.ar[1][1] = 5.0782;
	A.ar[1][2] = 0.48036;
	A.ar[2][0] = -2.0646;
	A.ar[2][1] = 0.48036;
	A.ar[2][2] = 4.0293;

	Vector b = Vector(3);
	b.ar[0] = 10.029;
	b.ar[1] = -0.28161;
	b.ar[2] = -10.0159;

	Vector b1 = Vector(3);
	b1.ar[0] = 1.71509;
	b1.ar[1] = -4.484721;
	b1.ar[2] = 0.09906;

	A.print();
	b1.print();

	Vector x0 = Vector(3);
	x0.ar[0] = 3;
	x0.ar[1] = 4;
	x0.ar[2] = 5;

	si(A, b1, x0).print();

	nm(A, b1, x0).print();

}