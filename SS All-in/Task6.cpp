#include "Tasks.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>

double h, t, T;
int k;

double ustar(double y, double x)
{
	return std::exp(-4 * y)*std::sin(2 * x) + std::exp(-y)*(1 - x*x);
	//return x + y;
}

double a(double x)
{
	return std::exp(-x);
	//return x;
}

double b(double x)
{
	return std::exp(-4 * x) * std::sin(2);
	//return 1 + x;
}

double g(double x)
{
	return std::sin(2 * x) + 1 - x*x;
	//return x;
}

double f(double y, double x)
{
	return std::exp(-y)*(x*x - 3);
	//return 1;
}

void ETS()
{
	std::cout << "\nExplicit three-point scheme:\n";
	int N1 = T / t + 1, N2 = 1 / h + 1;
	Matrix A(N1, N2);
	
	for (int i = 0; i < N1; ++i)
	{
		A.ar[i][0] = a(t*i);
		A.ar[i][N2-1] = b(t*i);
	}
	for (int i = 0; i < N2; ++i)
	{
		A.ar[0][i] = g(i*h);
	}
	for (int i = 1; i < N1; ++i)
	{
		for (int j = 1; j < N2 - 1; ++j)
		{
			A.ar[i][j] = (1 - 2 * t / (h*h))*A.ar[i - 1][j] + (t / (h*h))*(A.ar[i - 1][j - 1] + A.ar[i - 1][j + 1]) + t*f((i - 1)*t, j*h);
		}
	}
	double delta = 0;
	std::cout << k << ' ' << k*t << '\n';
	for (int i = 0; i < N2; ++i)
	{
		//std::cout << std::scientific << A.ar[k][i]  << ' ';
		delta = std::max(delta, std::abs(ustar(k*t, i*h) - A.ar[k][i]));
	}
	std::cout << "\ndelta = " << delta << "\n\n";
}

void ITS()
{
	std::cout << "\nImplicit three-point scheme:\n";
	int N1 = T / t + 1, N2 = 1 / h + 1;
	Matrix A(N1, N2);
	Matrix B(N2, N2 + 1);
	for (int i = 0; i < N2; ++i)
	{
		A.ar[0][i] = g(i*h);
	}
	for (int i = 0; i < N2; ++i)
	{
		for (int j = 0; j < N2; ++j)
		{
			B.ar[i][j] = 0;
		}
	}
	B.ar[0][0] = 1;
	B.ar[N2 - 1][N2 - 1] = 1;
	for (int i = 1; i < N2 - 1; ++i)
	{
		B.ar[i][i - 1] = -t / (h*h);
		B.ar[i][i] = 1 + 2 * t / (h*h);
		B.ar[i][i + 1] = -t / (h*h);
	}
	Vector r;
	for (int i = 1; i < N1; ++i)
	{
		B.ar[0][N2] = a(i*t);
		for (int j = 1; j < N2-1; ++j)
		{
			B.ar[j][N2] = t*f(i*t, j*h) + A.ar[i - 1][j];
		}
		B.ar[N2 - 1][N2] = b(i*t);
		r = gaussC(B);
		for (int j = 0; j < N2; ++j)
		{
			A.ar[i][j] = r.ar[j];
		}
	}
	double delta = 0;
	std::cout << k << ' ' << k*t << '\n';
	for (int i = 0; i < N2; ++i)
	{
		//std::cout << std::scientific << A.ar[k][i] << ' ';
		delta = std::max(delta, std::abs(ustar(k*t, i*h) - A.ar[k][i]));
	}
	std::cout << "\ndelta = " << delta << "\n\n";
}

void Task6()
{
	h = 0.01;
	int n = 1 / h;
	t = h*h / 10;
	T = t*n;
	k = 1;
	std::cout << std::scientific << "Net steps: " << n <<'\n';
	std::cout << std::scientific << "Comparasion: t <= " << ' ' << h*h / 2 << '\n';
	std::cout << std::scientific << "Chosen: t = " << t << '\n';
	std::cout << "Theoretical error: " << h*h + t << '\n';
	ETS();
	//ITS();
}