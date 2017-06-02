#include "Tasks.h"
#include <iostream>
#include <iomanip>

const int N = 11;
const double PI = 3.14159;

Vector SLE(Vector a, Vector b, Vector c, Vector d)
{
	Matrix A(N, N);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A.ar[i][j] = 0;
		}
		if (a.ar[i] != 0)
			A.ar[i][i - 1] = a.ar[i];
		A.ar[i][i] = b.ar[i];
		if (c.ar[i] != 0)
			A.ar[i][i + 1] = c.ar[i];
	}
	A.printE(d);
	double m[N], k[N];
	m[0] = -c.ar[0] / b.ar[0];
	k[0] = d.ar[0] / b.ar[0];
	std::cout << "\n0 " << std::scientific << m[0] << ' ' << k[0] << '\n';
	for (int i = 1; i <= N - 1; ++i)
	{
		m[i] = -c.ar[i - 1] / (a.ar[i - 1] * m[i - 1] + b.ar[i - 1]);
		k[i] = (d.ar[i-1]-k[i-1]*a.ar[i-1])/ (a.ar[i - 1] * m[i - 1] + b.ar[i - 1]);
		std::cout << i << std::scientific << ' ' << m[i] << ' ' << k[i] << '\n';
	}
	Vector y(N);
	y.ar[N - 1] = (d.ar[N - 1] - a.ar[N - 1] * k[N - 1]) / (a.ar[N - 1] * m[N - 1] + b.ar[N - 1]);
	for (int i = N - 2; i >= 0; --i)
	{
		y.ar[i] = m[i + 1] * y.ar[i + 1] + k[i + 1];
	}
	std::cout << "Result:\n";
	y.printH();
	std::cout << '\n';
	std::cout << "Error:\n";
	std::cout << std::scientific <<  (A*y - d).norminf();
	std::cout << '\n';
	return Vector(y);
}

void BVP(double o, double w, double alpha, double beta)
{
	Vector a(N), b(N), c(N), d(N);
	double p[N], q[N], r[N], f[N];
	for (int i = 0; i < N-1; i++)
	{
		double xi = o + i*(w - o) / (N - 1);
		p[i] = 1;
		q[i] = log(xi+2);
		r[i] = -xi;
		f[i] = xi + 2;
	}
	double h = (w - 0) / (N - 1);
	a.ar[0] = 0;
	b.ar[0] = alpha*h-1;
	c.ar[0] = 1;
	d.ar[0] = 0;
	a.ar[N - 1] = -1;
	b.ar[N - 1] = h*beta + 1;
	c.ar[N - 1] = 0;
	d.ar[N - 1] = 0;
	for (int i = 1; i < N - 1; ++i)
	{
		a.ar[i] = p[i] - h*q[i] / 2;
		b.ar[i] = -2 * p[i] + h*h*r[i];
		c.ar[i] = p[i] + h*q[i] / 2;
		d.ar[i] = h*h*f[i];
	}
	Vector res[2];
	res[0] = SLE(a, b, c, d);

	b.ar[0] = 2 *alpha* h + (a.ar[1] / c.ar[1] - 3);
	c.ar[0] = (b.ar[1] / c.ar[1] + 4);
	d.ar[0] = d.ar[1] / c.ar[1];
	a.ar[N - 1] = -(4 + b.ar[N - 2] / a.ar[N - 2]);
	b.ar[N - 1] = 2 * beta*h + (3 - c.ar[N - 2] / a.ar[N - 2]);
	d.ar[N - 1] = -d.ar[N - 2] / a.ar[N - 2];

	res[1] = SLE(a, b, c, d);
	
	std::cout << "\nResults:\n";
	res[0].printH();
	res[1].printH();
	std::cout << "\nError: ";
	std::cout << (res[0] - res[1]).norminf() << '\n';
}

void Task4()
{
	BVP(0, 1, -0.6, 0.4);
}
