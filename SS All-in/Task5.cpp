#include "Tasks.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

Vector gaussC(Matrix A)
{
	//std::cout << "GaussC:" << std::endl;
	int n = A.m;
	Vector res = Vector(n);
	//A.print();
	double tmp, det = 1;
	Vector b = Vector(n);

	for (int i = 0; i < n; i++)
	{
		b.ar[i] = A.ar[i][n];
	}

	int * order;
	order = new int[n];
	for (int i = 0; i < n; i++)
	{
		order[i] = i;
	}

	for (int q = 0; q < n; q++)
	{
		double lead = 0;
		int ilead;
		for (int j = q; j < n; j++)
		{
			double a = abs(A.ar[q][j]);
			if (a > lead)
			{
				lead = a;
				ilead = j;
			}
		}
		if (ilead != q)
		{
			A.swapC(q, ilead);
			det *= -1;
			order[q] = ilead;
			order[ilead] = q;
		}

		tmp = A.ar[q][q];
		det *= tmp;

		for (int j = q; j < n + 1; j++)
		{
			A.ar[q][j] /= tmp;
		}
		for (int i = q + 1; i < n; i++)
		{
			tmp = A.ar[i][q];
			for (int j = q; j < n + 1; j++)
			{
				A.ar[i][j] -= A.ar[q][j] * tmp;
			}
		}
		//A.print();
	}

	for (int q = n - 1; q >= 0; q--)
	{
		tmp = A.ar[q][n];
		for (int w = q + 1; w < n; w++)
		{
			tmp -= A.ar[q][w] * res.ar[w];
		}
		res.ar[q] = tmp;
	}

	for (int i = 0; i < n; i++)
	{
		if (order[i] != i)
		{
			tmp = res.ar[i];
			res.ar[i] = res.ar[order[i]];
			res.ar[order[i]] = tmp;
			order[order[i]] = order[i];
			order[i] = i;
		}
	}

	return Vector(res);
}

Vector gauss4t()
{
	Vector res(4);
	res.ar[0] = -std::sqrt(3.0 / 7 + 2.0 / 7 * std::sqrt(6.0 / 5));
	res.ar[1] = -std::sqrt(3.0 / 7 - 2.0 / 7 * std::sqrt(6.0 / 5));
	res.ar[2] = std::sqrt(3.0 / 7 - 2.0 / 7 * std::sqrt(6.0 / 5));
	res.ar[3] = std::sqrt(3.0 / 7 + 2.0 / 7 * std::sqrt(6.0 / 5));
	return res;
}
Vector gauss4A()
{
	Vector res(4);
	res.ar[0] = (18 - std::sqrt(30.0)) / 36;
	res.ar[1] = (18 + std::sqrt(30.0)) / 36;
	res.ar[2] = (18 + std::sqrt(30.0)) / 36;
	res.ar[3] = (18 - std::sqrt(30.0)) / 36;
	return res;
}
Vector gauss5t()
{
	Vector res(5);
	res.ar[0] = -std::sqrt(5 + 2 * std::sqrt(10.0 / 7))/3;
	res.ar[1] = -std::sqrt(5 - 2 * std::sqrt(10.0 / 7))/3;
	res.ar[2] = 0;
	res.ar[3] = std::sqrt(5 - 2 * std::sqrt(10.0 / 7))/3;
	res.ar[4] = std::sqrt(5 + 2 * std::sqrt(10.0 / 7))/3;
	return res;
}
Vector gauss5A()
{
	Vector res(5);
	res.ar[0] = (322 - 13 * sqrt(70.0)) / 900;
	res.ar[1] = (322 + 13 * sqrt(70.0)) / 900;
	res.ar[2] = 128.0 / 225;
	res.ar[3] = (322 + 13 * sqrt(70.0)) / 900;
	res.ar[4] = (322 - 13 * sqrt(70.0)) / 900;
	return res;
}

double K(double x, double y)
{
	return 1 / (2 + x*x + y*y);
}

double f(double x)
{
	return  1 - x + x*x;
}

Vector ufour(4), ufive(5);

double u4(double x)
{
	double res = 0;
	res += f(x);
	for (int i = 0; i < 4; ++i)
	{
		double xi = gauss4t().ar[i];
		res -= 0.5*gauss4A().ar[i] * K(x, xi) * ufour.ar[i];
	}
	return res;
}

double u5(double x)
{
	double res = 0;
	res += f(x);
	for (int i = 0; i < 5; ++i)
	{
		double xi = gauss5t().ar[i];
		res -= 0.5*gauss5A().ar[i] * K(x, xi) * ufive.ar[i];
	}
	return res;
}
void MMK()
{
	std::cout << "MMK\n";
	Matrix A(4, 5);
	Vector x = gauss4t();
	for (int i = 0; i < 4; ++i)
		x.ar[i] = x.ar[i] * 0.5 + 0.5;
	for (int i = 0; i<4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			A.ar[i][j] = gauss4A().ar[j] * 0.5*K(x.ar[i], x.ar[j]) + (i==j);
		}
	for (int i = 0; i < 4; ++i)
		A.ar[i][4] = f(x.ar[i]);
	ufour = gaussC(A);
	ufour.printH();

	Matrix B(5, 6);
	x = gauss5t();
	for (int i = 0; i < 5; ++i)
		x.ar[i] = x.ar[i] * 0.5 + 0.5;
	for (int i = 0; i<5; ++i)
		for (int j = 0; j < 5; ++j)
		{
			B.ar[i][j] = gauss5A().ar[j] * 0.5*K(x.ar[i], x.ar[j]) + (i == j);
		}
	for (int i = 0; i < 5; ++i)
		B.ar[i][5] = f(x.ar[i]);
	ufive = gaussC(B);

	for (int i = 0; i <= 10; ++i)
	{
		std::cout << std::scientific << i*0.1 << ' ' << u4(i*0.1) << ' ' << u5(i*0.1) << '\n';
	}
};

Vector A3, A4;

double a0(double x)
{
	return 1;
}

double a1(double x)
{
	return -x*x;
}

double a2(double x)
{
	return x*x*x*x;
}

double a3(double x)
{
	return -x*x*x*x*x*x;
}

double U3(double x)
{
	double res = f(x);
	res -= a0(x)*A3.ar[0];
	res -= a1(x)*A3.ar[1];
	res -= a2(x)*A3.ar[2];
	return res;
}

double U4(double x)
{
	double res = f(x);
	res -= a0(x)*A4.ar[0];
	res -= a1(x)*A4.ar[1];
	res -= a2(x)*A4.ar[2];
	res -= a3(x)*A4.ar[3];
	return res;
}

double b0(double x)
{
	return 1 / (x*x + 2);
}

double b1(double x)
{
	return 1 / ((x*x + 2)*(x*x + 2));
}

double b2(double x)
{
	return 1 / ((x*x + 2)*(x*x + 2)*(x*x + 2));
}

double b3(double x)
{
	return 1 / ((x*x + 2)*(x*x + 2)*(x*x + 2)*(x*x + 2));
}

double iSimpson(double(*f1)(double x), double(*f2)(double x), double a, double b)
{
	return (b - a) / 6 * (f1(a)*f2(a) + 4 * f1((b - a) / 2)*f2((b - a) / 2) + f1(b)*f2(b));
}

void KS()
{
	std::cout << "\nKernel Substitution\n";
	Matrix A(3, 4);
	A.ar[0][3] = iSimpson(&f, &b0, 0, 1);
	A.ar[1][3] = iSimpson(&f, &b1, 0, 1);
	A.ar[2][3] = iSimpson(&f, &b2, 0, 1);
	A.ar[0][0] = 1 + iSimpson(&a0, &b0, 0, 1);
	A.ar[1][0] = iSimpson(&a0, &b1, 0, 1);
	A.ar[2][0] = iSimpson(&a0, &b2, 0, 1);
	A.ar[0][1] = iSimpson(&a1, &b0, 0, 1);
	A.ar[1][1] = 1 + iSimpson(&a1, &b1, 0, 1);
	A.ar[2][1] = iSimpson(&a1, &b2, 0, 1);
	A.ar[0][2] = iSimpson(&a2, &b0, 0, 1);
	A.ar[1][2] = iSimpson(&a2, &b1, 0, 1);
	A.ar[2][2] = 1 + iSimpson(&a2, &b2, 0, 1);
	A3 = gaussC(A);

	Matrix B(4, 5);
	B.ar[0][4] = iSimpson(&f, &b0, 0, 1);
	B.ar[1][4] = iSimpson(&f, &b1, 0, 1);
	B.ar[2][4] = iSimpson(&f, &b2, 0, 1);
	B.ar[3][4] = iSimpson(&f, &b3, 0, 1);
	B.ar[0][0] = 1 + iSimpson(&a0, &b0, 0, 1);
	B.ar[1][0] = iSimpson(&a0, &b1, 0, 1);
	B.ar[2][0] = iSimpson(&a0, &b2, 0, 1);
	B.ar[3][0] = iSimpson(&a0, &b3, 0, 1);
	B.ar[0][1] = iSimpson(&a1, &b0, 0, 1);
	B.ar[1][1] = 1 + iSimpson(&a1, &b1, 0, 1);
	B.ar[2][1] = iSimpson(&a1, &b2, 0, 1);
	B.ar[3][1] = iSimpson(&a1, &b3, 0, 1);
	B.ar[0][2] = iSimpson(&a2, &b0, 0, 1);
	B.ar[1][2] = iSimpson(&a2, &b1, 0, 1);
	B.ar[2][2] = 1 + iSimpson(&a2, &b2, 0, 1);
	B.ar[3][2] = iSimpson(&a2, &b3, 0, 1);
	B.ar[0][3] = iSimpson(&a3, &b0, 0, 1);
	B.ar[1][3] = iSimpson(&a3, &b1, 0, 1);
	B.ar[2][3] = iSimpson(&a3, &b2, 0, 1);
	B.ar[3][3] = 1 + iSimpson(&a3, &b3, 0, 1);
	A4 = gaussC(B);

	double delta = 0;
	for (int i = 0; i <= 10; ++i)
	{
		double res1 = U3(i*0.1), res2 = U4(i*0.1);
		delta = std::max(delta, std::abs(res1 - res2));
		std::cout << std::scientific << i*0.1 << ' ' << res1 << ' ' << res2 << '\n';
	}
	std::cout << "delta = " << std::scientific << delta << '\n';
};

void Task5()
{
	MMK();
	KS();
}