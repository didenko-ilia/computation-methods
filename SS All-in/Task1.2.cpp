#include "Tasks.h"
#include <iostream>
#include <iomanip>

void gauss(Matrix A)
{
	std::cout << "Gauss:" << std::endl;
	int n = A.m;
	Vector res = Vector(n);
	A.print();
	double tmp, det = 1;
	Vector b = Vector(n);
	for (int i = 0; i < n; i++)
	{
		b.ar[i] = A.ar[i][n];
	}
	for (int q = 0; q < A.m; q++)
	{	
		tmp = A.ar[q][q];
		det *= tmp;

		//if (tmp < eps) std::cout << "Small leading element" << endl;

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
		A.print();
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

	Matrix a = Matrix(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a.ar[i][j] = A.ar[i][j];

	std::cout << "Result:" << std::endl;
	res.print();

	std::cout << "Error:" << std::endl;
	(a*res - b).print();

	//std::cout << "Determinant: " << det << std::endl;
}

void reverse(Matrix A)
{
	std::cout << "Reverse:" << std::endl;
	//std::cout << "GaussC:" << std::endl;
	Matrix Ap = Matrix(A.m, 2 * A.m);
	int n = A.m;
	Vector res = Vector(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Ap.ar[i][j] = A.ar[i][j];
		}
		for (int j = n; j < 2 * n; j++)
		{
			if (i == j-n)
			{
				Ap.ar[i][j] = 1;
			}
			else
			{
				Ap.ar[i][j] = 0;
			}
		}
	}
	Matrix rev = Matrix(n, n);
	Ap.print();
	double tmp;
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
			order[q] = ilead;
			order[ilead] = q;
		}

		tmp = Ap.ar[q][q];

		for (int j = q; j < 2*n; j++)
		{
			Ap.ar[q][j] /= tmp;
		}
		for (int i = q + 1; i < n; i++)
		{
			tmp = Ap.ar[i][q];
			for (int j = q; j < 2*n; j++)
			{
				Ap.ar[i][j] -= Ap.ar[q][j] * tmp;
			}
		}
		Ap.print();
	}

	for (int q = n - 1; q >= 0; q--)
	{
		for (int i = 0; i < n; i++)
		{
			rev.ar[q][i] = Ap.ar[q][n + i];
		}
		for (int w = q + 1; w < n; w++)
		{
			for (int i = 0; i < n; i++)
				rev.ar[q][i] -= Ap.ar[q][w] * rev.ar[w][i];
		}
	}

	rev.print();

	for (int i = 0; i < n; i++)
	{
		if (order[i] != i)
		{
			for (int j = 0; j<n; j++)
			{
				tmp = rev.ar[i][j];
				rev.ar[i][j] = rev.ar[order[i]][j];
				rev.ar[order[i]][j] = tmp;
			}
			order[order[i]] = order[i];
			order[i] = i;
		}
	}

	Matrix a = Matrix(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a.ar[i][j] = A.ar[i][j];

	res = rev*b;
	std::cout << "Result:" << std::endl;
	res.print();

	std::cout << "Error:" << std::endl;
	(a*res - b).print();

	
	std::cout << "Condition number:" << std::endl;
	std::cout << "Using 1 norm:" << std::endl;
	double mu1 = a.norm1()*rev.norm1();
	std::cout << mu1 << std::endl;
	std::cout << "Using inf norm:" << std::endl;
	double muinf = a.norminf()*rev.norminf();
	std::cout << muinf << std::endl;
	std::cout << std::endl;
}

void Task12()
{
	std::cout << "Task 1.2:" << std::setprecision(8) << std::endl;
	Matrix A = Matrix(3, 4);
	A.ar[0][0] = 6.6872;
	A.ar[0][1] = 0.80267;
	A.ar[0][2] = -2.0646;
	A.ar[1][0] = 0.80267;
	A.ar[1][1] = 5.0782;
	A.ar[1][2] = 0.48036;
	A.ar[2][0] = -2.0646;
	A.ar[2][1] = 0.48036;
	A.ar[2][2] = 4.0293;
	A.ar[0][3] = 1.71509;
	A.ar[1][3] = -4.484721;
	A.ar[2][3] = 0.09906;
	//gauss(A);

	//return Vector(gaussC(A));

	reverse(A);
}