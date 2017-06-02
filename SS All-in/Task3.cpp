#include "Tasks.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
Matrix E = Matrix(3, 3);

const double eps1 = 1e-7;
const double eps2 = 1e-5;

double powersMethod(Matrix A, Vector y, int q)
{
	double l, l1;
	Vector y1;
	y1 = A*y;
	l1 = y1.ar[q] / y.ar[q];
	y1 = y1 * (1 / y1.norminf());
	cout << setprecision(9) << scientific;
	cout << "Powers Method" << endl;
	cout << "First approximation: ";
	y.printH();
	cout << endl;
	cout << scientific << "Precision: " << eps1 << endl;
	cout << 1 << ' '  << l1 << ' ';
	y1.printH();
	cout << endl;
	int N = 1;
	do
	{
		N++;
		y = y1;
		l = l1;
		y1 = A*y;
		l1 = y1.ar[q] / y.ar[q];
		y1 = y1 * (1 / y1.norminf());
		cout << N << ' ' << scientific << l1 << ' ';
		y1.printH();
		cout << endl;
	} while (abs(l1 - l) > eps1);
	cout << "Error: " << ((A*y1) - (y1*l1)).norminf() << endl;
	cout << endl;
	return l1;
}

double scalarMethod(Matrix A, Vector y)
{
	double l, l1;
	Vector y1;
	y1 = A*y;
	l1 = (y1*y)/(y*y);
	y1 = y1 * (1 / y1.norminf());
	cout << setprecision(9);
	cout << "Scalar Method" << endl;
	cout << "First approximation: ";
	y.printH();
	cout << endl;
	cout << scientific << "Precision: " << eps1 << endl;
	cout << 1 << ' ' << scientific << l1 << ' ';
	y1.printH();
	cout << endl;
	int N = 1;
	do
	{
		N++;
		y = y1;
		l = l1;
		y1 = A*y;
		l1 = (y1*y) / (y*y);
		y1 = y1 * (1 / y1.norminf());
		cout << N << ' ' << scientific << l1 << ' ';
		y1.printH();
		cout << endl;
	} while (abs(l1 - l) > eps1);
	cout << "Error: " << ((A*y1) - (y1*l1)).norminf() << endl;
	cout << endl;
	return l1;
}

double oppositeEndOfSpectre(Matrix A, double a)
{
	cout << "Opposite Side of Spectre" << endl;
	Matrix B = A - E*a;
	Vector u = Vector(3);
	u.ar[0] = 1;
	u.ar[1] = 1;
	u.ar[2] = 1;
	double b = powersMethod(B, u, 2) + a;
	cout << "Opposite Side: " << b << endl;
	return b;
}

int sgn(double a)
{
	return ((a > 0) - (a < 0));
}

double lastEigenvalue(Matrix A, double a, double b)
{
	double sum = 0;
	for (int i = 0; i < A.n; i++)
		sum += A.ar[i][i];
	sum -= a + b;
	cout << "Last Eigenvalue: " << sum << endl;
	return sum;
}

void jacobi(Matrix A)
{
	Matrix X = E;
	Matrix S = A;
	A.print();
	int maxi, maxj;
	double max;
	int k = 0;
	do
	{
		k++;
		max = 0;
		for (int i = 0; i < A.n; i++)
			for (int j = i + 1; j < A.n; j++)
				if (abs(A.ar[i][j]) > max)
				{
					max = abs(A.ar[i][j]);
					maxi = i;
					maxj = j;
				}
		if (max < eps2)
			break;
		double d = sqrt((A.ar[maxi][maxi] - A.ar[maxj][maxj])*(A.ar[maxi][maxi] - A.ar[maxj][maxj]) + 4 * A.ar[maxi][maxj] * A.ar[maxi][maxj]);
		double c = sqrt(0.5*(1 + abs(A.ar[maxi][maxi] - A.ar[maxj][maxj]) / d));
		double s = sgn(A.ar[maxi][maxj] * (A.ar[maxi][maxi] - A.ar[maxj][maxj]))*sqrt(0.5*(1 - abs(A.ar[maxi][maxi] - A.ar[maxj][maxj]) / d));
		Matrix V = Matrix(3, 3);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
					V.ar[i][j] = 1;
				else
					V.ar[i][j] = 0;
			}
		}
		V.ar[maxi][maxi] = c;
		V.ar[maxj][maxj] = c;
		V.ar[maxi][maxj] = -s;
		V.ar[maxj][maxi] = s;
		X = X*V;
		Matrix A1 = A;
		//A1.print();
		for (int i = 0; i < A.n; i++)
		{
			if (i != maxi)
			{
				A.ar[i][maxi] = c*A1.ar[i][maxi] + s*A1.ar[i][maxj];
				A.ar[maxi][i] = c*A1.ar[i][maxi] + s*A1.ar[i][maxj];
			}
			if (i != maxj)
			{
				A.ar[i][maxj] = -s*A1.ar[i][maxi] + c*A1.ar[i][maxj];
				A.ar[maxj][i] = -s*A1.ar[i][maxi] + c*A1.ar[i][maxj];
			}
		}
		A.ar[maxi][maxi] = c*c*A1.ar[maxi][maxi] + 2 * c*s *A1.ar[maxi][maxj] + s*s*A1.ar[maxj][maxj];
		A.ar[maxj][maxj] = s*s*A1.ar[maxi][maxi] - 2 * c*s*A1.ar[maxi][maxj] + c*c*A1.ar[maxj][maxj];
		A.ar[maxi][maxj] = 0;
		A.ar[maxj][maxi] = 0;
		cout << "Step " << k << endl;
		A.print();


	} while (max > eps2);

	Vector x1 = Vector(3), x2 = Vector(3), x3 = Vector(3);
	for (int i = 0; i < 3; i++)
	{
		x1.ar[i] = X.ar[i][0];
		x2.ar[i] = X.ar[i][1];
		x3.ar[i] = X.ar[i][2];
	}
	cout << "Eigenvalue: " << A.ar[0][0] << "\tEigenvector: ";
	x1.printH();
	cout << endl;
	cout << "Error: " << ((S*x1) - (x1*A.ar[0][0])).norminf() << endl;
	cout << endl;
	cout << "Eigenvalue: " << A.ar[1][1] << "\t Eigenvector: ";
	x2.printH();
	cout << endl;
	cout << "Error: " << ((S*x2) - (x2*A.ar[1][1])).norminf() << endl;
	cout << endl;
	cout << "Eigenvalue: " << A.ar[2][2] << "\tEigenvector: ";
	x3.printH();
	cout << endl;
	cout << "Error: " << ((S*x3) - (x3*A.ar[2][2])).norminf() << endl;
}

void Task3()
{
	Matrix A = Matrix(3,3);
	A.ar[0][0] = -0.820053;
	A.ar[0][1] = -0.135416;
	A.ar[0][2] = 0.269479;
	A.ar[1][0] = -0.135416;
	A.ar[1][1] = 0.514864;
	A.ar[1][2] = 0.027061;
	A.ar[2][0] = 0.269479;
	A.ar[2][1] = 0.027061;
	A.ar[2][2] = -0.833651;

	Vector y = Vector(3);
	y.ar[0] = 1;
	y.ar[1] = 1;
	y.ar[2] = 1;

	E.ar[0][0] = 1;
	E.ar[0][1] = 0;
	E.ar[0][2] = 0;
	E.ar[1][0] = 0;
	E.ar[1][1] = 1;
	E.ar[1][2] = 0;
	E.ar[2][0] = 0;
	E.ar[2][1] = 0;
	E.ar[2][2] = 1;

	double a;
	a = powersMethod(A, y, 2);

	scalarMethod(A, y);

	double b;
	b = oppositeEndOfSpectre(A, a);

	double c;
	c = lastEigenvalue(A, a, b);

	jacobi(A);
}