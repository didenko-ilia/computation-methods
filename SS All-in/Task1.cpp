#include "Tasks.h"
#include <iostream>

void Task1()
{
	Matrix A = Matrix(2,2), Ainv = Matrix(2,2);
	A.ar[0][0] = 1.00;
	A.ar[0][1] = 0.99;
	A.ar[1][0] = 0.99;
	A.ar[1][1] = 0.98;
	Vector b = Vector(2), db = Vector(2);
	b.ar[0] = 1.99;
	b.ar[1] = 1.97;
	db.ar[0] = -0.000097;
	db.ar[1] = 0.000106;
	Vector x = Vector(2), xd = Vector(2);
	x.ar[0] = 1;
	x.ar[1] = 1;
	xd.ar[0] = 3.0000;
	xd.ar[1] = -1.0203;

	std::cout << "x:" << std::endl;
	x.print();
	std::cout << "x':" << std::endl;
	xd.print();
	std::cout << "r:" << std::endl;
	Vector r = Vector(2);
	r = (A*x) -b;
	r.print();
	std::cout << "r':" << std::endl;
	Vector rd = Vector(2);
	rd = (A*xd) - b - db;
	rd.print();

	std::cout << "Real relative error:" << std::endl;
	std::cout << "Using 1 norm:" << std::endl;
	std::cout << (xd - x).norm1() / x.norm1() << std::endl;
	std::cout << "Using 2 norm:" << std::endl;
	std::cout << (xd - x).norm2() / x.norm2() << std::endl;
	std::cout << "Using inf norm:" << std::endl;
	std::cout << (xd - x).norminf() / x.norminf() << std::endl;
	std::cout << std::endl;

	double dA = A.ar[0][0] * A.ar[1][1] - A.ar[0][1] * A.ar[1][0];
	dA = 1 / dA;
	Ainv.ar[0][0] = A.ar[1][1];
	Ainv.ar[1][1] = A.ar[0][0];
	Ainv.ar[1][0] = -A.ar[0][1];
	Ainv.ar[0][1] = -A.ar[1][0];
	Ainv = Ainv * dA;

	std::cout << "Condition number:" << std::endl;
	std::cout << "Using 1 norm:" << std::endl;
	double mu1 = A.norm1()*Ainv.norm1();
	std::cout << mu1 << std::endl;
	std::cout << "Using inf norm:" << std::endl;
	double muinf = A.norminf()*Ainv.norminf();
	std::cout << muinf << std::endl;
	std::cout << std::endl;

	std::cout << "Theoretical relative error:" << std::endl;
	std::cout << "Using 1 norm:" << std::endl;
	std::cout << mu1*db.norm1() / b.norm1() << std::endl;
	std::cout << "Using inf norm:" << std::endl;
	std::cout << muinf*db.norminf() / b.norminf() << std::endl;

	int q;
	std::cin >> q;
}