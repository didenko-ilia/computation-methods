#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "Vector.h"

class Matrix
{
public:
	int m, n;
	double ** ar;

	Matrix(int q, int w);

	Matrix(){m = 0; n = 0;};

	~Matrix();

	Matrix(Matrix &other);

	Matrix operator+(Matrix &other);

	Matrix operator*(Matrix &other);

	Matrix operator*(double &k);

	Vector operator*(Vector &vec);

	Matrix operator-(Matrix &other);

	Matrix &operator=(Matrix &other);

	Matrix transpose();

	void print();

	double norm1();

	double norminf(); 

	void init();

	void swapC(int a, int b);

	void printE(Vector q);
};

#endif