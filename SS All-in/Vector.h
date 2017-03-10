#ifndef _VECTOR_H_
#define _VECTOR_H_

class Vector
{
public:
	int n;
	double * ar;

	Vector(int q);

	Vector() { n = 0;};

	Vector(Vector &other);

	~Vector();

	Vector operator+(Vector &other);

	double operator*(Vector &other);

	Vector operator*(double k);

	Vector operator-(Vector &other);

	Vector& operator=(Vector &other);

	void print();

	void init();

	double norm1();

	double norm2();

	double norminf();
};

#endif