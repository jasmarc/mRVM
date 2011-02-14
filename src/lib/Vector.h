#ifndef VECTOR_H_
#define VECTOR_H_

#include <gsl_matrix.h>
#include "Matrix.h"

namespace jason {

class Matrix;

class Vector {
public:
	Vector(double *data, size_t size);
	virtual ~Vector();
	Matrix *RepmatVert(size_t k);
	Matrix *RepmatHoriz(size_t k);
	double Multiply(Vector *other);
	void Print();
	friend class Matrix;
private:
	Vector(gsl_vector *v);
	gsl_vector *v;
};

}

#endif /* VECTOR_H_ */
