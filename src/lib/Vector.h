#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <gsl_matrix.h>

namespace jason {

class Vector {
public:
	Vector(double *data, size_t size);
	virtual ~Vector();
	void Print();
	friend class Matrix;
private:
	Vector(gsl_vector *v);
	gsl_vector *v;
};

}

#endif /* VECTOR_H_ */
