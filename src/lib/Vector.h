/*
 * Vector.h
 *
 *  Created on: Feb 14, 2011
 *      Author: jason
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <gsl_linalg.h>
#include <gsl_matrix.h>

namespace jason {

class Vector {
public:
	Vector(double *data, size_t size);
	virtual ~Vector();
	void Print();
	friend class Matrix;
private:
	gsl_vector *v;
};

}

#endif /* VECTOR_H_ */
