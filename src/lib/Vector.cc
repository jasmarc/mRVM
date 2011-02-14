/*
 * Vector.cpp
 *
 *  Created on: Feb 14, 2011
 *      Author: jason
 */

#include "Vector.h"

namespace jason {

Vector::Vector(double *data, size_t size) {
	this->v = gsl_vector_alloc(size);
	for (size_t i = 0; i < size; ++i) {
		gsl_vector_set(this->v, i, data[i]);
	}
}

void Vector::Print() {
	for (size_t j = 0; j < v->size; j++) {
		printf("%.2f ", gsl_vector_get(v, j));
	}
	printf("\n");
}

Vector::~Vector() {
	gsl_vector_free(this->v);
}

}
