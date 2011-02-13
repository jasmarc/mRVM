#include "Matrix.h"
#include <iostream>

using namespace std;

namespace jason {

Matrix::Matrix(const char* filename) {
	cout << "read " << filename << "\n";
	int rows, cols;
	FILE *f;
	f = fopen(filename, "r");
	rows = NumberOfRows(f);
	cols = NumberOfColumns(f);
	this->m = gsl_matrix_alloc(rows, cols);
	gsl_matrix_fscanf(f, this->m);
	fclose(f);
}

Matrix::~Matrix() {
	gsl_matrix_free(this->m);
}

void Matrix::Print() {
	for (size_t i = 0; i < this->m->size1; ++i) {
		for (size_t j = 0; j < this->m->size2; ++j) {
			printf("%.2f ", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
}

int Matrix::NumberOfRows(FILE *f) {
	char lastChar = '\n';
	char currentChar = NULL;
	int count = 0;
	while ((currentChar = fgetc(f)) != EOF) {
		if (lastChar == '\n' && currentChar != '\n') {
			++count;
		}
		lastChar = currentChar;
	}
	fseek(f, 0, SEEK_SET);
	return count;
}

int Matrix::NumberOfColumns(FILE *f) {
	char lastChar = ' ';
	char currentChar = NULL;
	int count = 0;
	while ((currentChar = fgetc(f)) != '\n') {
		if (isspace(lastChar) && !isspace(currentChar)) {
			++count;
		}
		lastChar = currentChar;
	}
	fseek(f, 0, SEEK_SET);
	return count;
}

}
