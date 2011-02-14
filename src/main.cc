#include <iostream>
#include "Matrix.h"

using namespace std;
using namespace jason;

int main ()
{
	cout << "Starting...\n";
	double m_arr[] = {1, 2, 9,
	                  4, 6, 5,
	                  3, 7, 4};
	double v_arr[] = {1, 2, 3};
	Vector *v = new Vector(v_arr, 3);
	v->Print();
	cout << "\n";
	Matrix *m = new Matrix(v);
	m->Print();
	cout << "End.\n";
	return 0;
}
