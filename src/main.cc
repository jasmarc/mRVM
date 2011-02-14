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
	Matrix *m = new Matrix(m_arr, 3, 3);
	m->Print();
	cout << "\n";
	Vector *v = m->Row(0);
	v->Print();
	cout << "\n";
	v = m->Column(0);
	v->Print();
	cout << "End.\n";
	return 0;
}
