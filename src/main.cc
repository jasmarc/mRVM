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
	Matrix *m1 = new Matrix(m_arr, 3, 3);
	Matrix *m2 = m1->Invert();
	m1->Print();
	cout << "\n";
	m2->Print();
	cout << "End.\n";
	return 0;
}
