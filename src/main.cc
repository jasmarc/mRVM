#include <iostream>
#include "Matrix.h"

using namespace std;
using namespace jason;

int main ()
{
	cout << "Starting...\n";
	double m1_arr[] = {1, 2, 9,
	                   4, 6, 5,
	                   3, 7, 4,
			  		   6, 2, 0};
	double m2_arr[] = {1, 2, 9,
	                   4, 6, 5,
	                   3, 7, 4};
	double v_arr[] = {1, 2, 3};
	Matrix *m1 = new Matrix(m1_arr, 4, 3);
	Vector *v1 = m1->Row(0);
	Vector *v2 = m1->Row(2);
	v1->Print();
	v2->Print();
	cout << v1->Multiply(v2);
	cout << "End.\n";
	return 0;
}
