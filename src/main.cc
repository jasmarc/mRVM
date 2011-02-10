#include <iostream>
#include "Matrix.h"

using namespace std;
using namespace jason;

int main ()
{
	cout << "Starting...\n";
	Matrix *m = new Matrix("./src/test/test.dat");
	m->Print();
	cout << "End.\n";
	return 0;
}
