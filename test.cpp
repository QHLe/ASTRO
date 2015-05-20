#include <iostream>
#include "SVector.hpp"
#include "MVector.hpp"

using namespace std;


int main(){

	SVector vector;
	MVector mvector;
	vector = linspace(1,0.1,10);
	vector.Print();	
//	vector(2) = 3;
//	vector.Print();
	cout<<endl;
	mvector = vector;
	mvector.Print();
	return 0;
}
