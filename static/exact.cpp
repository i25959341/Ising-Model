#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int main()
{

	double Z,aE,aE2,aM2,aM;
	double C,X;
	int i;
	double x=0, y=0;

	ofstream ans ("exact.txt");
  //a simple loop that sums the squares of the counter
  for(double beta=0.1; beta<=1.0;beta=beta+0.001){

	Z=2.0*exp(8.0*beta)+12.0+2.0*exp(-8.0*beta);
	aE=-(1/Z)*(16.0*exp(8.0*beta)-16.0*exp(-8.0*beta));
	aE2= (1/Z)*(128.0*exp(8.0*beta)+128.0*exp(-8.0*beta));
	aM= (1/Z)*(8.0*exp(8.0*beta)+16.0);
	aM2= (1/Z)*(32.0*exp(8.0*beta)+32.0);

	C=(1.0/4.0)*beta*beta*(aE2-aE*aE);
	X=(1.0/4.0)*beta*beta*(aM2-aM*aM);

    ans << 1.0/beta << "\t" << aE << "\t" << C<<  "\t" << X << endl;

  }

  //exit
  return 0;
}
