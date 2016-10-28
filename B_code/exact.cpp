// This is use to calculate the exact result for the 2x2 Ising model and compare with the Monte Carlo method
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
	//Introduce the Stat Mec Variables
	double Z,aE,aE2,aM2,aM;
	double C,X;

	int i;
	double x=0, y=0;

	ofstream ans ("exact.txt");
	//For loop over different beta
	for(double beta=0.1; beta<=1.0;beta=beta+0.001){
	//Caculate the exact result for Magnetisation and Susceptability
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
