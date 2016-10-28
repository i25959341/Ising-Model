// This use the normal ising model as a backbone but we insert an external magnetic field in the array
//Please read the isingN100.ccp for the comments on metropolis aglo and normal ising model comments
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h> //GSL Library
#include <gsl/gsl_math.h>
#include <sstream>

using namespace std;

int mod(int , int );

int main(){

	stringstream ss[30]; //Initialise a stringstream array used for naming the files

	//Initialie the seed and GSL random generator
	unsigned long s=1234;
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	//Defining lattice size and some counters used for debugging
	int n,m,y,z,P;
	unsigned long k,c;
	double N,Cc;
	//Defining Thermodynamics variables
	double Esum=0.0,E2sum=0.0,Ssum=0.0,S2sum=0.0;
	double x,E,dE,aE=0,aE2=0,M,C,X,S,dS,aS=0.0,aS2=0.0,aM=0.0,aM2=0.0;

	int q=0;

	ofstream ans[30];

	n=30; // Number of arrays
	N=static_cast<double>(n);//make the int to double, easier for calculation

	vector <vector<double> > vec(n , vector <double>(n)); //Construct a new n x n array;

        double trackE[15000]; //Array used to track Energy
        double trackS[15000]; //Array used to track S

for (double uB=1.0; uB>=0; uB=uB-0.25 ) {
	//Initialise the energy and spin
        E=-2.0*N*N-N*N*uB;
        S=N*N;
        //Initialise the array into ordered state
        for(int i=0; i<n; i++){
        for(int j=0; j<n;j++){
                vec[i][j]=1.0;
                }
        }

	B=static_cast<int>((uB+1.0)*10.0); //Cast an interger value for different magnetics for file name

	ss[B] <<"mag"<<uB<<".txt" ; //Construct a file name using string

        ans[B].open(ss[B].str().c_str()); // outstream different file for different B

        ans[B] << "1/beta" << '\t' << "E" <<'\t'<< "M" << '\t'<< "C"<<'\t'<<"X"<<endl;

for (double beta=2.0; beta>0.1; beta-=0.005){
	//clearn the variables
	M=0.0;
	c=0;
	q=0;
	P=0;

	cout<< beta << endl;
	gsl_rng_set(r, k);
	k=0;

	for(int i=0; i<15000; i ++){trackE[i]=0.0; trackS[i]=0.0;} // Initialise random double value in the tracking array

	double lowE=0.0;//First 50 Energy average in trackE
	double midE=0.0;//50-100 Energy average in trackE
	double higE=0.0;// 100-150 energy avrage in trackE
	//clean the variables
	aE=0.0;
	aE2=0.0;
	aS=0.0;
	aS2=0.0;
	aM=0.0;
	aM2=0.0;

	for(int T=0;T<1000000000000000000; T++){	//Monte Carlo method starts

		q++;
		k++;
		y=gsl_rng_uniform_int(r,n);	//randomly pick a random site to flip
		z=gsl_rng_uniform_int(r,n);

		vec[y][z]=(-1.0)*vec[y][z];	//Flip the spin and calculate the change
		dS=2.0*vec[y][z];
		dE=-2.0*vec[y][z]*(vec[mod((y-1),n)][z]+vec[(y+1)%n][z]+vec[y][mod((z-1),n)]+vec[y][(z+1)%n])-uB*dS;
		//Metropolis Algorithm starts
		if(dE>0.0){ x=gsl_rng_uniform(r);
			if(x>=exp(-dE*beta)){
						vec[y][z]=(-1.0)*vec[y][z];
						dE=0;
						dS=0;
			}
		}

		for (int p=0 ;p<n ; p++) // Calculate Magnetisation
                {
                for (int q=0 ;q<n ; q++)
                {
                M+=vec[p][q];}}
		//update the spin and energy
		S=S+dS;
		E=E+dE;
		//Track the energy and spin using array
		trackE[(T%15000)]=E;
		trackS[(T%15000)]=S;

		if(T>100000){//Equilibrium condition starts after certain amount of Monte carlo step
		P++;
		//clean the variables
		lowE=0.0;
		midE=0.0;
		higE=0.0;
		//Equilibrium condition algo begins
		for (int i=0;i<5000;i++){
					//split the Energy array into three sections
					lowE=lowE+trackE[i];
                                        midE=midE+trackE[(i+5000)];
                                        higE=higE+trackE[(i+10000)];}
				//if the differetn of energy at each sections is smaller than a certain amount and happens a multiple times, BREAK
				if(fabs(lowE-midE)/5000.0<0.05&&fabs(higE-midE)/5000.0<0.05 && c >500){break;}else{c++;}

		Cc= static_cast<double>(P)-1.0 ;//keep tracking of the number of steps used to calculate average
		//Average calculating
		aE=(aE*Cc+E)/(Cc+1.0);
		aE2=(aE2*Cc+(E*E))/(Cc+1.0);
		aS=(aS*Cc+S)/(Cc+1.0);
		aS2=(aS2*Cc+S*S)/(Cc+1.0);

		aM=(aM*Cc+M)/(Cc+1.0);
		aM2=(aM2*Cc+M*M)/(Cc+1.0);

	}

	}
		cout<<q<<endl;
		//calculate the therodynamics variables
		C=(aE2-aE*aE)*beta*beta/(N*N);
		X=(aS2-aS*aS)*beta/(N*N);

		ans[B] <<1.0/beta << '\t' << aE <<'\t'<< aM <<'\t'<< aM2 <<'\t'<< C<<'\t'<<X<<endl;

		for(int y=0;y<n;y++)
			{
			for(int x=0;x<n;x++){
			if(vec[x][y]<0){
				cout<<" − ";}
			else{
				cout<<" + ";}
			}
			cout<<endl;
		}
}
}
return 0;

}

int mod(int a, int b){ //Mod Calculator

	int m;
        if (a<0)
        { m = b+ (a%b) ; }
        else
        {m =a%b; }
        return m;
}

