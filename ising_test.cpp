#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

using namespace std;

int mod(int , int );

int main(){

	unsigned long s=1234;
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	int n,m,y,z,p;
	unsigned long k,c;
	double N,Cc;
	double Esum=0.0,E2sum=0.0,Ssum=0.0,S2sum=0.0;
	double x,E,dE,aE=0,aE2=0,M,C,X,S,dS,aS=0.0,aS2=0.0,aM=0.0,aM2=0.0;

	int q=0;

	ofstream ans ("static_data50.txt");
	ans << "1/beta" << '\t' << "E" <<'\t'<< "M" << '\t'<< "C"<<'\t'<<"X"<<endl;
	n=50; // Number of arrays
	N=static_cast<double>(n);

	vector <vector<double> > vec(n , vector <double>(n)); //Construct a new n x n array;

	E=-2.0*N*N;
	S=N*N;
	//Initialise the array into ordered state
	for(int i=0; i<n; i++){
	for(int j=0; j<n;j++){
		vec[i][j]=1.0;
		}
	}

        double trackE[15000]; //Array used to track Energy
        double trackS[15000]; //Array used to track S

double beta=0.9;

	M=0.0;
	c=0;
	q=0;
	p=0;

	cout<< beta << endl;
	gsl_rng_set(r, k);
	k=0;

	for(int i=0; i<15000; i ++){trackE[i]=0.0; trackS[i]=0.0;} // Initialise random double value in the tracking array

	double lowE=0.0;//First 50 Energy average in trackE
	double midE=0.0;//50-100 Energy average in trackE
	double higE=0.0;// 100-150 energy avrage in trackE

	aE=0.0;
	aE2=0.0;
	aS=0.0;
	aS2=0.0;
	aM=0.0;
	aM2=0.0;

	for(int T=0;T<1000000000000; T++){

		q++;
		k++;
		y=gsl_rng_uniform_int(r,n);
		z=gsl_rng_uniform_int(r,n);

		vec[y][z]=(-1.0)*vec[y][z];
		dE=-2.0*vec[y][z]*(vec[mod((y-1),n)][z]+vec[(y+1)%n][z]+vec[y][mod((z-1),n)]+vec[y][(z+1)%n]);
		dS=2.0*vec[y][z];

		if(dE>0.0){ x=gsl_rng_uniform(r);
			if(x>=exp(-dE*beta)){
						vec[y][z]=(-1.0)*vec[y][z];
						dE=0;
						dS=0;
			}
		}

		E=E+dE;
		S=S+dS;

		trackE[(T%15000)]=E;
		trackS[(T%15000)]=S;

		if(T>100000){
		p++;
		lowE=0.0;
		midE=0.0;
		higE=0.0;

		for (int i=0;i<5000;i++){
					lowE=lowE+trackE[i];
                                        midE=midE+trackE[(i+5000)];
                                        higE=higE+trackE[(i+10000)];}

				if(fabs(lowE-midE)/5000.0<0.05&&fabs(higE-midE)/5000.0<0.05 && c >500){break;}else{c++;}

		Cc= static_cast<double>(p)-1.0 ;

		aE=(aE*Cc+E)/(Cc+1.0);
		aE2=(aE2*Cc+(E*E))/(Cc+1.0);
		aS=(aS*Cc+S)/(Cc+1.0);
		aS2=(aS2*Cc+S*S)/(Cc+1.0);

		M=0.0;

		for (int p=0 ;p<n ; p++) // Calculate Magnetisation
                {
                for (int q=0 ;q<n ; q++)
                {
                M+=vec[p][q];}}

                M=M/(N*N);

		aM=(aM*Cc+M)/(Cc+1.0);
		aM2=(aM2*Cc+M*M)/(Cc+1.0);

		ans<<T<< '\t' << M<< endl;

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

