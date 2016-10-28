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


	int n,m,y,z;
	unsigned long k;
	double N;
	double Esum=0.0,E2sum=0.0,Ssum=0.0,S2sum=0.0;
	double x,E,dE,aE=0,aE2=0,M,C,X,S,dS,aS=0.0,aS2=0.0;

	ofstream ans ("test_data.txt");
	ans << "1/beta" << '\t' << "E" <<'\t'<< "M" << '\t'<< "C"<<'\t'<<"X"<<endl;
	n=2; // Number of arrays
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

double beta=1.0;

	M =0.0;
	k=0;
	cout<< beta << endl;
	gsl_rng_set(r, k);

	double trackE[1500]; //Array used to track Energy
	double trackS[1500]; //Array used to track S

	for(int i=0; i<1500; i ++){trackE[i]=0.0; trackS[i]=0.0;} // Initialise random double value in the tracking array

	double lowE=0.0;//First 50 Energy average in trackE
	double midE=0.0;//50-100 Energy average in trackE
	double higE=0.0;// 100-150 energy avrage in trackE

	for(int T=0;T <1000000000; T++){

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

		k++;

		E=E+dE;
		S=S+dS;

		if(T<=15000){
		trackE[(T%1500)]=E;
		trackS[(T%1500)]=S;}

		if(T>15000){

		for (int i=0; i<1500;i++){trackE[i%1500]=trackE[(i+1)%1500];
					trackS[i%1500]=trackS[(i+1)%1500];}

		trackE[1499]=E;
		trackS[1499]=S;

		for (int i=0;i<500;i++){

					lowE=lowE+trackE[i%1500];
                                        midE=midE+trackE[(i+500)%1500];
                                        higE=higE+trackE[(i+1000)%1500];}

				if(fabs(lowE-midE)/500.0<0.01||fabs(higE-midE)/500.0<0.01){break;}
			}

		for (int p=0 ;p<n ; p++) // Calculate Magnetisation
                {
                for (int q=0 ;q<n ; q++)
                {
                M+=vec[p][q];}}

                M=M/(N*N);

	ans << T << '\t' << M << endl;
	}
	        cout << k << endl;
		                for(int a=0;a<1500;a++){Esum=Esum+trackE[a];
                                        E2sum=E2sum+trackE[a]*trackE[a];
                                        Ssum=Ssum+trackS[a];
                                        S2sum=S2sum+trackS[a]*trackS[a];
                                        }

		aE=Esum/1500.0;
		aE2=E2sum/1500.0;
		aS=Ssum/1500.0;
		aS2=S2sum/1500.0;

		C=(aE2-aE*aE)/((N*N)*(beta*beta));
		X=(aS2-aS*aS)*beta/(N*N);

		for (int p=0 ;p<n ; p++) // Calculate Magnetisation
		{
		for (int q=0 ;q<n ; q++)
		{
		M+=vec[p][q];}}

		M=M/(N*N);

//		ans << beta << '\t' << E <<'\t'<< M <<'\t'<< C<<'\t'<<X<<endl;


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

