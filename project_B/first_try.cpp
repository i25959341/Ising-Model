#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_math.h>

using namespace std;

int main()
{
	ofstream outfile("N=10.txt");
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,unsigned long(10));
	
	int N=10, c=0, a;
	double E,M,S,dE,dS,dEc=0.0, dEr=0.0,avgE=0.0,avgE2=0.0,avgS=0.0,avgS2=0.0,C,X,s,x;
	E=-2*N*N;		 // set initial Energy
	S=N*N;			 // set initial Spin

	int vec[100];
	for (int i=0; i<N*N; i++)			 // set array of N*N
	{
		vec[i] = 1;
	}
	for(double b=0.0; b<=1.0 ; b+=0.01)
	{ 
		M=0.0;
		a=0;
		
		for(int T=1;T<1000;T++) //run until equilibrium reach
		{
			s = gsl_rng_uniform (r,N*N);
			
			vec[s]=(-1)*vec[s];
			
			if(s%N==0)              //  Consider the first column
			{
				dEc = -(vec[s]*vec[s+N-1]+vec[s]*vec[s+1]);
			}
			else if (s%N==(N-1))	// Consider the last column
			{
				dEc = -(vec[s]*vec[s-N+1]+vec[s]*vec[s-1]);
			}
			else					//other columns
			{
				dEc = -(vec[s]*vec[s+1]+vec[s]*vec[s-1]);
			}

			if (s<(N-1))				//Consider the top row
			{
				dEr = -(vec[s]*vec[s+N]+vec[s]*vec[s+(N*N-N)]);
			}
			else if (s>=(N*N-N))		 //consider the bottom row
			{
				dEr = -(vec[s]*vec[s-N]+vec[s]*vec[s-(N*N-N)]);
			}
			else					//other rows
			{
				dEr = -(vec[s]*vec[s+N]+vec[s]*vec[s-N]);
			}
			dE = dEc+dEr;					 //total energy contributed from 4 neighbouring site
			dS = 2.0*vec[s];
			
			if(dE>0)
			{
				x=gsl_rng_uniform (r);
				if ( x>= exp(-dE*b) )
				{ 
					vec[s]=(-1)*vec[s];		//spin not flipped, start the loop again
					dE=0;
					dS=0;
				}
			}
			E += dE;
			S += dS;
		}

		for(int t=1;t<1000;t++)		//equilibrium reached
		{
			s = gsl_rng_uniform (r,N*N);
			
			vec[s]=(-1)*vec[s];
			
			if(s%N==0)              //  Consider the first column
			{
				dEc = -(vec[s]*vec[s+N-1]+vec[s]*vec[s+1]);
			}
			else if (s%N==(N-1))	// Consider the last column
			{
				dEc = -(vec[s]*vec[s-N+1]+vec[s]*vec[s-1]);
			}
			else					//other columns
			{
				dEc = -(vec[s]*vec[s+1]+vec[s]*vec[s-1]);
			}

			if (s<(N-1))				//Consider the top row
			{
				dEr = -(vec[s]*vec[s+N]+vec[s]*vec[s+(N*N-N)]);
			}
			else if (s>=(N*N-N))		 //consider the bottom row
			{
				dEr = -(vec[s]*vec[s-N]+vec[s]*vec[s-(N*N-N)]);
			}
			else					//other rows
			{
				dEr = -(vec[s]*vec[s+N]+vec[s]*vec[s-N]);
			}
			dE = dEc+dEr;					 //total energy contributed from 4 neighbouring site
			dS = 2.0*vec[s];
			
			if(dE>0)
			{
				x=gsl_rng_uniform (r);
				if ( x>= exp(-dE*b) )
				{ 
					vec[s]=(-1)*vec[s];		//spin not flipped, start the loop again
					dE=0;
					dS=0;
				}
			}
			E += dE;
			S += dS;

			avgE = (avgE*a+E)/(a+1); //Average Energy <E>
			avgE2= (avgE2*a+(E*E))/(s+1); //Average Energy <E^2>
			
			avgS = (avgS*a+S)/(a+1); //Average Spin <S>
			avgS2= (avgS2*)+(S*S))/(a+1); //Average Spin <S^2>

			C = (avgE2-avgE*avgE)*b*b/(N*N) ; 
			X = (avgS2-avgS*avgS)*b/(N*N); 
			a++;
		}
		for (int c=0; c<N*N; i++)			 
		{
			M += vec[c]/(N*N);
		}
	}
return 0;
}
