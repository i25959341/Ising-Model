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
int mod(int a,int b) // Modulus calculator
{ 
int m;
if (a<0)
{
m=b +(a % b);
}
else
{
m=a%b;
}
return m;
} 
int main ()
{
const gsl_rng_type * T;
gsl_rng * r;
T = gsl_rng_default;
r = gsl_rng_alloc (T);
gsl_rng_set(r,unsigned long(5));
int n,m,y,z,k=1;
double x,E,dE,aE=0,aE2=0,M,C,X,S,dS,aS=0,aS2=0;
ofstream ans("ans.txt");
n=100;
vector<vector<double>> vec(n+1,vector<double>(n+1)); // Construct a new(n times n)square grid 
E=-20000; // Intial Energy
S=10000; //Intial Spin
for (int r=0 ;r<=n-1 ; r++) // Initial configuration
{
for (int s=0 ;s<=n-1 ; s++) 
{
vec[r][s]=1;
}
}
for(double B=0.1; B<2.01 ; B+=0.01) // Loop different Temperature
{ 
M=0.0;
k=0;
cout<<B<<endl;
for(int T=1;T<1000000;T++) //Make sure equilibrium reach
{
y=gsl_rng_uniform_int(r,n);
z=gsl_rng_uniform_int(r,n);
vec[y][z]=(-1)*vec[y][z];
dE=-2*(vec[y][z]*vec[y][(z+1)%(n)]+vec[y][z]*vec[y][mod((z-1),n)]+vec[y][z]*vec[mod(y-1,n)][z]+vec[y][z]*vec[(y+1)%(n)][z]);
dS=2*vec[y][z];
if(dE>0)
{
x=gsl_rng_uniform (r);
if ( x>= exp(-dE*B) )
{ 
vec[y][z]=(-1)*vec[y][z];
dE=0;
dS=0;
}
}
E+=dE;
S+=dS;
}
for (int t=1; t<200000 ; t++) // loop after equilibrium reach
{
y=gsl_rng_uniform_int(r,n);
z=gsl_rng_uniform_int(r,n);
vec[y][z]=(-1)*vec[y][z];
dE=-2*(vec[y][z]*vec[y][(z+1)%(n)]+vec[y][z]*vec[y][mod((z-1),n)]+vec[y][z]*vec[mod(y-1,n)][z]+vec[y][z]*vec[(y+1)%(n)][z]);
dS=2*vec[y][z];
if(dE>0)
{
x=gsl_rng_uniform (r);
if ( x>= exp(-dE*B) )
{ 
vec[y][z]=(-1)*vec[y][z];
dE=0;
dS=0;
}
}
E+=dE;
S+=dS;
aE= (aE*(k)+(E))/(k+1); //Calculate Average Energy <E>
aE2= (aE2*(k)+(E)*(E))/(k+1); //Calculate Average Energy <E^2>
C=(aE2-aE*aE)*B*B/(n*n) ; //Calculate specific heat capacity
aS = (aS*(k)+(S))/(k+1); //Calculate Average Spin <S>
aS2= (aS2*(k)+(S)*(S))/(k+1); //Calculate Average Spin <S^2>
X=(aS2-aS*aS)*B/(n*n); //Calculate susceptibility <E>
k++;
}
for (int p=0 ;p<=n-1 ; p++) // Calculate Magnetisation
{
for (int q=0 ;q<=n-1 ; q++) 
{
M+=vec[p][q];
}
}
ans << 1/B << '\t' << E <<'\t'<< M/(10000.0) <<'\t'<< aS /(10000.0)<< '\t'<< C<<'\t'<<X<<endl;
}
return 0;
}
