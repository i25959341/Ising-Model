//define class for 20x20 array of spins

#ifndef SPINARRAY_H
#define SPINARRAY_H
#include "spinarray.h"

class spinarray
{
private:
	//list of all 400 members (0101,0102,0103, ... ,0120;0201,0202,0203, ... ,0220; ... ; 2001,2002,2003, ... ,2020)
	double list[400]; 

public:

//default constructor
	spinarray()
	{
		for (int i=0;i<400;i++){
			list[i]=1.0; //initialise all members to 1.0
		}
	}

//Access method for members
	double get (int i, int j) //by x, y coords
	{
		int n=(i-1)*20+(j-1);
		return list[n];
	}

//Modifier method for members
	void set (int i, int j, double s) //by x, y coords
	{
		int n=(i-1)*20+(j-1);
		list[n]=s;
	}

}; 

#endif


