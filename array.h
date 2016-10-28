//define a class for the array of spins

#ifndef ARRAY_H
#define ARRAY_H
#include "array.h"

class array{
private:
	//list of all NXN members
	double list[400];

public:

//default constructor
	initialisearray(){
		for(int i=0; i<400; i++){
		list[i]=1.0;
	}
	}

//Access method for members
	double get (int i, int j) //by x, y coords
	{
		int n=(i-1)*20+(j-1);
		return list[n];
	}

// Modifier method
	void set (int i, int j, double s){

		int n=(i-1)*20+(j-1);
		list[n]=s;

	}

}

}

#endif
