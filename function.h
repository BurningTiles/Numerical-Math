/**
 * Author  : BurningTiles
 * Created : 2020-09-07 17:49:16
 * Link    : GitHub.com/BurningTiles
 * Header  : Function
**/

/* 
 * Change the functions as it is required. 
 * Be careful before change it and follow the syntax else it can produce wrong result.
 */


#ifndef FUNCTION_H
#define FUNCTION_H

// include math header for complex operation
// for operation with e you need to define it.
#include <cmath>

// Factorial function.
inline int fact(int x){
	return (x<1)? 1 : x * fact(x-1);
}

// Creating mod function used in while loop testing.
// It will simply find the difference between numbers and output it in positive value number.
inline double mod(double in, double out) {
	double temp = in - out;
	return (temp < 0)? -temp : temp;
}

// Change the function as you want. it is function used in many headers.
// This function is highly required in bisection and these type of headers.
inline double f(double x){
	return (x*x-97);
}

// It is differentiated function of above function.
// It is required for newtonraphson method.
inline double f1(double x){
	return (2*x);
}

// Function of x and y.
// Used in methods like euler, RKmethod and other.
inline double fxy(double x, double y){
	return (x+y);
}

#endif
