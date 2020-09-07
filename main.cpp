/**
 * Author  : BurningTiles
 * Created : 2020-09-07 17:48:20
 * Link    : GitHub.com/BurningTiles
 * Program : Numerical Math Main
**/

// Set the accuracy as you want.
#define accuracy 0.00001

// Header included
#include "nmath.h"

using namespace std;

int main(){
	Newtonraphson object;  // Change object type according to methods. make sure first character is capital.
	return object.main();  // Calling main function of object which will do all necessary operations.
}

/**

Example : 
Function set for x*2-97
(so it will give square root of 97 using Newtonraphson method)

Enter approximated root : 8
Iteration   Approx_Root   Value_of__f(x)
1           10.06250000     4.2539062500
2            9.85112578     0.0446790624
3            9.84885806     0.0000051425
4            9.84885780     0.0000000000


Root of given function is 9.84885780

**/