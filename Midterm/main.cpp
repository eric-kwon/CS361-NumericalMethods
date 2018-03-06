/*
 *      CSCI361 Midterm
 *      DUE: MARCH 2, 2018
 *      NAME: ERIC KWON
 */

#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

// Combinatorics - n choose r
int nCr(int n, int r) {
    
    // Validation - if r is 0 return 1
    if (r == 0) return 1;

    // Since nCr output values are symmetric - use lower bound to reduce overhead
    if (r > n/2) return nCr(n, n-r);

    // Result placeholder
    long result = 1; 

    // Loop to calculate the value
    for(int i = 1 ; i <= r ; i++) {
        result *= n - i + 1;
        result /= i;
    }

    // Return the result
    return result;
}

// Fuss-Catalan function
int FussCatalan(int m, int p, int r) {

    // Validation 1 - if m or p < 0 - return 0
    if (m < 0 || p < 0) return 0;
    
    // Validation 1.1 - if m is 1 - return 1
    if (m == 1) return 1;

    // Placeholder for denominator
    int denom = (m * p) + r;

    // Validation 2 - if (mp + r) < m - return 0
    if (denom < m) return 0;

    // Placeholder for nCr
    long combo = nCr(denom, m);

    // Return result
    return (r * combo) / denom;
}

// Sum of array function
double absum(int n, const double *a, const double *b, double x) {

    // Validation
    // Case 1 - empty array
    // Case 2 - if the last index of array A is less than x
    if (n <= 0 || a[n-1] < x)
        return 0;

    // Placeholder for result
    double result = 0;

    // Index placeholder
    int i = n - 1;

    // Loop and sum
    // Check 1 - check for valid index
    // Check 2 - comparison check
    while (i >= 0 && a[i] > x) {
        result += b[i];
        i--;
    }

    // Return the result
    return result;
}

// Summation for S(x)
double sum_f(double x, int n) {

    // Constants
    const double pi = 4.0 * atan2(1.0, 1.0);
    const double frac = 1.0 / n;

    // Looping result placeholder
    double sum = 0.0;

    // Loop for summation
    for (int i = 0 ; i <= n - 1 ; i++) {
        double sintemp = sin((i * pi) / n);
        sum += cos(x * sintemp);
    }

    // Return result
    return frac * sum;
}

// Summation for S'(x)
double sum_fprime(double x, int n) {

    // Constants
    const double pi = 4.0 * atan2(1.0, 1.0);
    const double frac = -1.0 / n;

    // Looping result placeholder
    double sum = 0.0;

    // Loop for summation
    for (int i = 0 ; i <= n - 1 ; i++) {
        double sintemp = sin((i * pi) / n);
        sum += sin(x * sintemp) * sintemp;
    }

    // Return result
    return frac * sum;
}


// Function to be tested under Newton-Raphson
void func(double x, double &f, double &fprime) {

    //Constants - set N at 20
    const int n = 20;

    // S(x) and S'(x)
    f = sum_f(x, n);
    fprime = sum_fprime(x, n);
}

// Newton Raphson function - return 0 on success, 1 on fail
int root_NR(double tol,
            int max_iter,
            double &x,
            double &function)
{
	// Initialization
    double f = 0.0;
    double fprime = 0.0;
    double x_diff = 0.0;
    int iteration = 1;    
    
    // NR Iteration
    do{
    	// Placeholder for x0
    	double x0 = x;
    	
    	// Function to be tested
    	func(x, f, fprime);
    	
    	// Finding the following x
    	x = x0 - (f / fprime);
    	
    	// Difference to compare tolerance
    	x_diff = fabs(x - x0);
    	
    	// If the difference is below the tolerance - return 0
    	if (x_diff <= tol)
    		return 0;
    	
    	// Save f to function placeholder
    	function = f;
    	
    	// Increment
    	iteration++;
	}while(iteration < max_iter);

    // Return 0 on success
    return 0;
}

int main() {

    // Test codes here
    return 0;
}
