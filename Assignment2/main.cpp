/*
 *      NAME : ERIC KWON
 *      CSCI361 HW #2
 */

#include <iostream>
#include <math.h>
using namespace std;

// Function to find Cumulative Normal Distribution
// Returns N(x) with given x
double cum_norm(double x) {
    const double root = sqrt(0.5);
    return 0.5 * (1.0 + erf(x * root));
}

// Function to be tested under bisection
// double func(double x) {
//     return x * x;
// }

// Function testing using Cumulative Normal - Bisection
double func(double x) {
    return cum_norm(x);
}

// Use int return type - since the function may not converge
// Return 0 on success, otherwise return 1
int root_bisection(double target,
                   double tol_f,        // Function tolerance
                   double tol_x,        // X-axis tolerance
                   int max_iter,
                   double x_low,
                   double x_high,
                   double &x,
                   long &num_iter)      // Upper limit
{
    // Initialization
    x = 0;
    num_iter = 0;

    // Calculation
    double y_low = func(x_low);
    double diff_y_low = y_low - target;

    // Validation
    // Value already within tolerance - succeed
    if (abs(diff_y_low) <= tol_f) {
        x = x_low;
        return 0;
    }

    // Calculation
    double y_high = func(x_high);
    double diff_y_high = y_high - target;

    // Validation
    // Value already within tolerance - succeed
    if (abs(diff_y_high) <= tol_f) {
        x = x_high;
        return 0;
    }

    // Valiation
    // Low and high have sames signs - fail
    if (diff_y_low * diff_y_high > 0.0) {
        x = 0;
        return 1;
    }

    // Bisection Iteration
    for (num_iter = 1 ; num_iter < max_iter ; num_iter++) {

        // Locate the midpoint and perform calculation
        x = (x_low + x_high) / 2.0;
        double y = func(x);
        double diff_y = y - target;

        // Validation
        // Function value within tolerance - succeed
        if (abs(diff_y) <= tol_f) {
            return 0;
        }

        // Check signs and reassign midpoint
        // Target lying on the upper half
        if (diff_y * diff_y_low > 0.0) {
            x_low = x;
        }
        // Target lying on the lower half
        else {
            x_high = x;
        }

        // Validation
        // X-axis converges - succeed
        if (abs(x_high - x_low) <= tol_x) {
            return 0;
        }
    }

    // If code reaches here - failed
    x = 0;
    num_iter = max_iter;
    return 1;
}

// Function to be tested under Newton-Raphson
// void func(double x, double &f, double &fprime) {
//     f = x * x;
//     fprime = 2.0 * x;
// }

// Function testing using Cumulative Normal - Newton-Raphson
void func(double x, double &f, double &fprime) {
    const double pi = 4.0 * atan2(1.0, 1.0);
    f = cum_norm(x);
    fprime = exp(-0.5 * x * x) / sqrt(2.0 * pi);
}

// Use int return type - since the function may not converge
// Return 0 on success, otherwise return 1
int root_NR(double target,
            double tol_f,               // Function tolerance
            double tol_x,               // X-axis tolerance
            int max_iter,
            double x0,
            double &x,
            long &num_iter)             // Upper limit
{
    // Initialization
    const double tol_fprime = 1.0e-12;
    double f = 0;
    double f_prime = 0;

    // Parameter to begin iteration
    x = x0;

    // NR Iteration
    for (num_iter = 1 ; num_iter < max_iter ; num_iter++) {
        
        // Function to be tested
        func(x, f, f_prime);

        // Calculate the difference
        double diff_f = f - target;

        // Validation
        // Function value within tolerance - succeed
        if (abs(diff_f) <= tol_f) {
            return 0;
        }

        // Validation
        // F' within tolerance - fail
        if (abs(f_prime) <= tol_fprime) {
            x = 0;
            return 1;
        }

        // Calculation for delta
        double delta_x = diff_f / f_prime;

        // Validation
        // Delta within x-axis tolerance - succeed
        if (abs(delta_x) <= tol_x) {
            return 0;
        }

        // Update the value of x
        x -= delta_x;
    }

    // If code reaches here - failed
    x = 0;
    num_iter = max_iter;
    return 1;
}

int main() {
    // Test Function Here
}