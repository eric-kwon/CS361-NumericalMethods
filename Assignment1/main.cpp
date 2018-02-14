/*
 *      NAME : ERIC KWON
 *      CSCI361 HW #1
 */

#include <iostream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <vector>
using namespace std;

// Problem 1.1.1 - GCD
long gcd_euclid(long a, long b) {

    if (a < b) 
        return gcd_euclid(b,a);

    // Validation
    // If either input is less than 0 return 0
    if (a < 0 || b < 0)
        return 0;

    // If either input is 0 return the other input
    if (a == 0 || b == 0)
        return abs(a - b);

    long c = a % b;

    if (c == 0) 
        return b;
    else if (c == 1) 
        return 1;
    else 
        return gcd_euclid(b,c);
}





// Problem 1.1.2 - LCM
long LCM(long a, long b) {

    // Validation
    // If either input is less than or equal to 0 return 0
    if (a <= 0 || b <= 0)
        return 0;

    long gcd = gcd_euclid(a, b);
    long mult = a * b;

    return mult / gcd;
}




// Problem 1.3 - Horner's Rule + Taylor Series
double exp_sum(double x) {

    // Parameters
    int n_max = 30;

    // Sum at n_max
    double sum = 1 + (x / n_max);

    // Loop for applied Horner's rule onto Taylor series
    for (int i = n_max-1 ; i > 0 ; i--) {
        sum = 1 + (x / i) * sum;
    }

    return sum;
}





// Problem 1.4 - Convex Hull
class Point {
    public:
        // Class variables
        double x;
        double y;
        double r;
        double theta;

        // Compute theta and r
        Point(double x1, double y1) : x(x1), y(y1) {
            r = sqrt(x * x + y * y);
            if (r > 0.0) 
                theta = atan2(y, x);
            else
                theta = 0;
        }
};

// Comparison function
bool vcomp(const Point &a, const Point &b) {
    if (a.theta < b.theta)
        return true;
    else if (a.theta > b.theta)
        return false;

    return (a.r < b.r);
}

// Convex hull function
void convex_hull_func(int n,                              // Number of points
                        const vector<double> &x,          // Array of x-coordinate
                        const vector<double> &y,          // Array of y-coordinate
                        vector<Point> &convex_hull)       // Array of point output
{
    // Tolerance
    const double tol = 1.0e-14;     

    // Clear the array                    
    convex_hull.clear();                                

    // Validation - n is less than or equal to zero - return
    if (n <= 0) return;                           

    // Testing for cases where there are less than or equal to 3 points
    // Populate the convex hull in such a case and return
    if (n <= 3) {                                      
        for (int i = 0 ; i < n ; i++) {
            convex_hull.push_back(Point(x[i], y[i]));
        }
        convex_hull.push_back(Point(x[0], y[0]));
        return;
    }

    // Finding the minimum value y_min, x_min, and i_min
    int imin = 0;
    double xmin = x[0];
    double ymin = y[0];

    // Searching loop
    for (int i = 0 ; i < n ; i++) {
        if (ymin > y[i]) {
            imin = i;
            xmin = x[i];
            ymin = y[i];
        }
        else if (ymin == y[i] && xmin > x[i]) {
            imin = i;
            xmin = x[i];
            ymin = y[i];
        }
        else
            continue;
    }

    // Placeholder for the cooridnate vector
    vector<Point> v;

    // Coordinate vector population
    for (int iv = 0 ; iv < n ; iv++) {
        double vx = x[iv] - xmin;
        double vy = y[iv] - ymin;
        if (iv == imin) {
            vx = 0;
            vy = 0;
        }
        v.push_back(Point(vx, vy));
    }

    // Sorting the vector
    sort(v.begin(), v.end(), vcomp);

    // Adding x_min and y_min to the sorted vector
    for (int ip = 0 ; ip < n ; ip++) {
        v[ip].x += xmin;
        v[ip].y += ymin;
    }

    // Referenced vector
    vector<Point> &p = v;

    // Stack initialization
    convex_hull.push_back(p[0]);
    convex_hull.push_back(p[1]);

    // Size parameter
    int i = convex_hull.size();

    // Loop through the points and updating the stack
    while (i < n) {
        int last = convex_hull.size() - 1;
        int second_last = last - 1;

        double ux = convex_hull[last].x - convex_hull[second_last].x;
        double uy = convex_hull[last].y - convex_hull[second_last].y;
        double wx = p[i].x - convex_hull[last].x;
        double wy = p[i].y - convex_hull[last].y;

        double cross_product = ux * wy - uy * wx;

        if (cross_product > tol) {
            convex_hull.push_back(p[i]);
            i++;
        }
        else if (fabs(cross_product) <= tol) {
            convex_hull.pop_back();
            convex_hull.push_back(p[i]);
            i++;
        }
        else {
            convex_hull.pop_back();
        }
    }

    // Close out the polygon
    convex_hull.push_back(p[0]);
    
    // Exit the function
    return;
}                     





int main() {

    // Test code here

    return 0;
}