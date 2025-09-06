#include <stdio.h>
#include <math.h>

#define MAX_ITERATIONS 100

// returns sign of a double. equivalent to (x > 0) ? 1 : ((x < 0) ? -1 : 0)
int sign(double x) {
    return (x > 0) - (x < 0);
}

// bisection method
double bisection(double (*function)(double), double initial, double final,
                 double precision) {
    double x, x1, x2, fx1, fx, e; // variables for the method and printing
    int n = 0; // iteration counter
    x1 = initial;
    x2 = final;
    do {
        x = (x1 + x2)/2; // midpoint xm
        fx1 = function(x1); // f(x1)
        fx = function(x); // f(xm)
        e = fabs(x2 - x1); // absolute difference between endpoints
        printf("n=%2d x1=%+lf x2=%+lf xm=%+lf f(x1)=%+lf f(xm)=%+lf e=%lf\n",
               n, x1, x2, x, fx1, fx, e);
        // update endpoints based on the sign of f(xm) and f(x1)
        if (sign(fx) == sign(fx1)) {
            x1 = x;
        }
        else {
            x2 = x;
        }
        n++; // increment iteration counter
        // halt condition due to precision or number of iterations
    } while (e > precision && n < MAX_ITERATIONS);
    return x;
}

// newton raphson method
double newton(double (*function)(double), double (*derivative)(double),
              double guess, double precision) {
    int n = 0; // iteration counter
    // variables for the method and printing
    double x = guess;
    double subtract, fx, fpx, abs_diff;
    do {
        fx = function(x); // f(x)
        fpx = derivative(x); // f'(x)
        subtract = function(x)/derivative(x); // f(x)/f'(x)
        abs_diff = fabs(subtract); // absolute difference between iterations
        printf("n=%d x=%+.12lf f(x)=%+.12lf f'(x)=%+.12lf e=%.12lf\n",
               n, x, fx, fpx, abs_diff);
        x -= subtract; // xnext = xcurr - f(x)/f'(x)
        n++; // increment iteration counter
        // halt condition due to precision or number of iterations
    } while (abs_diff > precision && n < MAX_ITERATIONS);
    return x;
}

// secant method
double secant(double (*function)(double), double x0, double x1, 
              double precision) {
    int n = 0; // iteration counter
    // variables for the method and printing
    double xprev = x0;
    double xcurr = x1;
    double xnext, slope, fxcurr, fxprev, subtract, e;
    printf("n=%2d x=%+.12lf f(x)=%+.12lf e=%.12lf\n", n++, x0, function(x0), fabs(x1 - x0));
    do {
        fxcurr = function(xcurr); // f(xcurr)
        fxprev = function(xprev); // f(xprev)
        slope = (fxcurr - fxprev)/(xcurr - xprev); // estimates f'(xcurr)
        subtract = fxcurr/slope; // estimates f(xcurr)/f'(xcurr)
        e = fabs(subtract); // absolute difference between iterations
        printf("n=%2d x=%+.12lf f(x)=%+.12lf e=%.12lf\n", n, xcurr, fxcurr, e);
        // computes xnext and updates xcurr and xprev
        xnext = xcurr - subtract;
        xprev = xcurr;
        xcurr = xnext;
        n++;
        // halt condition due to precision or number of iterations
    } while (e > precision && n < MAX_ITERATIONS);
    return xnext;
}

// function for items a and b
double a(double x) {
    return pow(x,3) + x + cos(x);
}

// its derivative for NR method
double aprime(double x) {
    return 3*pow(x, 2) + 1 - sin(x);
}

// function for item c
#define V0 9
double c(double x) {
    return tan(sqrt(x + V0)) - sqrt(-x/(x + V0));
}

int main() {
    // items a and b
    double bisection_root = bisection(a,-1, 1, 1e-6);
    printf("initial estimate with bisection method: %.6lf -> %.6lf\n",
           bisection_root, a(bisection_root));
    double nr_root = newton(a, aprime, -1, 1e-12);
    printf("root found by Newton-Raphson method: %.12lf -> %.12lf\n",
           nr_root, a(nr_root));

    // item c
    double energy = secant(c, -8, -7, 1e-12);
    printf("root found by secant method: %.12lf -> %.12lf", energy, c(energy));
    return 0;
}
