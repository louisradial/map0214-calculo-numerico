#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Simpson's rule integration
double test_function(double x) {
    return 5 - 5*pow(x,4);
}

float test_function_float(float x) {
    return 5 - 5*powf(x,4);
}

double simpson_integration(double (*f)(double), double x0, double x1, 
                          int log2N) {
    int intervals = 1 << log2N; // 2^p with p < 31
    double h = (x1 - x0)/intervals;
    double x = x0 + h;
    double integral = f(x0) + f(x1);
    for (int i = 0; i < intervals; i += 2) {
        integral += 4 * f(x) + 2 * f(x + h);
        x += 2*h;
    }
    return integral * h/3;
}

float simpson_integration_float(float (*f)(float), float x0, float x1,
                                int log2N) {
    int intervals = 1 << log2N; // 2^p with p < 31
    float h = (x1 - x0)/intervals;
    float x = x0 + h;
    float integral = f(x0) + f(x1);
    for (int i = 1; i < intervals; i += 2) {
        integral += 4 * f(x) + 2*f(x + h);
        x += 2*h;
    }
    return integral * h/3;
}

void simpson_integration_study(double (*f)(double), float (*f_float)(float),
                               double x0, double x1, double value, int p,
                               char* filename) {
    FILE* fp = fopen(filename, "w");
    int n = 2;
    double err, err_float, result;
    float result_float;
    for (int i = 1; i <= p; ++i) {
        result = simpson_integration(f, x0, x1, i);
        err = fabs(result - value);
        result_float = simpson_integration_float(f_float, x0, x1, i);
        err_float = fabsl(result_float - value);
        fprintf(fp, "%2d,%8d,%.12f,%.12lf,%16.12lf,%.12lf,%.12lf,%16.12lf\n",
               i, n,
               result_float, err_float, log2(err_float),
               result, err, log2(err)
               );
        n <<= 1;
    }
    fclose(fp);
    return;
}

// Trapezoidal rule integration

#define PI acos(-1.0)

double integrand(double t, double k) {
    return pow(1 - k*pow(sin(t), 2), -0.5);
}

double pendulum_period(double t0, int n) {
    double x0 = 0;
    double x1 = PI/2;
    double h = (x1 - x0)/n;
    double x = h;
    double k = (1 - cos(t0))/2;
    double integral = integrand(x0, k) + integrand(x1, k);
    for (int i = 0; i < n; ++i) {
        integral += 2 * integrand(x, k);
        x += h;
    }
    return (2 / PI) * integral * h/2;
}

void pendulum_study(int n, int divisions, char* filename) {
    FILE* fp = fopen(filename, "w");
    double t0;
    double h = PI/(n+1);
    double t;
    for (int i = 0; i < n; ++i) {
        t0 += h;
        t = pendulum_period(t0, divisions);
        fprintf(fp, "%2d, %6.2lf, %16.12lf\n", i+1, 180*t0/PI, t);
    }
    fclose(fp);
    return;
}

// Random number generation (LCG)

#define RANDOM_SEED 8992822
#define RANDOM_A 1664525
#define RANDOM_B 1013904223
#define RANDOM_M 2147483647

typedef struct {
    unsigned long long a;
    unsigned long long b;
    unsigned long long m;
    unsigned long long seed;
} Random;

Random* initialize_lgc(unsigned long long a, unsigned long long b,
                       unsigned long long m, unsigned long long seed) {
    Random* r = (Random*)malloc(sizeof(Random));
    r->a = a;
    r->b = b;
    r->m = m;
    r->seed = seed;
    return r;
}

double random_lcg(Random* r) {
    double x = (double)r->seed/r->m;
    r->seed = (r->a*r->seed + r->b) % r->m;
    return x;
}

void random_test(unsigned long long a, unsigned long long b,
                 unsigned long long m, unsigned long long seed, int n) {
    printf("testing lgc\n");
    Random* r = initialize_lgc(a,b,m,seed);
    double x;
    for (int i = 0; i <= n; ++i) {
        seed = r->seed;
        x = random_lcg(r);
        printf("%2d %10lld %.12lf\n", i, seed, x);
    }
    free(r);
    return;
}

// Monte Carlo integration

double parabola(double x) {
    return 1 - x*x;
}

double integration_montecarlo_trial(double(*f)(double), double x0, double x1,
                                    Random *r, int points) {
    // assumption: f(x) takes values between 0 and 1 for x0 <= x <= x1
    double x,y;
    int inside = 0;
    for (int i = 0; i < points; ++i) {
        x = x0 + (x1 - x0)*random_lcg(r);
        y = random_lcg(r);
        if (y <= f(x)) {
            inside += 1;
        }
    }
    return (x1 - x0)*inside/points;
}

#define POINTS 100

void parabola_quadrature(int max_attempts, char* filename) {
    FILE* fp = fopen(filename, "w");
    Random *r = initialize_lgc(RANDOM_A, RANDOM_B, RANDOM_M, RANDOM_SEED);
    double *trial = (double*)malloc(max_attempts*sizeof(double));
    double integral, mean, variance, std, mean_std;
    for (int attempts = 2; attempts <= max_attempts; attempts *= 2) {
        mean = 0;
        for (int i = 0; i < attempts; ++i) {
            integral = integration_montecarlo_trial(parabola, 0, 1, r, POINTS);
            mean += integral;
            trial[i] = integral;
        }
        mean /= attempts;
        variance = 0;
        for (int i = 0; i < attempts; ++i) {
            variance += pow(trial[i] - mean,2);
        }
        variance = variance/(attempts - 1);
        std = sqrt(variance);
        mean_std = std/sqrt(attempts);
        fprintf(fp, "%6d, %.12lf, %.12lf, %.12lf\n",
                attempts, mean, std, mean_std);
    }
    free(r);
    free(trial);
    fclose(fp);
}

int main() {
    // 1
    simpson_integration_study(test_function, test_function_float,
                              0, 1, 4.0, 25, "1.csv");
    // 2
    pendulum_study(20, 1 << 16, "2.csv");
    // 3
    random_test(RANDOM_A, RANDOM_B, RANDOM_M, RANDOM_SEED, 10);
    printf("\nparabola quadrature\n");
    parabola_quadrature(1 << 17, "3.csv");
    return 0;
}
