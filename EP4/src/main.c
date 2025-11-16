#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NO_ARGS NULL
typedef char* filename;
typedef double (*ode_function)(double t, double x, double xdot, double* params);

typedef struct {
    ode_function function;
    double* params;
} ode;

double compute(ode *g, double t, double x, double xdot) {
    return g->function(t, x, xdot, g->params);
}

void euler(ode *g, double *ti, double *yi, double *zi, double h) {
    // y' = z, z' = g(t,y,z)
    double z = *zi;
    double y = *yi;
    double t = *ti;
    *yi += h*z;
    *zi += h*compute(g, t, y, z);
    *ti += h;
}

void solve_euler(ode g, double t0, double y0, double z0, double h, double tf,
                 filename output) {
    FILE* fp = fopen(output, "w");
    if (fp == NULL) {
        printf("File error");
        return;
    }
    double t = t0;
    double y = y0;
    double z = z0;
    fprintf(fp, "%.6lf, %.6lf, %.6lf\n", t, y, z);
    while (fabs(t - tf) > h/2) {
        euler(&g, &t, &y, &z, h);
        fprintf(fp, "%.6lf, %.6lf, %.6lf\n", t, y, z);
    }
    fclose(fp);
}

void rk4(ode *g, double *ti, double *yi, double *zi, double h) {
    // y' = z, z' = g(t,y,z)
    double z = *zi;
    double y = *yi;
    double t = *ti;
    double k1y = h * z;
    double k1z = h * compute(g, t, y, z);
    double k2y = h * (z + k1z/2);
    double k2z = h* compute(g, t + h/2, y + k1y/2, z + k1z/2);
    double k3y = h * (z + k2z/2);
    double k3z = h * compute(g, t + h/2, y + k2y/2, z + k2z/2);
    double k4y = h * (z + k3z);
    double k4z = h * compute(g, t + h, y + k3y, z + k3z);
    y += (k1y + 2*k2y + 2*k3y + k4y)/6;
    z += (k1z + 2*k2z + 2*k3z + k4z)/6;
    t += h;
    *ti = t;
    *yi = y;
    *zi = z;
    return;
}

void solve_rk4(ode g, double t0, double y0, double z0, double hi, double ti,
               double hf, double tf, filename output) {
    FILE* fp = fopen(output, "w");
    if (fp == NULL) {
        printf("File error");
        return;
    }
    double t = t0;
    double y = y0;
    double z = z0;
    // evolve without printing, used to skip transient
    while (fabs(t - ti) > hi/2) {
        rk4(&g, &t, &y, &z, hi);
    }
    fprintf(fp, "%.6lf, %.6lf, %.6lf\n", t, y, z);
    while (fabs(t - tf) > hf/2) {
        rk4(&g, &t, &y, &z, hf);
        fprintf(fp, "%.6lf, %.6lf, %.6lf\n", t, y, z);
    }
    fclose(fp);
}

void poincaré(ode g, double t0, double x0, double v0, double period,
               int t_periods, int s_periods, double f, FILE* file) {
    double t = t0;
    double x = x0;
    double v = v0;
    double t_step = period*1e-2;
    double s_step = period*1e-3;
    // transient evolution computes t_periods periods
    for (int i = 0; i < t_periods*100; i++) {
        rk4(&g, &t, &x, &v, t_step);
    }
    // steady state evolution computes s_periods periods
    for (int i = 0; i < s_periods; i++) {
        // print after each period
        for (int j = 0; j < 1000; j++) {
            rk4(&g, &t, &x, &v, s_step);
        }
        fprintf(file, "%lf, %lf, %lf, %lf\n", t, x, v, f);
    }
}

double ex1(double t, double y, double ydot, double* params) {
    return ydot + y - pow(t,3) - 4*pow(t, 2) + 4*t + 2;
}

#define PI acos(-1.0)

double double_well(double t, double x, double xdot, double* params) {
    double g = params[0]; // 2γ
    double f = params[1]; // F
    double w = params[2]; // ω
    return 0.5 * x * (1 - x*x) - g*xdot + f*cos(w*t);
}

#define DOUBLE_WELL(g,f,w) (ode){double_well, (double[3]){g,f,w}}

int main() {
    // 1
    ode ode_1 = {ex1, NO_ARGS};
    solve_rk4(ode_1, 0,0,0,1e-2, 0, 1e-2, 6, "rk4.csv");
    solve_euler(ode_1, 0,0,0,1e-2, 6, "euler.csv");
    // 2a change initial conditions
    double t0 = 0;
    double x0 = -1;
    double v0 = 1;
    double h = 1e-3;
    double tf = 40;
    solve_rk4(DOUBLE_WELL(0,0,0), t0, x0, 0.1, h, t0, h, tf, "2a01.csv");
    solve_rk4(DOUBLE_WELL(0,0,0), t0, x0, 0.5, h, t0, h, tf*3, "2a05.csv");
    solve_rk4(DOUBLE_WELL(0,0,0), t0, x0,  v0, h, t0, h, tf, "2a10.csv");
    // 2b change γ
    double g = 0.25;
    solve_rk4(DOUBLE_WELL(0.80,0,0), t0, x0, v0, h, t0, h, tf, "2b80.csv");
    solve_rk4(DOUBLE_WELL(g,0,0), t0, x0, v0, h, t0, h, tf, "2b25.csv");
    // 2c γ = 0.25, ω = 0.95, change F
    double w = 0.95;
    double hi = h/10;
    double ti = 400;
    tf = 500;
    solve_rk4(DOUBLE_WELL(g,0.190,w), t0, x0, v0, hi, ti, h, tf, "2c190.csv");
    solve_rk4(DOUBLE_WELL(g,0.203,w), t0, x0, v0, hi, ti, h, tf, "2c203.csv");
    solve_rk4(DOUBLE_WELL(g,0.240,w), t0, x0, v0, hi, ti*3, h, tf*3, "2c240.csv");
    solve_rk4(DOUBLE_WELL(g,0.330,w), t0, x0, v0, hi, ti, h, tf, "2c330.csv");
    solve_rk4(DOUBLE_WELL(g,0.600,w), t0, x0, v0, hi, ti, h, tf, "2c600.csv");
    // 3 poincaré section
    FILE* poincaré_section = fopen("poincaré_section.csv", "w");
    ode duffing = DOUBLE_WELL(g,0,w);
    double fmax = 0.7;
    double df = 5e-4;
    double period = 2*PI/w;
    for (double f = 0; fabs(f - fmax) > df/2; f += df) {
        duffing.params[1] = f;
        poincaré(duffing, t0, x0, v0, period, 2000, 100, f, poincaré_section);
    }
    fclose(poincaré_section);
    // feigenbaum constant
    FILE* feigenbaum = fopen("feigenbaum.csv", "w");
    df = 5e-5;
    fmax = 0.205;
    for (double f = 0.18; fabs(f - 0.207) > df/2; f += df) {
        duffing.params[1] = f;
        poincaré(duffing, t0, x0, v0, period, 2000, 100, f, feigenbaum);
    }
    fclose(feigenbaum);
    // 4 poincaré map
    FILE* poincaré_map = fopen("poincaré_map.csv", "w");
    duffing.params[1]=0.24;
    poincaré(duffing, t0, x0, v0, period, 2000, 20000, 0.24, poincaré_map);
    fclose(poincaré_map);
    return 0;
}
