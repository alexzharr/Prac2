#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-10
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double f1(double x, double p, double alpha);
double f2(double x, double p, double alpha);
double find_ddu0(double* x, int n, double y0, double y1, double dy0);
void RungeKutta(double* t, double* x, double *p, int n, double x0, double p0, double alpha);

double f1(double x, double p, double alpha)
{
    double t =  M_PI / (4 + alpha * x * x);
    return p + sin(t);
}

double f2(double x, double p, double alpha)
{
    double t =  M_PI / (4 + alpha * x * x);
    return (2 * alpha * p * t * cos(t)) / (4 + alpha * x * x) + 0.5 * x;
}

/*double find_p0(double *x, int n, double y0, double y1, double dy0)
{
    double f_cent, f_right;
    double right_diff;
    double result = 0;
    double *u;

    u = (double *)malloc(n * sizeof(double));
    
    RungeKutta(x, u, n, y0, dy0, result);
    f_cent = u[n - 1];
    RungeKutta(x, u, n, y0, dy0, result + 1);
    f_right = u[n - 1];
    right_diff = (f_right - f_cent);
    result += (y1 - f_cent) / right_diff;  

    free(u);
    return result;
}*/

void RungeKutta(double* t, double* x, double *p, int n, double x0, double p0, double alpha)
{
    double h = 1. / (n - 1);
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];

    x[0] = x0;
    p[0] = p0;

    for (int i = 0; i < n - 1; i++)
    {
        k1[0] = f1(x[i], p[i], alpha);
        k1[1] = f2(x[i], p[i], alpha);
        k2[0] = f1(x[i] + (h / 2) * k1[0], p[i] + (h / 2) * k1[1], alpha);
        k2[1] = f2(x[i] + (h / 2) * k1[0], p[i] + (h / 2) * k1[1], alpha); 
        k3[0] = f1(x[i] + (h / 2) * k2[0], p[i] + (h / 2) * k2[1], alpha);
        k3[1] = f2(x[i] + (h / 2) * k2[0], p[i] + (h / 2) * k2[1], alpha);
        k4[0] = f1(x[i] + h * k3[0], p[i] + h * k3[1], alpha);
        k4[1] = f2(x[i] + h * k3[0], p[i] + h * k3[1], alpha);

        x[i + 1] = x[i] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        p[i + 1] = p[i] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);   
            
    }
    free(x);
    free(p);
}

int main()
{
    int n;
    double h;
    double alpha;
    double *t, *x, *p;

    printf("Enter alpha: ");
    scanf("%lf", &alpha);

    h = 1. / (n - 1);
    t = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));
    p = (double *)malloc(n * sizeof(double));
    for (int i = 1; i < n; i++)
        t[i] = t[i - 1] + h;

    RungeKutta(t, x, p, n, 0, 0, alpha);

    free(t);
    free(x);
    free(p);
    return 0;
}