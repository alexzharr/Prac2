#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define X0 0
#define P1 1
#define EPS 1e-6
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef struct Point
{
    double t;
    double x;
    double p;
    struct Point *next;
} Point;

void pushBack(Point *last, double t_, double x_, double p_) {
    Point *pt = (Point*)malloc(sizeof(Point));
    pt->t = t_;
    pt->x = x_;
    pt->p = p_;
    last->next = pt;
}


double f1(double x, double p, double alpha);
double f2(double x, double p, double alpha);
double find_p0(double *t, double *x, double *p, int n, double x0, double alpha);
double RungeKutta(double* t, double* x, double *p, int n, double x0, double p0, double alpha);

/*double f1(double x, double p, double alpha)
{
    return p + x;
}

double f2(double x, double p, double alpha)
{
    return x;
}*/

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

double find_p0(double *t, double *x, double *p, int n, double x0, double alpha)
{
    double f_left, f_right;
    double df = 1;
    double p0 = 0;
    double p1 = 0;
    int a;

    do
    {
        p0 = p0 - p1 / df;
        p1 = RungeKutta(t, x, p, n, X0, p0, alpha);
        f_left = RungeKutta(t, x, p, n, X0, p0 - EPS, alpha);
        f_right = RungeKutta(t, x, p, n, X0, p0 + EPS, alpha);
        df = (f_right - f_left) / (2 * EPS);
    } while (fabs(p1) > EPS);


    return p0;
}

double RungeKutta(double* t, double* x, double *p, int n, double x0, double p0, double alpha)
{
    double h;
    double loc_err;
    double E[2];
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];

    t[0] = 0;
    x[0] = x0;
    p[0] = p0;

    for (int i = 0; i < n - 1; i++)
    {
        h = 1 - t[i];
        do
        {
            h /= 2;
            k1[0] = f1(x[i], p[i], alpha);
            k1[1] = f2(x[i], p[i], alpha);
            k2[0] = f1(x[i] + (h / 2) * k1[0], p[i] + (h / 2) * k1[1], alpha);
            k2[1] = f2(x[i] + (h / 2) * k1[0], p[i] + (h / 2) * k1[1], alpha); 
            k3[0] = f1(x[i] + (h / 2) * k2[0], p[i] + (h / 2) * k2[1], alpha);
            k3[1] = f2(x[i] + (h / 2) * k2[0], p[i] + (h / 2) * k2[1], alpha);
            k4[0] = f1(x[i] + h * k3[0], p[i] + h * k3[1], alpha);
            k4[1] = f2(x[i] + h * k3[0], p[i] + h * k3[1], alpha);
            E[0] = (k1[0] - 4 * k2[0] + 2 * k3[0] + k4[0]) / 6;
            E[1] = (k1[1] - 4 * k2[1] + 2 * k3[1] + k4[1]) / 6;
            loc_err = sqrt(E[0] * E[0] + E[1] * E[1]);
        } while (fabs(loc_err) > EPS);
        printf("%le\n", h);
        t[i + 1] = t[i] + h;
        x[i + 1] = x[i] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        p[i + 1] = p[i] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    }
    return p[n - 1];
}

int main()
{
    int n = 250;
    double h, alpha, p0, s;
    double *t, *x, *p;
    double error;
    Point *head = NULL;

    printf("Enter alpha: ");
    scanf("%lf", &alpha);
    printf("Enter s: ");
    scanf("%lf", &s);

    t = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));
    p = (double *)malloc(n * sizeof(double));

    /*p0 = find_p0(t, x, p, n, X0, alpha);
    printf("IIIIIIIII %le\n", p0);*/
    RungeKutta(t, x, p, n, X0, 1, alpha);


    FILE *fp;
    fp = fopen("test.txt", "w");
    for(int i = 0; i < n; i++)
    {
        if (error < fabs(x[i] - solution(t[i])))
            error = fabs(x[i] - solution(t[i]));
        fprintf(fp, "%lf %lf\n", t[i], p[i]);
    }
    fclose(fp);

    printf("Error = %le\n", error);

    free(t);
    free(x);
    free(p);
    return 0;
}