#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define EPS 1e-8
#define X0 1
#define P1 0
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define NEWTON_METHOD_START_VALUE -0.57670298529662

typedef struct Point
{
    double t;
    double x;
    double p;
    struct Point *next;
} Point;

void pushBack(Point **last, double t_, double x_, double p_);
void clearList(Point **head);
void printList(const Point *head);
void fprintList(const Point *head, FILE *fp);
double f1(double x, double p, double alpha);
double f2(double x, double p, double alpha);
double find_p0(double x0, double alpha);
double RungeKutta(Point** head, double x0, double p0, double alpha);
double Integral(Point* head);

void pushBack(Point **last, double t_, double x_, double p_) {
    Point *pt = (Point*)malloc(sizeof(Point));
    pt->t = t_;
    pt->x = x_;
    pt->p = p_;
    pt->next = NULL;
    (*last)->next = pt;
    (*last) = pt;
}
void clearList(Point **head) {
    Point* prev = NULL;
    while ((*head)->next) {
        prev = (*head);
        (*head) = (*head)->next;
        free(prev);
    }
    free(*head);
}
void printList(const Point *head) {
    while (head) {
        printf("%le %le %le\n", head->t, head->x, head->p);
        head = head->next;
    }
}
void fprintList(const Point *head, FILE *fp) {
    while (head) {
        fprintf(fp, "%lf %lf %lf\n", head->t, head->x, head->p);
        head = head->next;
    }
}

double f1(double x, double p, double alpha)
{
    double t =  M_PI / (4 + alpha * x * x);
    return p + sin(t);
}
double f2(double x, double p, double alpha)
{
    double t =  4 + alpha * x * x;
    return (2 * alpha * M_PI * x * p * cos(M_PI / t)) / (t * t) + 0.5 * x;
}
double find_p0(double x0, double alpha)
{
    double f_left, f_right;
    double df = 1;
    double p0 = NEWTON_METHOD_START_VALUE;
    double p1 = 0;
    Point *head;

    do
    {
        p0 = p0 - p1 / df;
        head = (Point*)malloc(sizeof(Point));
        p1 = RungeKutta(&head, X0, p0, alpha);
        clearList(&head);
        head = (Point*)malloc(sizeof(Point));
        f_left = RungeKutta(&head, X0, p0 - EPS, alpha);
        clearList(&head);
        head = (Point*)malloc(sizeof(Point));
        f_right = RungeKutta(&head, X0, p0 + EPS, alpha);
        clearList(&head);
        df = (f_right - f_left) / (2 * EPS);
        printf("%le\n", p1);
    } while (fabs(p1) > EPS);
    printf("p(0) = %.11f\n", p0);
    return p0;
}
double RungeKutta(Point** head, double x0, double p0, double alpha)
{
    int a;
    double h, loc_err, xh;
    double next_t, next_x, next_p;
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    double k5[2];
    double k6[2];
    Point* last;

    Point *pt = (Point*)malloc(sizeof(Point));
    pt->t = 0;
    pt->x = x0;
    pt->p = p0;
    pt->next = NULL;

    (*head) = pt;
    last = pt;

    h = 0.5;
    while (1 - last->t > EPS)
    {
        if (h > 1 - last->t)
        {
            h = 1 - last->t;
        }
        h = h * 2;
        next_x = DBL_MAX;
        do
        {
            xh = next_x;
            k1[0] = h * f1(last->x, last->p, alpha);
            k1[1] = h * f2(last->x, last->p, alpha);
            k2[0] = h * f1(last->x + (1. / 2) * k1[0], last->p + (1. / 2) * k1[1], alpha);
            k2[1] = h * f2(last->x + (1. / 2) * k1[0], last->p + (1. / 2) * k1[1], alpha); 
            k3[0] = h * f1(last->x + (1. / 4) * (k1[0] + k2[0]), last->p + (1. / 4) * (k1[1] + k2[1]), alpha);
            k3[1] = h * f2(last->x + (1. / 4) * (k1[0] + k2[0]), last->p + (1. / 4) * (k1[1] + k2[1]), alpha);
            k4[0] = h * f1(last->x + (-k2[0] + 2 * k3[0]), last->p + (-k2[1] + 2 * k3[1]), alpha);
            k4[1] = h * f2(last->x + (-k2[0] + 2 * k3[0]), last->p + (-k2[1] + 2 * k3[1]), alpha);
            k5[0] = h * f1(last->x + (1. / 27) * (7 * k1[0] + 10 * k2[0] + k4[0]), last->p + (1. / 27) * (7 * k1[1] + 10 * k2[1] + k4[1]), alpha);
            k5[1] = h * f2(last->x + (1. / 27) * (7 * k1[0] + 10 * k2[0] + k4[0]), last->p + (1. / 27) * (7 * k1[1] + 10 * k2[1] + k4[1]), alpha);
            k6[0] = h * f1(last->x + (1. / 625) * (28 * k1[0] - 125 * k2[0] + 546 * k3[0] + 54 *k4[0] - 378 * k5[0]), last->p + (1. / 625) * (28 * k1[1] - 125 * k2[1] + 546 * k3[1] + 54 * k4[1] - 378 * k5[1]), alpha);
            k6[1] = h * f2(last->x + (1. / 625) * (28 * k1[0] - 125 * k2[0] + 546 * k3[0] + 54 *k4[0] - 378 * k5[0]), last->p + (1. / 625) * (28 * k1[1] - 125 * k2[1] + 546 * k3[1] + 54 * k4[1] - 378 * k5[1]), alpha);
            //printf("%le ^^^^ %le\n", last->t, h);
            next_t = last->t + h;
            next_x = last->x + k1[0] / 24 + 5 * k4[0] / 48 + 27 * k5[0] / 56 + 125 * k6[0] / 336;
            next_p = last->p + k1[1] / 24 + 5 * k4[1] / 48 + 27 * k5[1] / 56 + 125 * k6[1] / 336;
            if (fabs(xh - next_x) / 31 > EPS)
            {
                h = h / 2;
            }
            else 
            {
                break;
            }
            /*printf("   %lf ||| %lf ||| %le ||| %le\n", xh, next_x, h, last->t);
            scanf("%d", &a);*/
        } while (1);

        pushBack(&last, next_t, next_x, next_p);
    }
    printf("x(1) =  %.11f\n", last->x);
    return last->p;
}
double Integral(Point* head)
{
    double f_a, f_b;
    double res = 0;
    while (head->next)
    {
        f_a = head->x * head->x + 2 * head->p * head->p;
        f_b = head->next->x * head->next->x + 2 * head->next->p * head->next->p;
        res = res + (f_a + f_b) * (head->next->t - head->t) / 2;
        head = head->next;
    }
    return res;
}

int main()
{
    double h, alpha, p0, B;
    double error;
    Point *head = (Point *)malloc(sizeof(head));

    printf("Enter alpha: ");
    scanf("%lf", &alpha);

    p0 = find_p0(X0, alpha);
    RungeKutta(&head, X0, p0, alpha);

    printf("Integral = %.11f\n", Integral(head));

    //printList(head);

    /*FILE *file;
    file = fopen("test.txt", "w");
    fprintList(head, file);
    fclose(file);*/

    clearList(&head);
    return 0;
}