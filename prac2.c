#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-9
#define K 32
#define X0 1
#define P1 0
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

void pushBack(Point **last, double t_, double x_, double p_);
void clearList(Point **head);
void printList(const Point *head);
void fprintList(const Point *head, FILE *fp);
double f1(double x, double p, double alpha);
double f2(double x, double p, double alpha);
double find_p0(double x0, double alpha);
double RungeKutta(Point** head, double x0, double p0, double alpha);

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
    double p0 = 0;
    double p1 = 0;
    Point *head;
    int a;

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
    printf("%le\n", p0);
    return p0;
}
double RungeKutta(Point** head, double x0, double p0, double alpha)
{
    int a;
    double h, loc_err;
    double E1, E2;
    double next_t, next_x, next_p;
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    Point* last;

    Point *pt = (Point*)malloc(sizeof(Point));
    pt->t = 0;
    pt->x = x0;
    pt->p = p0;
    pt->next = NULL;

    (*head) = pt;
    last = pt;

    while (1 - last->t > EPS)
    {
        h = 1;
        do
        {
            k1[0] = f1(last->x, last->p, alpha);
            k1[1] = f2(last->x, last->p, alpha);
            k2[0] = f1(last->x + (h / 2) * k1[0], last->p + (h / 2) * k1[1], alpha);
            k2[1] = f2(last->x + (h / 2) * k1[0], last->p + (h / 2) * k1[1], alpha); 
            k3[0] = f1(last->x + (h / 2) * k2[0], last->p + (h / 2) * k2[1], alpha);
            k3[1] = f2(last->x + (h / 2) * k2[0], last->p + (h / 2) * k2[1], alpha);
            k4[0] = f1(last->x + h * k3[0], last->p + h * k3[1], alpha);
            k4[1] = f2(last->x + h * k3[0], last->p + h * k3[1], alpha);
            E1 = (k1[0] - 4 * k2[0] + 2 * k3[0] + k4[0]) / 6;
            E2 = (k1[1] - 4 * k2[1] + 2 * k3[1] + k4[1]) / 6;
            loc_err = sqrt(E1 * E1 + E2 * E2);
            if (loc_err > EPS)
            {
                h = h / 2;
            }
            if (loc_err < EPS / K)
            {
                h = h * 2;
            }
        } while (fabs(loc_err) > EPS);
        next_t = last->t + h;
        next_x = last->x + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        next_p = last->p + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
        pushBack(&last, next_t, next_x, next_p);
    }
    return last->p;
}

int main()
{
    double h, alpha, p0;
    double error;

    Point *head = (Point *)malloc(sizeof(head));

    printf("Enter alpha: ");
    scanf("%lf", &alpha);

    p0 = find_p0(X0, alpha);
    RungeKutta(&head, X0, p0, alpha);

    //printList(head);

    FILE *file;
    file = fopen("test.txt", "w");
    fprintList(head, file);
    fclose(file);

    clearList(&head);
    return 0;
}