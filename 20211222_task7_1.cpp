#include <stdio.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

#include "glfw3.h"
#include <GL/gl.h>
#include <iomanip>

#define PI 3.1415926

using namespace std;

const int n = 1000;
double x[100000];
double u[100000];
double z[100000];
double S(double p)
{
    vector<double> step(n+1), s(n+1), F(n+1);
    vector<vector<double>> m(n+1);
    for (int i = 0; i <= n; i++)
    {
        m[i].resize(n+1);
    }
    for (int i = n; i > 0; i--)
    {
        F[i] = (u[i] - u[i - 1]) / (x[i] - x[i - 1]);
        step[i - 1] = x[i] - x[i - 1];
    }
    for (int i = 1; i <= n - 1; i++)
    {
        m[i][i] = 2 * (step[i - 1] + step[i]);
        if (i != 1)
        {
            m[i][i - 1] = step[i - 1];
            m[i - 1][i] = step[i - 1];
        }
        m[i][n - 1] = 6 * (F[i + 1] - F[i]);
    }

    for (int i = 1; i < n - 1; i++)
    {
        double temp = (m[i + 1][i] / m[i][i]);
        for (int j = 1; j <= n; j++)
            m[i + 1][j] -= temp * m[i][j];
    }

    for (int i = n - 1; i > 0; i--)
    {
        double sum = 0;
        for (int j = i; j <= n - 1; j++)
            sum += m[i][j] * s[j];
        s[i] = (m[i][n - 1] - sum) / m[i][i];
    }
    double res = 0;
    for (int i = 0; i <= n - 1; i++)
    {
        if (x[i] <= p && p <= x[i + 1])
        {
            double a1 = (s[i + 1] - s[i]) / (6 * step[i]);
            double b1 = s[i] / 2;
            double c1 = (u[i + 1] - u[i]) / step[i] - (2 * step[i] * s[i] + s[i + 1] * step[i]) / 6;
            double d1 = u[i];
            res = a1 * pow((p - x[i]), 3) + b1 * pow((p - x[i]), 2) + c1 * (p - x[i]) + d1;
        }
    }
    return res;
}

void displayBackground()
{
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    glVertex3f(-30, 0, 0);
    glVertex3f(30, 0, 0);
    glVertex3f(0, -30, 0);
    glVertex3f(0, 30, 0);
    glEnd();
}

void displaySolution()
{
    glBegin(GL_POINTS);
    glColor3d(1, 0, 0);
    for (int i = 0; i <= n; i++)
    {
        glVertex2d(x[i], u[i]);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3d(0, 1, 0);
    //for (int i = 0;i < n;i++) {        glVertex2d(x[i],u[i]);    }
    for (double x = 0; x <= 4 * PI; x += 0.01)
    {
        glVertex2d(x, S(x));
    }
    glEnd();
}
/*
double f(double tn, double un)
{
    return sqrt(2 * cos(un) - 2 * cos(1) );
}
*/
double f(double tn,double zn) {
    return zn;
}
double g(double tn,double un) {
    return -sin(un);
}

void RungeKutta(double step, double l, double r, double *x, double *u, int tol)
{
    double kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4;
    for (int i = 0;i < n;i++) {
        kf1 = f(x[i], z[i]);
        kg1 = g(x[i], u[i]);
        kf2 = f(x[i] + step / 2.0, z[i] + step * kg1 / 2.0);
        kg2 = g(x[i] + step / 2.0, u[i] + step * kf1 / 2.0);
        kf3 = f(x[i] + step / 2.0, z[i] + step * kg2 / 2.0);
        kg3 = g(x[i] + step / 2.0, u[i] + step * kf2 / 2.0);
        kf4 = f(x[i] + step, z[i] + step * kg3);
        kg4 = g(x[i] + step, u[i] + step * kf3);
        
        u[i + 1] = u[i] + step * (kf1 + 2 * kf2 + 2 * kf3 + kf4) / 6.0;
        z[i + 1] = z[i] + step * (kg1 + 2 * kg2 + 2 * kg3 + kg4) / 6.0;
    }
}

int main()
{
    double l, r;
    /*
    memset(x, 0, sizeof(x));
    memset(u, 0, sizeof(u));
    */
    l = 0;
    r = 4 * PI;
    x[0] = l;
    u[0] = 1;
    z[0] = 0;
    double step = (r - l) / n;
    
    for (int i = 1; i <= n; i++)
    {
        x[i] = x[i-1] + step;
    }
    RungeKutta(step, l, r, x, u, n);

    GLFWwindow *window;
    if (!glfwInit())
    {
        return -1;
    }
    window = glfwCreateWindow(1000, 1000, "Grapstepic", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glOrtho(-20, 20, -20, 20, 0, 100);
    while (!glfwWindowShouldClose(window))
    {
        displayBackground();
        displaySolution();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();

    return 0;
}