#include <stdio.h>
#include <stdlib.h>
#include <cmath>

//#include <GL/gl.h>
//#include "glfw3.h"

#include <iostream>
#include <graphics.h>
#include <conio.h>

#include <gl/glut.h>

using namespace std;

//////////////
/*
    epsilon <= max([a,b])|f^(4)(x)| * (b-a)^5 / 180*n^4
    =>
    epsilon = 
*/
/////////////
const double epsilon = 1 / 1000.0;

double f(double x) {
    return x / sqrt(1 - x * x);
}

double f1(double x) {
    return pow(2, x) + 1;
}
double f2(double x) {
    return pow(x, 5);
}
double f3(double x) {
    return (1 - x) / 3;
}

double MinusF(double f(double), double g(double), double x) {
    return f(x) - g(x);
}

double Derivative(double f(double), double g(double), double x) {
    double m = MinusF(f, g, x + epsilon);
    return m / epsilon;
}

double Xopd(double f(double), double g(double), double a, double b) {
    double m1 = MinusF(f, g, a);
    double m2 = MinusF(f, g, b);
    return (b * m1 - a * m2) / (m1 - m2);
}

double Kaca(double f(double), double g(double), double x) {
    return x - MinusF(f, g, x) / Derivative(f, g, x);
}

double root(double f(double), double g(double), double a, double b) {
    while (fabs(b - a) > epsilon) {
        a = Xopd(f, g, a, b);
        b = Kaca(f, g, b);
    }
    return a;
}

double integral(double f(double), double a, double b) {
    double I = 0;
    double tmp;
    double h = (b - a) * epsilon;
    for (int i = 0; i < 1 / epsilon; i++) {
        tmp = f(a + (i + 0.5) * h);
        I += tmp;
    }
    return h * I;
}
/*
double Simpson(double f(double), double a, double b, int n) {
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; i += 2)
        s += 4 * f(a + i * h);
    for (int i = 2; i < n - 1; i += 2)
        s += 2 * f(a + i * h);
    return s * h / 3.0;
}

double Simpson2n(double f(double), double a, double b, int n) {
    const double h = (b - a) / n;

    double s = f(a) + f(b);
    for (int i = 1; i < n; i += 2)
        s += 4 * f(a + i * h);
    for (int i = 2; i < n - 1; i += 2)
        s += 2 * f(a + i * h);
    return s * h / 3.0;
}
*/
double simpson(double f(double), double a, double b) {
    double c = (a + b) / 2.0;
    return (f(a) + f(b) + 4.0 * f(c)) * (b - a) / 6.0;
}
double ars(double f(double), double a, double b, double eps) {
    double c = (a + b) / 2.0;
    double mid = simpson(f, a, b), l = simpson(f, a, c), r = simpson(f, c, b);
    if (fabs(l + r - mid) <= epsilon)
        return l + r;
    return ars(f, a, c, eps / 2.0) + ars(f, c, b, eps / 2.0);
}

void myDisplay() {
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    glVertex3f(-1, 0, 0);
    glVertex3f(1, 0, 0);
    glVertex3f(0, -1, 0);
    glVertex3f(0, 1, 0);
    glEnd();

    glBegin(GL_LINE_STRIP);

    glColor3f(1, 0, 0);
    for (double x = -10.00; x <= 10.00; x += 0.001)
        glVertex2f((float)x / 5, (float)f1(x) / 10);

    glEnd();

    glBegin(GL_LINE_STRIP);

    glColor3f(0, 1, 0);
    for (double x = -10.00; x <= 10.00; x += 0.001)
        glVertex2f((float)x / 5, (float)f2(x) / 10);

    glEnd();

    glBegin(GL_LINE_STRIP);

    glColor3f(0, 0, 1);
    for (double x = -10.00; x <= 10.00; x += 0.001)
        glVertex2f((float)x / 5, (float)f3(x) / 10);

    glEnd();
    
    glBegin(GL_LINE_LOOP);
    glColor3f(1, 1, 0);
    for (double x = -2.00; x < 10.00; x += 0.01)
    {
        for (double y = -10.00; y <= 10.00; y += 0.01)
        {
            if (y<f1(x) / 10 && y>f2(x) / 10 && y > f3(x) / 10) glVertex2f((float)x / 5, (float)y);
        }
    }
    glEnd();
    
    glFlush();
}
/*
void displayBackground() {
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    glVertex3f(-10, 0, 0);
    glVertex3f(10, 0, 0);
    glVertex3f(0, -10, 0);
    glVertex3f(0, 10, 0);
    glEnd();
}
void displayFunction() {
    glBegin(GL_LINE_STRIP);
    glColor3d(1, 0, 0);
    for (double x = -10; x <= 10; x += 0.01) {
        glVertex2d(x, f1(x));
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3d(0, 1, 0);
    for (double x = -10; x <= 10; x += 0.01) {
        glVertex2d(x, f2(x));
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3d(0, 0, 1);
    for (double x = -10; x < 0; x += 0.01) {
        glVertex2d(x, f3(x));
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3d(0, 0, 1);
    for (double x = 0.01; x <= 10; x += 0.01) {
        glVertex2d(x, f3(x));
    }
    glEnd();
}
void displayArea(double a, double b, double c) {
    glBegin(GL_LINE_STRIP);
    glColor3d(1, 1, 0);
    for (double x = a; x <= b; x += 0.01) {
        for (double y = f3(x); y <= f1(x); y += 0.01) {
            glVertex2d(x, y);
        }
    }
    for (double x = b; x <= c; x += 0.01) {
        for (double y = f2(x); y <= f1(x); y += 0.01) {
            glVertex2d(x, y);
        }
    }
    glEnd();
}
*/

int main(int argc, char* argv[]) {
    double PI = 3.1415926;
    double ans = ars(f, 0, 0.9999999999, epsilon);
    cout << ans << endl;
    system("pause");

    double  p1 = root(f1, f2, 0, 2),
        p2 = root(f1, f3, -3, -2),
        p3 = root(f2, f3, 0, 1);
    /*
    initgraph(800, 700);
    setbkcolor(WHITE);
    setlinecolor(BLACK);
    cleardevice();
    setorigin(400, 350);
    line(-400, 00, 400, 00);
    line(0, 350, 0, -350);
    int s = 80;
    for (int i = -300; i <= 300; i++) {
        line(s * i, 0, s * i, -15);
        if (i % 5 == 0) line(s * i, 0, s * i, -15);
        line(0, s * i, 10, s * i);
        if (i % 5 == 0) line(0, s * i, 15, s * i);
    }
    double x, y;
    for (x = -100; x <= 100; x += epsilon) {
        y = pow(2, x) + 1;
        putpixel(s * x, -s * y, GREEN);
        y = pow(x, 5);
        putpixel(s * x, -s * y, BLUE);
        y = (1 - x) / 3;
        putpixel(s * x, -s * y, RED);
    }*/
    system("pause");

    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutCreateWindow("OpenGL graphs");

    glutDisplayFunc(&myDisplay);
    glutMainLoop();
   
    system("pause");
    
    /*
    int epsilon1 = GetEpsilon(p2, p1),
        epsilon2 = GetEpsilon(),
        epsilon3 = GetEpsilon();
    */
    double in1 = ars(f1, p2, p1, 1 / epsilon),
        in2 = ars(f2, p3, p1, 1 / epsilon),
        in3 = ars(f3, p2, p3, 1 / epsilon);


    cout << "Intergal: " << in1 - in2 - in3 << endl;
    /*
    GLFWwindow* window;
    if (!glfwInit()) {
        return -1;
    }
    window = glfwCreateWindow(1000, 1000, "Graphic", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glOrtho(-10, 10, -10, 10, 0, 100);
    while (!glfwWindowShouldClose(window)) {
        displayBackground();
        displayFunction();
        displayArea(p1, p2, p3);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
    */
    return 0;
}
