#include <iostream>
#include <cmath>
#include <gl/glut.h>
#include <vector>

#include <graphics.h>
#include <conio.h>

using namespace std;

const int num = 25;
double u[num][num];

double yy[10000];
double xx[10000];


double S(double p) {
    int n = 100;
    vector<double> h(n), s(n), F(n);
    vector<vector<double> > m(n);
    for (int i = 0; i < n; i++) {
        m[i].resize(n);
    }
    for (int i = n - 1; i > 0; i--) {
        F[i] = (yy[i] - yy[i - 1]) / (xx[i] - xx[i - 1]);
        h[i - 1] = xx[i] - xx[i - 1];
    }
    for (int i = 1; i < n - 1; i++) {
        m[i][i] = 2 * (h[i - 1] + h[i]);
        if (i != 1) {
            m[i][i - 1] = h[i - 1];
            m[i - 1][i] = h[i - 1];
        }
        m[i][n - 1] = 6 * (F[i + 1] - F[i]);
    }

    for (int i = 1; i < n - 2; i++) {
        double temp = (m[i + 1][i] / m[i][i]);
        for (int j = 1; j <= n - 1; j++)
            m[i + 1][j] -= temp * m[i][j];
    }

    for (int i = n - 2; i > 0; i--) {
        double sum = 0;
        for (int j = i; j <= n - 2; j++)
            sum += m[i][j] * s[j];
        s[i] = (m[i][n - 1] - sum) / m[i][i];
    }
    double res = 0;
    for (int i = 0; i < n - 1; i++) {
        if (xx[i] <= p && p <= xx[i + 1]) {
            double a1 = (s[i + 1] - s[i]) / (6 * h[i]);
            double b1 = s[i] / 2;
            double c1 = (yy[i + 1] - yy[i]) / h[i] - (2 * h[i] * s[i] + s[i + 1] * h[i]) / 6;
            double d1 = yy[i];
            res = a1 * pow((p - xx[i]), 3) + b1 * pow((p - xx[i]), 2) + c1 * (p - xx[i]) + d1;
        }
    }
    return res;
}

/*
double getij(int indexI, int indexJ, double step)
{
    cout << indexI << ", " << indexJ << endl;
    //if (indexI == 0 || indexI == num - 1 || indexJ == 0 || indexJ == num - 1)
      //  return u[indexI][indexJ];
    if (u[indexI][indexJ] != -1) {
        cout << "return " << indexI << ", " << indexJ << endl;
        return u[indexI][indexJ];
    }

    double aij = 1 / (4 - step / 10) * (step * step
        + (1 - step / 10) * getij(indexI + 1, indexJ, step)
        + getij(indexI, indexJ + 1, step)
        + getij(indexI - 1, indexJ, step)
        + getij(indexI, indexJ - 1, step));
    u[indexI][indexJ] = aij;
    return aij;
}
*/
int Solve(double omega, double b) {
    double epsilon = 1e-6;
    double num_of_iteration = 0;
    double step = (1. - 0) / (double)num;
    double d = 2 / (step * step) + 2 / (step * step);
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            u[i][j] = 0;
        }
    }
    /*  
    for (int i = 0; i < num; i++) {  
        u[0][i] = 0;
        u[num - 1][i] = 0;
        u[i][0] = 0;
        u[i][num - 1] = 0;
    }
    */
    
    /*
    for (int i = 1; i < num - 1; i++)
    {
        for (int j = 1; j < num - 1; j++)
        {
            u[i][j] = getij(i, j, step);
            cout << u[i][j] << " ";
        }
        cout << endl;
    }
    */
    double flag = 10;
    while (sqrt(flag) > epsilon) {
        flag = 0;
        for (int j = 1; j < num - 1; j++)
        {
            for (int i = 1; i < num - 1; i++)
            {
                double temp = -(u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]) / (step * step) - (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]) / (step * step) + b * (u[i + 1][j] - u[i - 1][j]) / (2 * step) - 1;
                flag += temp * temp;
                u[i][j] -= omega * temp / d;
            }
            flag *= step * step;
        }
        num_of_iteration++;
    }
    /*
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            cout.width(10);
            cout.precision(3);
            cout << u[i][j] << " ";
        }
        cout << endl;
    }*/

    cout << "omega: " << omega << endl << "num_of_iteration: " << num_of_iteration << endl << endl;
    return num_of_iteration;
}

double Gold(double left, double right, double b)
{
    double epsilon = 1e-7, tau = (sqrt(5) - 1) / 2;
    double a_gss = left, b_gss = right, f_lambda = 0, f_mu = 0;
    double lambda = a_gss + (1.0 - tau) * (b_gss - a_gss), mu = a_gss + tau * (b_gss - a_gss);
    f_lambda = Solve(lambda, b);
    f_mu = Solve(mu, b);
    while (fabs(mu - lambda) > epsilon)
    {
        if (f_lambda < f_mu)
        {
            b_gss = mu;
            mu = lambda;
            f_mu = f_lambda;
            lambda = a_gss + (1 - tau) * (b_gss - a_gss);
            f_lambda = Solve(lambda, b);
        }
        else
        {
            a_gss = lambda;
            lambda = mu;
            f_lambda = f_mu;
            mu = a_gss + tau * (b_gss - a_gss);
            f_mu = Solve(mu, b);
        }
    }
    return 0.5 * (lambda + mu);
}

void Display()
{
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    glVertex3f(1, 0, 0);
    glVertex3f(2, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1000, 0);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glColor3f(1, 0, 0);
    for (double i = 1; i <= 2; i += 1e-2)
        glVertex2f((float)i, (float)Solve(i, 0));
    glEnd();


    glBegin(GL_LINE_STRIP);
    glColor3f(0, 1, 0);
    for (double i = 1; i <= 2; i += 1e-2)
        glVertex2f((float)i, (float)Solve(i, 10));
    /*
    double step = 1e-2;
    double Index = 1;
    for (int i = 0; i <= 100; i++) {
        yy[i] = Solve(Index);
        xx[i] += step;
        Index += step;
    }
    
    for (double i = 1; i <= 2; i += 1e-2)
        glVertex2f(i, S(i));
    */
    glEnd();
    cout << "testtttttttttttttttttttt" << endl;

    glFlush();

}

int main(int argc, char* argv[])
{
    //cout << Gold(1, 2) << endl;    
    cout << Gold(1, 2, 0) << endl;
    cout << Gold(1, 2, 10) << endl;

    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(700, 700);
    glutCreateWindow("graph");

    glClearColor(1.0, 1.0, 1.0, 0.0);
    glMatrixMode(GL_PROJECTION | GL_MODELVIEW);
    gluOrtho2D(1, 2, 0, 1000);
    
    glutDisplayFunc(&Display);
    glutMainLoop();

}