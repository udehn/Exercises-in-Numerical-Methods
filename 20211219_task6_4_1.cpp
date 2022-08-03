#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>

#include <cmath>
#include <time.h>

//#include <GL/glut.h>

#include "glfw3.h"
#include <GL/gl.h>
#include <iomanip>

#include <complex>

#define PI 3.14159265358979

double epsilon = 0.01;

const int n = 201;
double yy[n];
double xx[n];

using namespace std;

double S(double p)
{
    vector<double> h(n), s(n), F(n);
    vector<vector<double>> m(n);
    for (int i = 0; i < n; i++)
    {
        m[i].resize(n);
    }
    for (int i = n - 1; i > 0; i--)
    {
        F[i] = (yy[i] - yy[i - 1]) / (xx[i] - xx[i - 1]);
        h[i - 1] = xx[i] - xx[i - 1];
    }
    for (int i = 1; i < n - 1; i++)
    {
        m[i][i] = 2 * (h[i - 1] + h[i]);
        if (i != 1)
        {
            m[i][i - 1] = h[i - 1];
            m[i - 1][i] = h[i - 1];
        }
        m[i][n - 1] = 6 * (F[i + 1] - F[i]);
    }

    for (int i = 1; i < n - 2; i++)
    {
        double temp = (m[i + 1][i] / m[i][i]);
        for (int j = 1; j <= n - 1; j++)
            m[i + 1][j] -= temp * m[i][j];
    }

    for (int i = n - 2; i > 0; i--)
    {
        double sum = 0;
        for (int j = i; j <= n - 2; j++)
            sum += m[i][j] * s[j];
        s[i] = (m[i][n - 1] - sum) / m[i][i];
    }
    double res = 0;
    for (int i = 0; i < n - 1; i++)
    {
        if (xx[i] <= p && p <= xx[i + 1])
        {
            double a1 = (s[i + 1] - s[i]) / (6 * h[i]);
            double b1 = s[i] / 2;
            double c1 = (yy[i + 1] - yy[i]) / h[i] - (2 * h[i] * s[i] + s[i + 1] * h[i]) / 6;
            double d1 = yy[i];
            res = a1 * pow((p - xx[i]), 3) + b1 * pow((p - xx[i]), 2) + c1 * (p - xx[i]) + d1;
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
    glVertex3f(-10, 0, 0);
    glVertex3f(10, 0, 0);
    glVertex3f(0, -10, 0);
    glVertex3f(0, 10, 0);
    glEnd();
}
void displaySolution()
{
    glBegin(GL_POINTS);
    glColor3d(1, 0, 0);
    for (int i = 0; i < n; i++)
    {
        glVertex2d(xx[i], yy[i]);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3d(0, 1, 0);
    //for (int i = 0;i < n;i++) {        glVertex2d(x[i],y[i]);    }
    for (double x = -1; x < 1; x += 0.01)
    {
        glVertex2d(x, S(x));
    }
    glEnd();
}

double K(double x, double s)
{
    double miu = 0.05;
    //double miu = 0.1;
    //double miu = 0.2;
    return miu / (pow(miu, 2) + pow(x - s, 2));
}

double G(double x, double s, double a, double b, int num)
{
    double ans = 0;
    double step = (b - a) / ((double)num - 1);
    double t = a;
    for (int i = 0; i < num; i++)
    {
        ans += K(t, x) * K(t, s) * step;
        t += step;
    }
    return ans;
}

double f_wave(double x)
{
    return cos(PI * x);
}

template <typename T>
class Matrix
{
private:
    vector<vector<T>> Matr;
    vector<vector<T>> L;
    vector<vector<T>> U;
    vector<T> f;
    vector<T> y;
    vector<T> x;
    int row, col;
    //double alpha1;
public:
    Matrix(int n);
    void Init();
    void input_f(int n);
    void Init_Fredholm(double a, double b, int num);

    void Init_Tikhonov(double a, double b, int num, double alpha);
    void input_f_Tikhonov(double a, double b, int num);
    void Input_xx_yy(double a, double b, int num);

    void Input_f();

    void get_LU();
    void calculate_x_y();
    void Print();
    void Print_x();

    void Gold(double a, double b, int num, double left, double right);
    double Solve(double a, double b, int num, double alpha);

    void Clean();
};

template <typename T>
Matrix<T>::Matrix(int n)
{
    row = n;
    col = n;
}

template <typename T>
void Matrix<T>::Init()
{
    /*       vector<T> mat_vet;
       for (double i = 0; i < row; i++) {
           for (double j = 0; j < col; j++) {
               mat_vet.push_back(1 / (i + j + 1));
           }
           Matr.push_back(mat_vet);
           mat_vet.clear();
       }
   */
    vector<T> mat_vet;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (j == col - 1)
                mat_vet.push_back(1);
            else if (i == j)
                mat_vet.push_back(1);
            else if (i < j)
                mat_vet.push_back(-1);
            else if (i > j && j != col - 1)
                mat_vet.push_back(0);
        }
        Matr.push_back(mat_vet);
        mat_vet.clear();
    }
}

template <typename T>
void Matrix<T>::Print()
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            cout << Matr[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    /*
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            cout << L[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            cout << U[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    */
}

template <typename T>
void Matrix<T>::Print_x()
{
    /*
    for (int i = 0; i < row; i++)
        cout << x[i] << " ";
    cout << endl;
    */
    for (int i = 0; i < row; i++)
        cout << x[i] << endl;
}

template <typename T>
void Matrix<T>::get_LU()
{
    row = Matr.size();
    col = row;
    T tmp = 0;
    L.resize(row);
    U.resize(row);
    for (auto &i : L)
        i.resize(row);
    for (auto &i : U)
        i.resize(row);

    for (int i = 0; i < row; i++)
        for (int j = 0; j <= i; j++)
            if (i == j)
                L[i][j] = 1;
            else
                U[i][j] = 0;

    for (int i = 0; i < row; i++)
    {
        U[0][i] = Matr[0][i];
        L[i][0] = Matr[i][0] / U[0][0];
    }
    for (int k = 1; k < row; k++)
    {
        for (int j = k; j < row; j++)
        {
            tmp = 0;
            for (int m = 0; m < k; m++)
                tmp += L[k][m] * U[m][j];
            U[k][j] = Matr[k][j] - tmp;
        }
        for (int i = k + 1; i < row; i++)
        {
            tmp = 0;
            for (int m = 0; m < k; m++)
                tmp += L[i][m] * U[m][k];
            L[i][k] = (Matr[i][k] - tmp) / U[k][k];
        }
    }
}

template <typename T>
void Matrix<T>::input_f(int n)
{
    for (int i = 0; i < n; i++)
        f.push_back(1);
}
/*
template <typename T>
void Matrix<T>::Input_f() {
    double temp = 0;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            temp += matr[i][j];
        }
        f.push_back(temp);
        temp = 0;
    }
}*/

template <typename T>
void Matrix<T>::calculate_x_y()
{
    int n = row;
    x.resize(n);
    y.resize(n);
    y[0] = f[0];
    for (int i = 1; i < n; i++)
    {
        y[i] = f[i];
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
    }
    x[n - 1] = y[n - 1] / U[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }
}

template <typename T>
void Matrix<T>::Init_Fredholm(double a, double b, int num)
{
    double step = (b - a) / ((double)num - 1);
    //cout << step << endl;
    vector<double> mat_vet;
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            /*
            if (j == col - 1)
                mat_vet.push_back(1);
            else if (i == j)
                mat_vet.push_back(1);
            else if (i < j)
                mat_vet.push_back(-1);
            else if (i > j && j != col - 1)
                mat_vet.push_back(0);
            */
            if (i == j)
                mat_vet.push_back(1 - step / (PI * (1 + pow((a + step * i) - (a + step * j), 2))));
            else
                mat_vet.push_back(-step / (PI * (1 + pow((a + step * i) - (a + step * j), 2))));
        }
        Matr.push_back(mat_vet);
        mat_vet.clear();
    }
}

template <typename T>
void Matrix<T>::Init_Tikhonov(double a, double b, int num, double alpha)
{
    double step = (b - a) / ((double)num - 1);
    vector<double> mat_vet;
    double x = a;
    double s = a;
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            if (i == j)
                mat_vet.push_back(alpha + step * G(x, s, a, b, num));
            else
                mat_vet.push_back(step * G(x, s, a, b, num));
            s += step;
            //cout << mat_vet[j] << "  ";
        }
        //cout << endl;
        Matr.push_back(mat_vet);
        mat_vet.clear();
        x += step;
        s = a;
    }
}

template <typename T>
void Matrix<T>::input_f_Tikhonov(double a, double b, int num)
{
    double temp = 0;
    double xi = a;
    double sj = a;
    double step = (b - a) / ((double)num - 1);
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            temp += step * K(sj, xi) * f_wave(sj);
            sj += step;
        }
        f.push_back(temp);
        temp = 0;
        sj = a;
        xi += step;
    }
}

template <typename T>
double Matrix<T>::Solve(double a, double b, int num, double alpha)
{
    Init_Tikhonov(a, b, num, alpha);
    input_f_Tikhonov(a, b, num);
    get_LU();
    calculate_x_y();

    vector<double> f_new;
    double temp = 0;
    double xi = a;
    double sj = a;
    double step = (b - a) / ((double)num - 1);

    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            temp += K(xi, sj) * x[j] * step;
            sj += step;
        }
        f_new.push_back(temp);
        temp = 0;
        sj = a;
        xi += step;
    }

    double ans = 0;
    xi = a;
    for (int i = 0; i < num; i++)
    {
        ans += pow(f_new[i] - f_wave(xi), 2);
        xi += step;
    }
    ans = sqrt(ans);
    /*
    cout << "alpha = " << alpha << endl;
    cout << "eps = " << ans << endl << endl;
    */
    //cout << ans << endl;
    return ans;
}

template <typename T>
void Matrix<T>::Gold(double a, double b, int num, double left, double right)
{
    input_f_Tikhonov(a, b, num);
    double epsilon = 1e-6, f_min, tau = (sqrt(5) - 1) / 2.0;
    double a_gss = left, b_gss = right, f_lambda = 0, f_mu = 0;
    double lambda = a_gss + (1.0 - tau) * (b_gss - a_gss), mu = a_gss + tau * (b_gss - a_gss);
    f_lambda = Solve(a, b, num, lambda);
    Clean();
    f_mu = Solve(a, b, num, mu);
    Clean();

    for (int i = 0; fabs(mu - lambda) > epsilon; i++)
    {
        if (f_lambda < f_mu)
        {
            b_gss = mu;
            mu = lambda;
            f_mu = f_lambda;
            lambda = a_gss + (1.0 - tau) * (b_gss - a_gss);
            f_lambda = Solve(a, b, num, lambda);
            Clean();
        }
        else
        {
            a_gss = lambda;
            lambda = mu;
            f_lambda = f_mu;
            mu = a_gss + tau * (b_gss - a_gss);
            f_mu = Solve(a, b, num, mu);
            Clean();
        }
    }
    //f_min = Solve(0.5 * (lambda + mu));  Clean();
    cout << "optimal alpha and epsilon:  " << endl;
    bool Compare = (f_lambda < f_mu);
    if (Compare)
    {
        cout << "alpha = " << lambda << endl;
        cout << "eps = " << f_lambda << endl
             << endl;
    }
    else
    {
        cout << "alpha = " << mu << endl;
        cout << "eps = " << f_mu << endl
             << endl;
    }
}

template <typename T>
void Matrix<T>::Clean()
{
    Matr.clear();
    L.clear();
    U.clear();
    y.clear();
    x.clear();
    f.clear();
}

template <typename T>
void Matrix<T>::Input_xx_yy(double a, double b, int num)
{
    double xxx = -1;
    double step = (b - a) / (num - 1);
    for (int i = 0; i < n; i++)
    {
        yy[i] = x[i];
        xx[i] = xxx;
        xxx += step;
        //cout << yy[i] << "   " << xx[i] << endl;
    }
}

int main()
{

    /////////////////2_1

    Matrix<double> m2(n);
    m2.Gold(-1, 1, n, 0, 1);
    m2.Input_xx_yy(-1, 1, n);

    /////////////////

    ///////////////// 1_4
/*    
    int n1 = 100;
    Matrix<double> m1(n1);
    m1.Init_Fredholm(-1, 1, n1);
    
    m1.input_f(n1);
    m1.get_LU();    
    m1.calculate_x_y();
    m1.Input_xx_yy(-1, 1, n1);
    //cout << "step is: " << (double)2 / n << endl;
    //m1.Print_x();
  */  
    /////////////////

    GLFWwindow *window;
    if (!glfwInit())
    {
        return -1;
    }
    window = glfwCreateWindow(1000, 1000, "Graphic", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glOrtho(-5, 5, -5, 5, 0, 100);
    while (!glfwWindowShouldClose(window))
    {
        displayBackground();
        displaySolution();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();

    /*
    int n2 = 201;
    Matrix<double> m2(n2);
    m2.Solve(-1, 1, 201, 0.140717);
    m2.Print_x();
      */
    ///////////////////

    ///////////////// 1_4
    /*
    int n1 = 100;
    Matrix<double> m1(n1);
    m1.Init_Fredholm(-1, 1, n1);
    
    m1.input_f(n1);
    m1.get_LU();
    
m1.calculate_x_y();
    
        //cout << "step is: " << (double)2 / n << endl;
    m1.Print_x();
    */
    //////////////////

    return 0;
}
