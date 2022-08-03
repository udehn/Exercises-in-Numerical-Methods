#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <cmath>
#include <time.h>
#include <complex>

using namespace std;
template <typename T>
class Matrix
{
private:
    vector<vector<T>> Matr;
    vector<vector<T>> L;
    vector<vector<T>> U;
    vector<T> b;
    vector<T> y;
    vector<T> x;
    int row, col;
public:
    Matrix(int n);
    void Init();
    void input_f(int n);
    void Input_f();
    void get_LU();
    void calculate_x_y();
    void Print();
    void Print_x();
};

template <typename T>
Matrix<T>::Matrix(int n) {
    row = n;
    col = n;
}

template <typename T>
void Matrix<T>::Init() {
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
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
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
}

template <typename T>
void Matrix<T>::Print_x() {
    for (int i = 0; i < row; i++)
        cout << x[i] << " ";
    cout << endl;
}

template <typename T>
void Matrix<T>::get_LU()
{
    row = Matr.size();
    col = row;
    T tmp = 0;
    L.resize(row);
    U.resize(row);
    for (auto& i : L)
        i.resize(row);
    for (auto& i : U)
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
void Matrix<T>::input_f(int n) {
    for (int i = 0; i < n; i++)
        b.push_back(1);
}
/*
template <typename T>
void Matrix<T>::Input_f() {
    double temp = 0;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            temp += matr[i][j];
        }
        b.push_back(temp);
        temp = 0;
    }
}*/

template <typename T>
void Matrix<T>::calculate_x_y()
{
    int n = row;
    x.resize(n);
    y.resize(n);
    y[0] = b[0];
    for (int i = 1; i < n; i++)
    {
        y[i] = b[i];
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

int main()
{
    int n = 10;
    Matrix<int> m(n);
    m.Init();
    m.input_f(n);
    m.get_LU();
    m.calculate_x_y();
    //m.Print();
    m.Print_x();
    return 0;
}
