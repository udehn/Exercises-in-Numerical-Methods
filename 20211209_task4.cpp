/*================================================================
*   Copyright (C) 2020 ChenXuandong Ltd. All rights reserved.
*
*   FileName:matrix.cpp
*   Author:ChenXuandong
*   Data:2020年12月18日
*   Description：
*
================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <string.h>

using namespace std;

class matrix {
private:
	int row;
	int col;
	double* matr;
public:
	matrix();
	matrix(int n, int m);
	matrix(const matrix& m1);
	matrix(double data);
	matrix(double* data, int m);
	matrix(int n, double* data);
	matrix(char* str);
	static matrix identity(int n);
	static matrix diagonal(double* vals, int n);

	void Init_matr();
	void Input_f(matrix& m);
	void Init_x();


	int rows() const;
	int columns() const;
	void InitRnd();
	double* getmatr();
	void Set(int nX, int nY, double r);
	double Get(int nX, int nY) const;
	~matrix() { delete[] matr; };
	matrix& operator *= (double data);
	matrix& operator *= (double* v);
	matrix& operator += (matrix& m1);
	matrix& operator -= (matrix& m1);
	matrix& operator *= (matrix& m1);
	matrix& operator - ();
	matrix operator [] (int i);
	friend ostream& operator<<(ostream& os, matrix& m);
	void print();
	void Divid(double R, int nVerh);
	void Minus(double R, int nNizz, int nVerh);
	void Tuda();
	void Obratno();
	double determinant();
	void inverse();
	void transpose();

	void ToB(matrix& A, double omega);
	matrix solve();
	void Copy(matrix& m1);
	double GetError();
};

matrix matrix::operator [](int i) {
	matrix mrow(1, col);
	for (int m = 1; m <= col; m++)
		mrow.Set(1, m, Get(i, m));
	return mrow;
}

matrix matrix::identity(int n) {
	if (n <= 0)
		throw("error matrix::identity m(n);\n");
	matrix tmp(n, n);
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			if (i == j) tmp.Set(i, j, 1.0);
			else tmp.Set(i, j, 0);
	return tmp;
}

matrix matrix::diagonal(double* vals, int n) {
	if (n <= 0 || vals == nullptr)
		throw("error matrix::diagonal m(vals, n);\n");
	matrix tmp(n, n);
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			if (i == j) tmp.Set(i, j, vals[i - 1]);
			else tmp.Set(i, j, 0);
	return tmp;
}

matrix::matrix() {
	row = 1;
	col = 1;
	matr = new double[1];
	Set(1, 1, 0);
}

matrix::matrix(int n, int m) {
	if (n <= 0 || m <= 0)
		throw("error matrix m(n, m);\n");
	row = n;
	col = m;
	matr = new double[row * col];
	/*
	matr = new double[row * col];
	for (int i = 1; i <= row; i++)
		for (int j = 1; j <= col; j++)
			Set(i, j, 0);
	*/
}

matrix::matrix(const matrix& m1) {
	row = m1.rows();
	col = m1.columns();
	if (matr != nullptr)
		delete[] matr;
	matr = new double[row * col];
	//matr = m1.matr;
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, m1.Get(i, j));
}

matrix::matrix(double data) {
	row = 1;
	col = 1;
	matr = new double[1];
	Set(1, 1, data);
}

matrix::matrix(double* data, int m) {
	if (m <= 0 || data == nullptr)
		throw("error: matrix m(double*, int);\n");
	row = 1;
	col = m;
	matr = new double[m];
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, data[j - 1]);
}

matrix::matrix(int n, double* data) {
	if (n <= 0 || data == nullptr)
		throw("error: matrix m(int, double*);\n");
	row = n;
	col = 1;
	matr = new double[n];
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, data[i - 1]);
}

matrix::matrix(char* str) {
	if (str == nullptr)
		throw("error: matrix m(char*);\n");
	int len = strlen(str);
	double tmp = 0, flag = 0.1;
	int m = 1, n = 1;
	int flagrow = 0, flagcol = 0;
	for (int i = 0; i < len; i++) {
		if (str[i] == '{')
			flagrow++;
		if (str[i] == ',')
			flagcol++;
	}
	row = flagrow - 1;
	col = (flagcol + 1) / row;
	matr = new double[row * col];
	for (int i = 0; i < len; i++)
		if (str[i] == ' ') continue;
		else if (str[i] == '{') continue;
		else if (str[i] == ',') continue;
		else if (str[i] >= '0' && str[i] <= '9') {
			while ((str[i] >= '0' && str[i] <= '9') || str[i] == '.') {
				if (str[i] == ' ') continue;
				else if (str[i] >= '0' && str[i] <= '9') {
					tmp = tmp * 10 + (double)(str[i] - '0');
					i++;
				}
				else if (str[i] == '.') {
					i++;
					while (str[i] >= '0' && str[i] <= '9' || str[i] == ' ') {
						if (str[i] == ' ') {
							i++;
							continue;
						}
						tmp = tmp + flag * (double)(str[i] - '0');
						flag *= 0.1;
						i++;
					}
					break;
				}
				else throw("error: matrix m(char*), wrong char\n");
			}
			Set(m, n, tmp);
			n++;
			tmp = 0;
			flag = 0.1;
			i--;
			continue;
		}
		else if (str[i] == '}') {
			m++;
			n = 1;
			continue;
		}
		else throw("error: matrix m(char*), wrong char\n");
}

int matrix::rows() const { return row; }

int matrix::columns() const { return col; }

double* matrix::getmatr() { return matr; }

void matrix::Set(int nX, int nY, double r) { matr[col * (nX - 1) + nY - 1] = r; }

double matrix::Get(int nX, int nY)const { return matr[col * (nX - 1) + nY - 1]; }

ostream& operator<<(ostream& os, matrix& m) {
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			cout.width(10);
			cout.precision(3);
			os << m.matr[i * m.col + j] << " ";
		}
		os << endl;
	}
	return os;
}

void matrix::InitRnd() {
	for (int i = 0; i < row; i++) for (int j = 0; j < col; j++)
		matr[i * col + j] = rand() % 10;
}

void matrix::print() {
	if (matr == nullptr)
		throw("error print, matr is null");
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++)
			cout << matr[i * col + j] << " ";
		cout << endl;
	}
}

matrix operator * (matrix& m1, double data) {
	matrix mnew(m1.rows(), m1.columns());
	for (int i = 1; i <= m1.rows(); i++)
		for (int j = 1; j <= m1.columns(); j++)
			mnew.Set(i, j, data * m1.Get(i, j));
	return mnew;
}

matrix& matrix::operator *= (double data) {
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, data * Get(i, j));
	return *this;
}

matrix operator * (matrix& m1, double* v) {
	int len = 0;
	for (int i = 0; v[i] != '\0'; i++)
		len++;
	if (m1.columns() != len)
		throw("error: *\n");
	matrix mnew(m1.rows(), 1);
	int i, j;
	double sum;
	for (i = 1; i <= m1.rows(); i++) {
		sum = 0;
		for (j = 1; j <= m1.columns(); j++)
			sum += m1.Get(i, j) * v[j - 1];
		mnew.Set(i, 1, sum);
	}
	return mnew;
}

matrix& matrix::operator *= (double* v) {
	int len = 0;
	for (int i = 0; v[i] != '\0'; i++)
		len++;
	if (columns() != len)
		throw("error: *\n");
	double* tmp = new double[row * col];
	int i, j;
	double sum;
	int flag = col;
	for (i = 0; i < row * col; i++) tmp[i] = matr[i];
	delete[] matr;
	matr = new double[row];
	col = 1;
	for (i = 1; i <= row; i++) {
		sum = 0;
		for (j = 1; j <= flag; j++)
			sum += tmp[(i - 1) * flag + j - 1] * v[j - 1];
		Set(i, 1, sum);
	}
	delete[] tmp;
	return *this;
}

matrix operator + (matrix& m1, matrix& m2) {
	if (!(m1.rows() == m2.rows() && m1.columns() == m2.columns()))	throw("error + \n");
	matrix mnew(m1.rows(), m1.columns());
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		mnew.Set(i, j, m1.Get(i, j) + m2.Get(i, j));
	return mnew;
}

matrix& matrix::operator += (matrix& m1) {
	if (!(rows() == m1.rows() && columns() == m1.columns())) throw("error += \n");
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		Set(i, j, Get(i, j) + m1.Get(i, j));
	return *this;
}

matrix operator - (matrix& m1, matrix& m2) {
	if (!(m1.rows() == m2.rows() && m1.columns() == m2.columns())) throw("error - \n");
	matrix mnew(m1.rows(), m1.columns());
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		mnew.Set(i, j, m1.Get(i, j) - m2.Get(i, j));
	return mnew;
}

matrix& matrix::operator -= (matrix& m1) {
	if (!(rows() == m1.rows() && columns() == m1.columns())) throw("error -= \n");
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		Set(i, j, Get(i, j) - m1.Get(i, j));
}

matrix operator * (matrix& m1, matrix& m2) {
	if (m1.columns() != m2.rows()) throw("error: *\n");
	int i, j, k;
	double* D = new double[m2.columns()];
	matrix mnew(m1.rows(), m2.columns());
	for (k = 1; k <= m1.rows(); k++) {
		for (i = 0; i < m2.columns(); i++) D[i] = 0;
		for (i = 1; i <= m2.columns(); i++) for (j = 1; j <= m2.rows(); j++)
			D[i - 1] += m1.Get(k, j) * m2.Get(j, i);
		for (i = 1; i <= m2.columns(); i++) mnew.Set(k, i, D[i - 1]);
	}
	delete[] D;
	return mnew;
}

matrix& matrix::operator *= (matrix& m1) {
	if (columns() != m1.rows()) throw("error: *=\n");
	int i, j, k;
	double* D = new double[m1.columns()];
	double* tmp = new double[row * col];
	for (int i = 0; i < row * col; i++) tmp[i] = matr[i];
	delete[] matr;
	matr = new double[row * m1.columns()];
	col = m1.columns();
	for (k = 1; k <= row; k++) {
		for (i = 0; i < m1.columns(); i++) D[i] = 0;
		for (i = 1; i <= m1.columns(); i++) for (j = 1; j <= m1.rows(); j++)
			D[i - 1] += tmp[(k - 1) * m1.rows() + j - 1] * m1.Get(j, i);
		for (i = 1; i <= m1.columns(); i++) Set(k, i, D[i - 1]);
	}
	delete[] tmp;
	delete[] D;
	return *this;
}

matrix& matrix::operator - () {
	if (matr == nullptr) throw("error: -m\n");
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, Get(i, j) - 1);
	return *this;
}

bool operator == (matrix& m1, matrix& m2) {
	if (m1.rows() != m2.rows() || m1.columns() != m2.columns()) return false;
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		if (m1.Get(i, j) == m2.Get(i, j)) continue;
		else return false;
	return true;
}

bool operator != (matrix& m1, matrix& m2) {
	if (m1.rows() != m2.rows() || m1.columns() != m2.columns()) return true;
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		if (m1.Get(i, j) != m2.Get(i, j)) return true;
		else continue;
	return false;
}

matrix operator | (matrix& m1, matrix& m2) {
	matrix mnew(m1.rows(), m1.columns() + m2.columns());
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		mnew.Set(i, j, m1.Get(i, j));
	for (int i = 1; i <= m2.rows(); i++) for (int j = 1 + m1.columns(); j <= m1.columns() + m2.columns(); j++)
		mnew.Set(i, j, m2.Get(i, j - m1.columns()));
	return mnew;
}

matrix operator / (matrix& m1, matrix& m2) {
	matrix mnew(m1.rows() + m2.rows(), m1.columns());
	for (int i = 1; i <= m1.rows(); i++) for (int j = 1; j <= m1.columns(); j++)
		mnew.Set(i, j, m1.Get(i, j));
	for (int i = 1 + m1.rows(); i <= m1.rows() + m2.rows(); i++) for (int j = 1; j <= m2.columns(); j++)
		mnew.Set(i, j, m2.Get(i - m1.rows(), j));
	return mnew;
}

void matrix::Divid(double R, int nVerh) {
	for (int j = 1; j <= col; j++) Set(nVerh, j, Get(nVerh, j) / R);
}

void matrix::Minus(double R, int nNizz, int nVerh) {
	for (int j = 1; j <= col; j++) Set(nNizz, j, Get(nNizz, j) - R * Get(nVerh, j));
}

void matrix::Tuda() {
	for (int j = 1; j < row; j++) {
		Divid(Get(j, j), j);
		for (int i = j; i < col; i++)
			Minus(Get(i + 1, j), i + 1, j);
	}
	Divid(Get(row, col), row);
}

void matrix::Obratno() {
	for (int j = row; j > 0; j--)
		for (int i = j - 1; i > 0; i--)
			Minus(Get(i, j), i, j);
}

void matrix::inverse() {
	if (row != col) throw("error: inverse\n");
	matrix mnew = matrix::identity(row);
	double tmp;
	for (int j = 1; j < row; j++) {
		tmp = Get(j, j);
		Divid(tmp, j);
		mnew.Divid(tmp, j);
		for (int i = j; i < col; i++) {
			tmp = Get(i + 1, j);
			Minus(tmp, i + 1, j);
			mnew.Minus(tmp, i + 1, j);
		}
	}
	tmp = Get(row, col);
	Divid(tmp, row);
	mnew.Divid(tmp, row);
	for (int j = row; j > 0; j--)
		for (int i = j - 1; i > 0; i--) {
			tmp = Get(i, j);
			Minus(tmp, i, j);
			mnew.Minus(tmp, i, j);
		}
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, mnew.Get(i, j));
}

double matrix::determinant() {
	if (row != col)
		throw("error: determinant, don't have determinant\n");
	double* tmp = new double[row * col];
	for (int i = 0; i < row * col; i++)
		tmp[i] = matr[i];
	double flag = 1;
	for (int j = 1; j < row; j++) {
		flag *= Get(j, j);
		Divid(Get(j, j), j);
		for (int i = j; i < col; i++)
			Minus(Get(i + 1, j), i + 1, j);
	}
	flag *= Get(row, col);
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		Set(i, j, tmp[(i - 1) * col + j - 1]);
	delete[] tmp;
	return flag;
}

void matrix::transpose() {
	int t = row;
	row = col;
	col = t;

	double* temp = new double[row * col];
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			temp[row * (j - 1) + i - 1] = Get(i, j);
		}
	}

	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			Set(i, j, temp[row * (i - 1) + j - 1]);
		}
	}
	delete[] temp;
}

matrix matrix::solve() {
	if (row != col - 1)
		throw("error: solve\n");
	matrix ml(row, col - 1);
	matrix mr(row, 1);
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col - 1; j++)
		ml.Set(i, j, Get(i, j));
	for (int i = 1; i <= row; i++)
		mr.Set(i, 1, Get(i, col));
	int mlrow = ml.rows(), mlcol = ml.columns();
	double tmp;
	for (int j = 1; j < mlrow; j++) {
		tmp = ml.Get(j, j);
		ml.Divid(tmp, j);
		mr.Divid(tmp, j);
		for (int i = j; i < mlcol; i++) {
			tmp = ml.Get(i + 1, j);
			ml.Minus(tmp, i + 1, j);
			mr.Minus(tmp, i + 1, j);
		}
	}
	tmp = ml.Get(mlrow, mlcol);
	ml.Divid(tmp, mlrow);
	mr.Divid(tmp, mlrow);
	for (int j = mlrow; j > 0; j--)
		for (int i = j - 1; i > 0; i--) {
			tmp = ml.Get(i, j);
			ml.Minus(tmp, i, j);
			mr.Minus(tmp, i, j);
		}
	return mr;
}


void matrix::Init_matr() {
	for (double i = 1; i <= row; i++) {
		for (double j = 1; j <= col; j++)
			Set(i, j, min(i, j) / max(i, j));
	}
}

void matrix::Input_f(matrix& m) {
	double temp = 0;
	for (int i = 1; i <= m.row; i++) {
		for (int j = 1; j <= m.col; j++) {
			temp += m.Get(i, j);
		}
		//cout << temp << endl;
		Set(i, 1, temp);
		temp = 0;
	}
}

void matrix::Init_x() {
	for (int i = 1; i <= row; i++) {
		Set(i, 1, 0);
	}
}

void matrix::ToB(matrix& A, double omega) {
	matrix A1(row, col);
	matrix A2(row, col);
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			if (i == j) {
				A1.Set(i, j, (double)A.Get(i, j) / 2);
				A2.Set(i, j, (double)A.Get(i, j) / 2);
			}
			else if (i < j) {
				A1.Set(i, j, A.Get(i, j));
				A2.Set(i, j, 0);
			}
			else {
				A1.Set(i, j, 0);
				A2.Set(i, j, A.Get(i, j));
			}
		}
	}

	matrix E = matrix::identity(row);
	A1 *= omega;
	A1 += E;
	A2 *= omega;
	A2 += E;
	matrix mtemp = A1 * A2;

	for (int i = 1; i <= row; i++)
		for (int j = 1; j <= col; j++)
			Set(i, j, mtemp.Get(i, j));
}

void matrix::Copy(matrix& m1) {
	row = m1.rows();
	col = m1.columns();
	if (matr != nullptr)
		delete[] matr;
	matr = new double[row * col];
	for (int i = 1; i <= row; i++)
		for (int j = 1; j <= col; j++)
			Set(i, j, m1.Get(i, j));
}

double matrix::GetError() {
	double e = 0;
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			e += pow(Get(i, j), 2);
		}
	}
	return sqrt(e);
}

void Jacobi(matrix& A) {
	double epsilon = 0.0001;
	int n = A.rows();
	double max_aij = fabs(A.Get(1, 2));
	int flagi = 1, flagj = 2;
	matrix P = matrix::identity(n);
	int num_iteration = 0;
	//cout << A << endl;
	while (1) {
		//cout << "P" << P << endl;
		//cout << "A" << A << endl;
		/*
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				if (i == j)
					P.Set(i, j, 1);
				else
					P.Set(i, j, 0);
			}
		}*/
		
		max_aij = fabs(A.Get(1, 2));
		flagi = 1;
		flagj = 2;
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				if (i != j && (fabs(A.Get(i, j)) > max_aij)) {
					max_aij = fabs(A.Get(i, j));
					flagi = i;
					flagj = j;
				}
			}
		}
		//cout << flagi << flagj << endl;
		//cout << "test" << endl;

		if (max_aij < epsilon) {
			cout << "Number of iteration:   " << num_iteration << endl << endl;
			cout << "Feature vector:" << endl;
			cout << P << endl;
			return;
		}
		
		/*
		double aii = A.Get(flagi, flagi);
		double ajj = A.Get(flagj, flagj);
		double aij = A.Get(flagi, flagj);
		double aji = A.Get(flagj, flagi);
		double sin_theta = sin(theta);
		double cos_theta = cos(theta);
		*/
		

		
		/*double theta;
		if (A.Get(flagi, flagi) == A.Get(flagj, flagj))
			theta = PI / 4;
		else
			theta = 0.5 * atan(2 * A.Get(flagi, flagj) / (A.Get(flagi, flagi) - A.Get(flagj, flagj)));
		*/
		double angle = 0.5 * atan(2 * A.Get(flagi, flagj) / (A.Get(flagi, flagi) - A.Get(flagj, flagj)));

		double aii = A.Get(flagi, flagi);
		double ajj = A.Get(flagj, flagj);
		double aij = A.Get(flagi, flagj);
		double sin_theta = sin(angle);
		double cos_theta = cos(angle);
		double sin_2angle = sin(2 * angle);
		double cos_2angle = cos(2 * angle);


		A.Set(flagi, flagi, aii * cos_theta * cos_theta + ajj * sin_theta * sin_theta + aij * sin_2angle);
		A.Set(flagj, flagj, aii * sin_theta * sin_theta + ajj * cos_theta * cos_theta - aij * sin_2angle);
		A.Set(flagi, flagj, 0.5 * (ajj - aii) * sin_2angle + aij * cos_2angle);
		A.Set(flagj, flagi, A.Get(flagi, flagj));
		for (int k = 1; k <= n; k++) {
			if (k != flagi && k != flagj) {
				double arowk = A.Get(flagi, k);
				double acolk = A.Get(flagj, k);
				A.Set(flagi, k, arowk * cos_theta + acolk * sin_theta);
				A.Set(k, flagi, A.Get(flagi, k));
				A.Set(flagj, k, acolk * cos_theta - arowk * sin_theta);
				A.Set(k, flagj, A.Get(flagj, k));
			}
		}
		
		/*
		A.Set(flagi, flagi,
			aii * cos_theta * cos_theta + ajj * sin_theta * sin_theta + aij * sin_theta * cos_theta + aji * sin_theta * cos_theta);
		A.Set(flagj, flagj,
			aii * sin_theta * sin_theta + ajj * cos_theta * cos_theta - aij * sin_theta * cos_theta - aji * sin_theta * cos_theta);
		A.Set(flagj, flagi,
			-aii * sin_theta * cos_theta + ajj * sin_theta * cos_theta - aij * sin_theta * sin_theta + aji * cos_theta * cos_theta);
		A.Set(flagi, flagj,
			-aii * sin_theta * cos_theta + ajj * sin_theta * cos_theta + aij * cos_theta * cos_theta - aji * sin_theta * sin_theta);
			
		cout << A << endl;
	
		
		A *= P;
		P.transpose();
		P *= A;
		A.Copy(P);
			*/

		//cout << A << endl;
		/*
		for (int k = 1; k <= n; k++) {
			if (k != flagi && k != flagj) {
				double arowk = A.Get(flagi, k);
				double acolk = A.Get(flagj, k);
				A.Set(flagi, k, arowk * cos_theta + acolk * sin_theta);
				A.Set(k, flagi, A.Get(flagi, k));
				A.Set(flagj, k, acolk * cos_theta - arowk * sin_theta);
				A.Set(k, flagj, A.Get(flagj, k));
			}
		}
		*/
		//cout << A << endl;
		
		if (num_iteration == 0) {
			P.Set(flagi, flagi, cos_theta);
			P.Set(flagi, flagj, -sin_theta);
			P.Set(flagj, flagi, sin_theta);
			P.Set(flagj, flagj, cos_theta);
		}
		else {
			for (int k = 1; k <= n; k++) {
				double pki = P.Get(k, flagi);
				double pkj = P.Get(k, flagj);;
				P.Set(k, flagi, pki * cos_theta + pkj * sin_theta);
				P.Set(k, flagj, pkj * cos_theta - pki * sin_theta);
			}
		}
		
		num_iteration++;
		//cout << ite_num << endl;

	}
}

int main() {
	int n = 10;
	matrix A(n, n);
	A.Init_matr();
	cout << "Matrix is:" << endl << A << endl;

	clock_t start, end;
	start = clock();

	Jacobi(A);

	end = clock();

	//cout << A << endl;

	cout << "Eigenvalues:" << endl;
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if (i == j) {
				cout.width(10);
				cout.precision(3);
				cout << A.Get(i, j) << " ";
			}
		}
	}

	cout << endl << endl << "The runtime is:" << (double)(end - start) / CLOCKS_PER_SEC << "s";

	return 0;
}

