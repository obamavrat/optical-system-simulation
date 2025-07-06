#include <cmath>
#include <iostream>

class vector2 {

public:
	double x = 0;
	double y = 0;
	void set(double x, double y) {
		this->x = x;
		this->y = y;
	}
	double sqrNorm() {
		double res = (this->x) * (this->x) + (this->y) * (this->y);
		return res;
	}
	double Norm() {
		return sqrt(this->sqrNorm());
	};
};

class vector {
public:
	double component[3];
	double &x= component[0];
	double &y= component[1];
	double &z= component[2];
	vector operator=(vector a) {
		for (int i = 0;i < 3;i++) {
			this->component[i] = a.component[i];
		}
	}
//vector() {
//			component[0] = 0;
//			component[1] = 0;
//			component[2] = 0;
//			x = 0;
//			y = 0;
//			z = 0;
//		 }
	double sqrNorm() {
		return (this->x) * (this->x) + (this->y) * (this->y) + (this->z) * (this->z);
	};
	double Norm() {
		return sqrt(this->sqrNorm());
	};
	vector operator + (vector a) {
		vector b;
		b.x = a.x + this->x;
		b.y = a.y + this->y;
		b.z = a.z + this->z;
		return b;
	};
	vector operator - (vector a) {
		vector b;
		b.x = this->x - a.x;
		b.y = this->y - a.y;
		b.z = this->z - a.z;
		return b;
	};
	vector operator * (double a) {
		vector b;
		b.x = this->x * a;
		b.y = this->y * a;
		b.z = this->z * a;
		return b;
	};
	vector operator / (double a) {
		vector b;
		b.x = this->x / a;
		b.y = this->y / a;
		b.z = this->z / a;
		return b;
	};
	double operator *(vector a) {
		double b;
		b = this->x * a.x + this->y * a.y + this->z * a.z;
		return b;
	}
	void normalize() {
		double N = this->Norm();
		this->x = this->x / N;
		this->y = this->y / N;
		this->z = this->z / N;
	}
	void setUnitVector(double phi, double psi) {
		this->x = cos(phi) * sin(psi);
		this->y = sin(phi) * sin(psi);
		this->z = cos(psi);
	}
	void set(vector a) {
		this->x = a.x;
		this->y = a.y;
		this->z = a.z;
	}
	void set(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	};
	void zero() {
		this->set(0, 0, 0);
	}
	vector2 proectionXY() {
		vector2 res;
		res.set(this->x, this->y);
		return res;
	}
	vector X (vector a) {
		vector res;
		res.x = this->y * a.z - this->z * a.y;
		res.y = a.x * this->z - a.z * this->x;
		res.z = this->x * a.y - this->y * a.x;
		return res;
	}
	void print() {
		std::cout << this->x << '\t' << this->y << '\t' << this->z << std::endl;
		
	}
};
int mod(int a, int m) {
	if (a < 0) {
		int b = (-a) % m;
		return (m - b) % m;
	}
	else {
		return a % m;
	}
}
class subMatrix {
public:
	double value[2][2];
	double det() {
		double res = this->value[0][0] * this->value[1][1] - this->value[1][0] * this->value[0][1];
		return res;
	}
	subMatrix operator * (double a) {
		subMatrix res;
		for (int m = 0;m < 2;m++) {
			for (int n = 0;n < 2;n++) {
				res.value[m][n] = this->value[m][n] * a;
			}
		}
		return res;
	}
};
class matrix {
public:
	double value[3][3];
	void set(double a, double b, double c,
		     double e, double f, double g,
		     double i, double j, double k) {
		this->value[0][0] = a;
		this->value[1][0] = b;
		this->value[2][0] = c;
		this->value[0][1] = e;
		this->value[1][1] = f;
		this->value[2][1] = g;
		this->value[0][2] = i;
		this->value[1][2] = j;
		this->value[2][2] = k;
	}
	double getValue(int i, int j) {
		int n = mod(i, 3);
		int m = mod(j, 3);
		return this->value[n][m];
	}
	double det() {
		double A1 = this->value[0][0] * (this->value[1][1] * this->value[2][2] - this->value[1][2] * this->value[2][1]);
		double A2 = this->value[0][1] * (this->value[0][1] * this->value[2][2] - this->value[2][1] * this->value[0][2]);
		double A3 = this->value[0][2] * (this->value[0][1] * this->value[1][2] - this->value[1][1] * this->value[0][2]);
		double res = A1 - A2 + A3;
		return res;
	}
	vector getCol(int i) {
		int j = mod(i, 3);
		vector res;
		res.set(this->value[i][0], this->value[i][1], this->value[i][2]);
		return res;
	}
	void setCol(int i, vector a) {
		int j = mod(i, 3);
		this->value[i][0] = a.x;
		this->value[i][1] = a.y;
		this->value[i][2] = a.z;
	}
	vector getLine(int i) {
		int j = mod(i, 3);
		vector res;
		res.set(this->value[0][i], this->value[1][i], this->value[2][i]);
		return res;
	}
	void setLine(int i, vector a) {
		int j = mod(i, 3);
		this->value[0][i] = a.x;
		this->value[1][i] = a.y;
		this->value[2][i] = a.z;
	}
	void setZero() {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				this->value[i][j] = 0;
			}
		}
	}
	void print() {
		vector a = this->getLine(0);
		vector b = this->getLine(1);
		vector c = this->getLine(2);
		a.print();
		b.print();
		c.print();
	}
	matrix transp() {
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = this->value[j][i];
			}
		}
	}
	matrix operator* (matrix a) {
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = 0;
				for (int k = 0;k < 3;k++) {
					res.value[i][j] += this->value[i][k] * a.value[k][j];
				}
				
			}
		}
		return res;
	}
	matrix operator* (double a) {
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = this->value[i][j]*a;
			}
		}
		return res;
	}
	matrix operator+ (matrix a) {
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = this->value[i][j] + a.value[i][j];
			}
		}
		return res;
	}
	matrix operator/ (double a) {
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = this->value[i][j] / a;
			}
		}
		return res;
	}
	vector operator* (vector a) {
		vector res;
		for (int i = 0; i < 3; i++) {
			res.component[i] = 0;
			for (int j = 0; j < 3;j++) {
				res.component[i] += this->value[i][j] * a.component[j];
			}
		}
		return res;
	}
	subMatrix minor(int i, int j) {
		subMatrix res;
		for (int m = 0;m < 2;m++) {
			for (int n = 0;n < 2;n++) {
				res.value[m][n] = this->getValue(i + m + 1, j + n + 1);
			}
		}
		return res;
	}
	matrix inv() {
		double det = this->det();
		matrix res;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3;j++) {
				res.value[i][j] = this->minor(i, j).det() / det * pow(-1, i + j);
			}
		}
	}


};


