
#include "NormalToeplitz.h"

using namespace std;

// generate a fbm acf vector
void fbm_acf(int N, double a, double* acf1) {
	double* acf2;
	acf2 = new double[N + 1];
	for (int ii = 0; ii <= N; ++ii) {
		acf2[ii] = pow(ii, a);
	}
	acf1[0] = acf2[1] - acf2[0];
	for (int ii = 1; ii < N; ++ii) {
		acf1[ii] = (acf2[ii + 1] + acf2[ii - 1]) / 2 - acf2[ii];
	}
}

// generate a arithmetic sequence
void z_func(int N, double a, double* z) {
	for (int ii = 0; ii < N; ++ii) {
		z[ii] = (double)(1 + ii) * a / 10;
	}
}

// print a vector
void printVector(double* x, int n) {
	for (int ii = 0; ii < n; ++ii) {
		std::cout << x[ii] << " ";
	}
	std::cout << std::endl;
	return;
}

// print a complex vector
void printComplex(std::complex<double>* x, int n) {
	for (int ii = 0; ii < n; ++ii) {
		std::cout << x[ii] << " ";
	}
	std::cout << std::endl << std::endl;
	return;
}

int main() {

	int N = 10;
	int p = 3;

	NormalToeplitz* Ta = new NormalToeplitz(N, p);
	double* acf1 = new double[N];
	double* z = new double[N];
	double* dacf = new double[N * p];
	double* dz = new double[N * p];
	double* d2acf = new double[p * p * N];
	double* d2z = new double[p * p * N];

	double* dl = new double[p];
	double* d2l = new double[p * p];

	fbm_acf(N, 0.8, acf1);
	z_func(N, 1, z);
	
	for (int ii = 0; ii < p; ++ii) {
		fbm_acf(N, (double)(ii + 1) / 10, &dacf[ii * N]);
		z_func(N, (double)(ii + 1) / 10, &dz[ii * N]);
	}
	
	for (int ii = 0; ii < p; ++ii) {
		for (int jj = 0; jj < p; ++jj) {
			fbm_acf(N, (double)(ii + jj + 2) / 10, &d2acf[(ii * p + jj) * N]);
			z_func(N, (double)(ii + jj + 2) / 10, &d2z[(ii * p + jj) * N]);
		}
	}

	cout << "print the log density: " << endl;
	cout << Ta->logdens(z, acf1) << endl << endl;

	cout << "print the gradient: " << endl;
	Ta->grad(dl, z, dz, acf1, dacf);
	printVector(dl, p);
	cout << endl << endl;

	cout << "print the hessian: " << endl;
	Ta->hess(d2l, z, dz, d2z, acf1, dacf, d2acf);
	printVector(d2l, p * p);
	cout << endl << endl;

	cout << "print the full gradient with respect to z: " << endl;
	Ta->grad_full(dz, dacf, z, acf1);
	printVector(dz, N);
	cout << endl << endl;

	cout << "print the full gradient with respect to acf: " << endl;
	printVector(dacf, N);
	cout << endl << endl;

	delete[] d2acf;
	delete[] d2z;
	delete[] d2l;
	delete[] dacf;
	delete[] dz;
	delete[] dl;
	delete[] acf1;
	delete[] z;

	return 0;
}