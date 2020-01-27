/// @file Toeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef NormalToeplitz_h
#define NormalToeplitz_h 1

#include "Toeplitz.h"
# define logPI 0.49714987269L // log of PI value

/// Compute \f$\text{tr}(L U)\f$, where \f$L\f$ is a lower triangular Toeplitz
/// matrix and \f$U\f$ is a upper triangular Toeplitz matrix. It is an
/// \f$O(n)\f$ algorithm.
///
/// @param[in] acf1 Real vector of length `n`, first column of \f$L\f$.
/// @param[in] acf2 Real vector of length `n`, first row of \f$U\f$.
/// @param[in] n Integer indicating the length of input.
inline double trace_LU(double* acf1, double* acf2, int n) {
	double trace = 0;
	for (int ii = 0; ii < n; ++ii) {
		trace += (n - ii) * acf1[ii] * acf2[ii];
	}
	return trace;
}

/// Function for vector multiplication
///
/// @param[in] v1 First input vector.
/// @param[in] v2 Second input vector.
/// @param[in] n Length of input vector.
/// @param[out] double number that equals v1 * v2.
inline double vecProd(double* v1, double* v2, int n) {
	double ans = 0;
	for (int ii = 0; ii < n; ++ii) {
		ans += v1[ii] * v2[ii];
	}
	return ans;
}

/// Class for Computation involving Toeplitz matrix.
///
/// Model: z ~ N_(0, V), where V is a Toeplitz matrix with parameter theta
class NormalToeplitz: public Toeplitz {
	int p; // number of parameters theta
	
	// storages for temporary vectors
	double* vec1;
	double* vec2;
	double* vec3;
	double* vec4;
	double* phi;

	// upper lower triangular Toeplitz matrices
	double* U1;
	double* U2;
	std::complex<double>* U1_fft;
	std::complex<double>* U2_fft;

public:
	/// Constructor.
	NormalToeplitz(int, int);

	/// Destructor
	~NormalToeplitz();

	/// Size
	int size();

	/// Dimension
	int dim();

	/// Function for a trace computation
	double trace_deriv(double*);

	/// Function for another trace computation
	double trace_hess(double*, double*);

	/// Log-Density.
	double logdens(double* z, double* acf);

	/// Log-likelihood Gradient.
	void grad(double* dldt, double* z, double* dzdt, double* acf, double* dacfdt);

	/// Hessian matrix of Log-likelihood.
	void hess(double* d2ldt, double* z, double* dzdt, double* d2zdt,
		double* acf, double* dacfdt, double* d2acfdt);

	/// Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
	void grad_full(double* dldz, double* dldacf, double* z, double* acf);

	/*
	/// Overloaded Log-Density.
	double logdens(double* z, Toeplitz* Tz);

	/// Overloaded Log-likelihood Gradient.
	void grad(double* dldt,
		double* z, double* dzdt,
		Toeplitz* Tz, double* dacfdt);

	/// Overloaded function for the Hessian matrix of Log-likelihood.
	void hess(double* d2ldt,
		double* z, double* dzdt, double* d2zdt,
		Toeplitz* Tz, double* dacfdt, double* d2acfdt);

	/// Overloaded function for the Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
	void grad_full(double* dldz, double* dldacf,
		double* z, Toeplitz* Tz);
	*/
};

/// Constructor.
///
/// @param N_ length of observation.
/// @param p_ number of parameters.
/// @param hasToep whether Toeplitz object Tz is imported from outside or generated within the class.
inline NormalToeplitz::NormalToeplitz(int N_, int p_) : Toeplitz(N_) {
	p = p_;

	vec1 = new double[n];
	vec2 = new double[n];
	vec3 = new double[n];
	vec4 = new double[n];
	phi = new double[n];	
	U1 = new double[2 * n];
	U2 = new double[2 * n];
	U1_fft = new std::complex<double>[2 * n];
	U2_fft = new std::complex<double>[2 * n];
}

/// Destructor
inline NormalToeplitz::~NormalToeplitz() {
	delete vec1;
	delete vec2;
	delete vec3;
	delete vec4;
	delete phi;

	delete[] U1;
	delete[] U2;
	delete[] U1_fft;
	delete[] U2_fft;
}

inline int NormalToeplitz::size() {
	return n;
}

inline int NormalToeplitz::dim() {
	return p;
}

inline double NormalToeplitz::trace_deriv(double* acf0) {
	double trace;
	if (n > 1) { // GSchur algorithm only supports N > 1 case.
		if (!has_solve) {
			solve_setup();
		}

		std::copy(acf0, acf0 + n, U1);
		std::fill(U1 + n, U1 + 2 * n, 0);
		V1->fft(U1_fft, U1); // U1
		std::fill(U2, U2 + 2 * n, 0);
		std::copy(acf0 + 1, acf0 + n, U2 + 1);
		V1->fft(U2_fft, U2);  // U2

		// tr{U1'L1L1'U1}
		Complex_Mult(y_fft, L1_fft, U1_fft, n + 1);
		V1->ifft(y, y_fft);
		trace = trace_LU(y, y, n);

		// tr{U1'L2L2'U1}
		Complex_Mult(y_fft, L2_fft, U1_fft, 2 * (n / 2 + 1));
		V1->ifft(y, y_fft);
		trace -= trace_LU(y, y, n);
		// tr{U2'L1L1'U2}
		Complex_Mult(y_fft, L1_fft, U2_fft, 2 * (n / 2 + 1));
		V1->ifft(y, y_fft);
		trace -= trace_LU(y, y, n);
		// tr{U2'L2L2'U2}
		Complex_Mult(y_fft, L2_fft, U2_fft, 2 * (n / 2 + 1));
		V1->ifft(y, y_fft);
		trace += trace_LU(y, y, n);

		// trace
		trace /= Gs->delta[0];
		trace /= acf0[0];
	}

	else {
		// N = 1 case.
		trace = acf0[0] / acf[0];
	}
	return trace;
}

/// Function for another trace computation
inline double NormalToeplitz::trace_hess(double* acf1, double* acf2) {
	double trace2 = 1;
	if (n > 1) {
		// GSchur algorithm only supports N > 1 case.
		if (!has_solve) {
			solve_setup();
		}
		
		// we store the negative derivative of delta in vector phi, where phi = solve(acf) * toep(acf1) * delta
		product(phi, Gs->delta, acf1);
		solve(phi, phi);
		trace2 = trace_deriv(acf2) * phi[0] * -1;

		double k1;
		std::copy(phi, phi + n, x);
		std::fill(x + n, x + 2 * n, 0);
		V1->fft(x_fft, x);
		std::copy(acf2, acf2 + n, z);
		std::fill(z + n, z + 2 * n, 0);
		V1->fft(z_fft, z);
		Complex_Mult(U1_fft, x_fft, z_fft, n + 1);
		V1->ifft(U1, U1_fft);
		std::copy(Gs->delta, Gs->delta + n, y);
		std::fill(y + n, y + 2 * n, 0);
		V1->fft(y_fft, y);
		Complex_Mult(U2_fft, y_fft, z_fft, n + 1);
		V1->ifft(U2, U2_fft);
		k1 = trace_LU(U1, U2, n) / acf2[0];
		z[0] = 0;
		V1->fft(z_fft, z);
		Complex_Mult(U1_fft, x_fft, z_fft, n + 1);
		V1->ifft(U1, U1_fft);
		Complex_Mult(U2_fft, y_fft, z_fft, n + 1);
		V1->ifft(U2, U2_fft);
		k1 -= trace_LU(U1, U2, n) / acf2[0];

		double k2;
		x[0] = 0;
		std::reverse(x + 1, x + n);
		V1->fft(x_fft, x);
		z[0] = acf2[0];
		V1->fft(z_fft, z);
		Complex_Mult(U1_fft, x_fft, z_fft, n + 1);
		V1->ifft(U1, U1_fft);
		y[0] = 0;
		std::reverse(y + 1, y + n);
		V1->fft(y_fft, y);
		Complex_Mult(U2_fft, y_fft, z_fft, n + 1);
		V1->ifft(U2, U2_fft);
		k2 = trace_LU(U1, U2, n) / acf2[0];
		z[0] = 0;
		V1->fft(z_fft, z);
		Complex_Mult(U1_fft, x_fft, z_fft, n + 1);
		V1->ifft(U1, U1_fft);
		Complex_Mult(U2_fft, y_fft, z_fft, n + 1);
		V1->ifft(U2, U2_fft);
		k2 -= trace_LU(U1, U2, n) / acf2[0];

		trace2 += 2 * (k1 - k2);
		trace2 /= Gs->delta[0];
		
	}
	else {
		// N = 1 case.
		trace2 = acf1[0] * acf2[0] / acf[0] / acf[0];
	}
	return trace2;
}

/// Log-Density.
///
/// @param[in] z length-N vector of observation.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[out] double number of the log-density.
inline double NormalToeplitz::logdens(double* z, double* acf) {
	double ldens = 0;
	setAcf(acf); // Tz = Toeplitz(acf)
	solve(vec1, z); // vec1 = Tz^{-1} * z
	ldens = vecProd(z, vec1, n); // ldens = t(z) * Tz^{-1} * z
	ldens += logDet() + n * logPI;
	ldens *= -0.5;
	return ldens;
}


/// Log-likelihood Gradient.
///
/// @param[out] dldt length-p vector of the gradient.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[in] dacfdt length-N*p vector of observation derivatives, elememts dacfdt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
inline void NormalToeplitz::grad(double* dldt, double* z, double* dzdt, 
	double* acf, double* dacfdt) {
	setAcf(acf); // Tz = Toeplitz(acf)
	solve(vec1, z); // vec1 = Tz^{-1} * z
	for (int ii = 0; ii < p; ++ii) {
		product(vec2, vec1, &dacfdt[ii * n]);
		dldt[ii] = vecProd(vec1, vec2, n) / 2;
		dldt[ii] -= vecProd(&dzdt[ii * n], vec1, n);
		dldt[ii] -= trace_deriv(&dacfdt[ii * n]) / 2;
	}
}



/// Hessian matrix of Log-likelihood.
///
/// @param[out] d2ldt length-p*p vector of the Hessian, elememts d2ldt[0:p + ii*p] is the (ii+1)-th column of Hessian matrix.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] d2zdt length-N*p*p vector of observation derivatives, elememts d2zdt[0:N + (ii*p+jj)*N] is the second derivative of z with respect to the ii-th parameter and then the jj-th parameter.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[in] dacfdt length-N*p vector of observation derivatives, elememts dacfdt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
/// @param[in] d2acfdt length-N*p*p vector of observation derivatives, elememts d2acfdt[0:N + (ii*p+jj)*N] is the second derivative of acf with respect to the ii-th parameter and then the jj-th parameter.
inline void NormalToeplitz::hess(double* d2ldt,
	double* z, double* dzdt, double* d2zdt,
	double* acf, double* dacfdt, double* d2acfdt) {
	setAcf(acf);
	solve(vec1, z);
	double ans;
	std::fill(d2ldt, d2ldt + p * p, 0);
	for (int ii = 0; ii < p; ++ii) {
		for (int jj = 0; jj <= ii; ++jj) {
			product(vec4, vec1, &dacfdt[jj * n]);
			product(vec3, vec1, &dacfdt[ii * n]);
			ans = vecProd(&d2zdt[(ii * p + jj) * n], vec1, n);

			solve(vec2, vec4);
			ans -= vecProd(&dzdt[ii * n], vec2, n);
			ans += vecProd(vec3, vec2, n);
			solve(vec2, vec3);
			ans -= vecProd(&dzdt[jj * n], vec2, n);
			solve(vec2, &dzdt[jj * n]);
			ans += vecProd(&dzdt[ii * n], vec2, n);
			ans *= 2;
			
			product(vec2, vec1, &d2acfdt[(ii * p + jj) * n]);
			ans -= vecProd(vec1, vec2, n);
			ans += trace_deriv(&d2acfdt[(ii * p + jj) * n]);
			ans -= trace_hess(&dacfdt[ii * n], &dacfdt[jj * n]);

			d2ldt[ii * p + jj] = -ans / 2;
		}
	}

	if (p > 1) {
		for (int ii = 0; ii < p; ++ii) {
			for (int jj = ii + 1; jj < p; ++jj) {
				d2ldt[ii * p + jj] = d2ldt[jj * p + ii];
			}
		}
	}

}


/// Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
///
/// @param[out] dldz length-N vector of the log-density with respect to the observation z.
/// @param[out] dldacf length-N vector of the log-density with respect to the auto-covariance acf.
/// @param[in] z length-N vector of observation.
/// @param[in] acf length-N vector of auto-covariance.
inline void NormalToeplitz::grad_full(double* dldz, double* dldacf, double* z, double* acf) {
	setAcf(acf);

	// gradient with respect to z
	solve(dldz, z);
	for (int ii = 0; ii < n; ii++) {
		dldz[ii] = -dldz[ii];
	}

	// gradient with respect to acf
	solve(vec1, z);	
	vec2[0] = 1;
	std::fill(vec2 + 1, vec2 + n, 0);
	solve(vec3, vec2);
	double tau1 = vec3[0];
	std::fill(phi, phi + n, 0);
	phi[0] = vec1[0];
	product(dldacf, vec1, phi, vec1); // dldacf = upper.toep(Vz) %*% Vz = ip
	vec4[0] = 0;
	for (int ii = 1; ii < n; ++ii) {
		vec4[ii] = vec3[n - ii];
	}
	for (int ii = 0; ii < n; ++ii) {
		vec2[ii] = (n - ii) * vec3[ii];
	} // vec2 = (n:1 * tau)

	phi[0] = vec3[0];
	product(vec1, vec2, phi, vec3); // vec1 = upper.toep(tau) %*% (n:1 * tau) = tr

	for (int ii = 0; ii < n; ++ii) {
		vec2[ii] = (n - ii) * vec4[ii];
	} // vec2 = (n:1 * tau2)

	phi[0] = vec4[0];
	product(vec3, vec2, phi, vec4); // vec3 = upper.toep(tau2) %*% (n:1 * tau2)

	for (int ii = 0; ii < n; ++ii) {
		vec1[ii] -= vec3[ii];
		vec1[ii] /= tau1;
		dldacf[ii] -= vec1[ii];
	} // vec1 = (vec1 - vec3) / tau[1] = tr, dldacf = ip - tr

	dldacf[0] /= 2;
}


/*

//////////////////////////////// overloaded version //////////////////////////////////////////

/// Overloaded Log-Density.
///
/// @param[in] z length-N vector of observation.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[out] double number of the log-density.
inline double NormalToeplitz::logdens(double* z, Toeplitz* Tz) {
	double ldens = 0;
	// acf is already stored in Tz
	Tz->solveVec(vec1, z); // vec1 = Tz^{-1} * z
	ldens = vecProd(z, vec1, N_); // ldens = t(z) * Tz^{-1} * z
	ldens += Tz->logDet() + N_ * logPI;
	ldens *= -0.5;
	return ldens;
}

/// Overloaded Log-likelihood Gradient.
///
/// @param[out] dldt length-p vector of the gradient.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[in] dacfdt length-N*p vector of observation derivatives, elememts dacfdt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
inline void NormalToeplitz::grad(double* dldt, double* z, double* dzdt,
	Toeplitz* Tz, double* dacfdt) {
	// acf is already stored in Tz
	Tz->solveVec(vec1, z); // vec1 = Tz^{-1} * z
	double tmp;
	for (int ii = 0; ii < p_; ++ii) {
		T2_->setAcf(&dacfdt[ii * N_]);
		T2_->multVec(vec2, vec1);
		dldt[ii] = vecProd(vec1, vec2, N_) / 2;
		dldt[ii] -= vecProd(&dzdt[ii * N_], vec1, N_);
		dldt[ii] -= Tz->traceProd(&dacfdt[ii * N_]) / 2;
	}
}

/// Overloaded function for the Hessian matrix of Log-likelihood.
///
/// @param[out] d2ldt length-p*p vector of the Hessian, elememts d2ldt[0:p + ii*p] is the (ii+1)-th column of Hessian matrix.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] d2zdt length-N*p*p vector of observation derivatives, elememts d2zdt[0:N + (ii*p+jj)*N] is the second derivative of z with respect to the ii-th parameter and then the jj-th parameter.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[in] dacfdt length-N*p vector of observation derivatives, elememts dacfdt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
/// @param[in] d2acfdt length-N*p*p vector of observation derivatives, elememts d2acfdt[0:N + (ii*p+jj)*N] is the second derivative of acf with respect to the ii-th parameter and then the jj-th parameter.
inline void NormalToeplitz::hess(double* d2ldt,
	double* z, double* dzdt, double* d2zdt,
	Toeplitz* Tz, double* dacfdt, double* d2acfdt) {
	// acf is already stored in Tz
	Tz->solveVec(vec1, z);
	double ans;
	for (int ii = 0; ii < p_; ++ii) {
		for (int jj = 0; jj <= ii; ++jj) {
			T2_->setAcf(&dacfdt[jj * N_]);
			T2_->multVec(vec4, vec1);
			T2_->setAcf(&dacfdt[ii * N_]);
			T2_->multVec(vec3, vec1);

			ans = vecProd(&d2zdt[(ii * p_ + jj) * N_], vec1, N_);
			// temp0 = solve(Tz, temp2)
			Tz->solveVec(vec2, vec4);
			ans -= vecProd(&dzdt[ii * N_], vec2, N_);
			ans += vecProd(vec3, vec2, N_);
			// temp0 = solve(Tz, temp1)
			Tz->solveVec(vec2, vec3);
			ans -= vecProd(&dzdt[jj * N_], vec2, N_);
			// temp0 = solve(Tz, dZ[,jj])
			Tz->solveVec(vec2, &dzdt[jj * N_]);
			ans += vecProd(&dzdt[ii * N_], vec2, N_);
			ans *= 2;

			T2_->setAcf(&d2acfdt[(ii * p_ + jj) * N_]);
			T2_->multVec(vec2, vec1);
			ans -= vecProd(vec1, vec2, N_);
			ans += Tz->traceProd(&d2acfdt[(ii * p_ + jj) * N_]);
			ans -= Tz->traceDeriv(&dacfdt[ii * N_], &dacfdt[jj * N_]);

			d2ldt[ii * p_ + jj] = -ans / 2;
		}
	}

	if (p_ > 1) {
		for (int ii = 0; ii < p_; ++ii) {
			for (int jj = ii + 1; jj < p_; ++jj) {
				d2ldt[ii * p_ + jj] = d2ldt[jj * p_ + ii];
			}
		}
	}

}

/// Overloaded function for the Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
///
/// @param[out] dldz length-N vector of the log-density with respect to the observation z.
/// @param[out] dldacf length-N vector of the log-density with respect to the auto-covariance acf.
/// @param[in] z length-N vector of observation.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
inline void NormalToeplitz::grad_full(double* dldz, double* dldacf,
	double* z, Toeplitz* Tz) {
	// acf is already stored in Tz

	// gradient with respect to z
	Tz->solveVec(dldz, z);
	for (int ii = 0; ii < N_; ii++) {
		dldz[ii] = -dldz[ii];
	}

	// gradient with respect to acf
	Tz->solveVec(vec1, z); // vec1 = Vz

	vec2[0] = 1;
	std::fill(vec2 + 1, vec2 + N_, 0); // vec2 = [1,0,0,...,0]
	Tz->solveVec(vec3, vec2); // vec3 = tau
	double tau1 = vec3[0];

	T2_->setAcf(vec1);
	T2_->mult0Vec(dldacf, vec1); // dldacf = upper.toep(Vz) %*% Vz = ip

	vec4[0] = 0;
	for (int ii = 1; ii < N_; ++ii) {
		vec4[ii] = vec3[N_ - ii];
	} // vec4 = tau2

	for (int ii = 0; ii < N_; ++ii) {
		vec2[ii] = (N_ - ii) * vec3[ii];
	} // vec2 = (N_:1 * tau)

	T2_->setAcf(vec3);
	T2_->mult0Vec(vec1, vec2); // vec1 = upper.toep(tau) %*% (N_:1 * tau) = tr

	for (int ii = 0; ii < N_; ++ii) {
		vec2[ii] = (N_ - ii) * vec4[ii];
	} // vec2 = (N_:1 * tau2)
	T2_->setAcf(vec4);
	T2_->mult0Vec(vec3, vec2); // vec3 = upper.toep(tau2) %*% (N_:1 * tau2)

	for (int ii = 0; ii < N_; ++ii) {
		vec1[ii] -= vec3[ii];
		vec1[ii] /= tau1;
		dldacf[ii] -= vec1[ii];
	} // vec1 = (vec1 - vec3) / tau[1] = tr
	  // dldacf = ip - tr

	dldacf[0] /= 2;
}

*/

#endif
