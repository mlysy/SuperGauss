/// @file NormalToeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef NormalToeplitz_h
#define NormalToeplitz_h 1

#include "Toeplitz.h"
# define PI 3.141592653589793238462643383279502884L // PI value

/// Class for Computation involving Toeplitz matrix.
///
/// Model: z ~ N_(0, V), where V is a Toeplitz matrix with parameter theta
class NormalToeplitz {
	Toeplitz* Tz_;
	Toeplitz* T2_;
	double* tempVec_;
	double* tempVec0_;
	double* tempVec1_;
	double* tempVec2_;
	bool hT_;
public:
	/// Constructor.
	/// Overload constructor doesnot work, currently add a FLAG hasToep
	NormalToeplitz(const int N, const int p, bool hasToep);

	/// Destructor
	~NormalToeplitz();

	int N_; // length of observation z
	int p_; // number of parameters theta

	/// Log-Density.
	double logdens(double* z, double* acf);
	/// this is the preferred form for R, 
	/// because we don't want to reallocate memory for every call.
	/// note that `Tz` can't be `const` because we potentially modify its internal structure.
	double logdens(double* z, Toeplitz* Tz);

	/// Gradient with respect to theta.
	/// The output `dldt` is the same length as `theta`.
	/// dzdt and dacfdz are p by N matrices.
	void grad(double* dldt,
		double* z, double** dzdt,
		double* acf, double** dacfdt);
	/// same thing with Toeplitz input
	void grad(double* dldt,
		double* z, double** dzdt,
		Toeplitz* Tz, double** dacfdt);

	/// Hessian with respect to theta.
	void hess(double** d2ldt,
		double* z, double** dzdt, double*** d2zdt,
		double* acf, double** dacfdt, double*** d2acfdt);
	/// same thing with Toeplitz input
	void hess(double** d2ldt,
		double* z, double** dzdt, double*** d2zdt,
		Toeplitz* Tz, double** dacfdt, double*** d2acfdt);

	/// Full gradient for auto-diff.
	/// The outputs `dldz` and `dldacf` are each the length of `z`.
	void grad_full(double* dldz, double* dldacf,
		double* z, double* acf);
	/// same thing with Toeplitz input
	void grad_full(double* dldz, double* dldacf,
		double* z, Toeplitz* Tz);
};

/// @param[in] N Size of Toeplitz matrix.
inline NormalToeplitz::NormalToeplitz(int N, int p, bool hasToep) {
	if (!hasToep) {
		Tz_ = new Toeplitz(N, 1);
	}
	hT_ = hasToep; // flag for destructor
	T2_ = new Toeplitz(N, 1);
	tempVec_ = new double[N];
	tempVec0_ = new double[N];
	tempVec1_ = new double[N];
	tempVec2_ = new double[N];
	N_ = N;
	p_ = p;
}

/// Overloaded Constructor of NormalToeplitz class, using internal pointer to existing Toeplitz matrix
/// The argument cannot be void, since we still need the N and p to create necessary memory
/*
inline NormalToeplitz::NormalToeplitz(int N, int p) {
	T2_ = new Toeplitz(N, 1);
	tempVec_ = new double[N];
	tempVec0_ = new double[N];
	tempVec1_ = new double[N];
	tempVec2_ = new double[N];
	N_ = N;
	p_ = p;
}
*/

inline NormalToeplitz::~NormalToeplitz() {
	if (!hT_) {
		delete Tz_;
	}
	delete T2_;
	delete tempVec_;
	delete tempVec0_;
	delete tempVec1_;
	delete tempVec2_;
}

/// Function for vector multiplication
///
/// @param[in] v1 First input vector.
/// @param[in] v2 Second input vector.
/// @param[in] n Length of input vector.
/// @param[out] double number that equals v1 * v2.
inline double vectorMult(double* v1, double* v2, int n) {
	double ans = 0;
	for (int ii = 0; ii < n; ++ii) {
		ans += v1[ii] * v2[ii];
	}
	return ans;
}

inline double NormalToeplitz::logdens(double* z, double* acf) {
	double ldens = 0;
	Tz_->setAcf(acf); // Tz_ = Toeplitz(acf)
	Tz_->solveVec(tempVec_, z); // tempVec_ = Tz_^{-1} * z
	ldens = vectorMult(z, tempVec_, N_); // ldens = t(z) * Tz_^{-1} * z
	ldens += Tz_->logDet() + N_ * log(PI);
	ldens *= -0.5;
	return ldens;
}

inline double NormalToeplitz::logdens(double* z, Toeplitz* Tz) {
	double ldens = 0;
	// acf is already stored in Tz
	Tz->solveVec(tempVec_, z); // tempVec_ = Tz^{-1} * z
	ldens = vectorMult(z, tempVec_, N_); // ldens = t(z) * Tz^{-1} * z
	ldens += Tz->logDet() + N_ * log(PI);
	ldens *= -0.5;
	return ldens;
}

inline void NormalToeplitz::grad(double* dldt, double* z, double** dzdt,
	double* acf, double** dacfdt) {
	Tz_->setAcf(acf); // Tz_ = Toeplitz(acf)
	Tz_->solveVec(tempVec_, z); // tempVec_ = Tz_^{-1} * z
	double tmp;
	for (int ii = 0; ii < p_; ++ii) {
		T2_->setAcf(dacfdt[ii]);
		T2_->multVec(tempVec0_, tempVec_);
		dldt[ii] = vectorMult(tempVec_, tempVec0_, N_) / 2;
		dldt[ii] -= vectorMult(dzdt[ii], tempVec_, N_);
		dldt[ii] -= Tz_->traceProd(dacfdt[ii]) / 2;
	}
}

inline void NormalToeplitz::grad(double* dldt, double* z, double** dzdt,
	Toeplitz* Tz, double** dacfdt) {
	// acf is already stored in Tz
	Tz->solveVec(tempVec_, z); // tempVec_ = Tz^{-1} * z
	double tmp;
	for (int ii = 0; ii < p_; ++ii) {
		T2_->setAcf(dacfdt[ii]);
		T2_->multVec(tempVec0_, tempVec_);
		dldt[ii] = vectorMult(tempVec_, tempVec0_, N_) / 2;
		dldt[ii] -= vectorMult(dzdt[ii], tempVec_, N_);
		dldt[ii] -= Tz->traceProd(dacfdt[ii]) / 2;
	}
}

inline void NormalToeplitz::hess(double** d2ldt,
	double* z, double** dzdt, double*** d2zdt,
	double* acf, double** dacfdt, double*** d2acfdt) {
	Tz_->setAcf(acf);
	Tz_->solveVec(tempVec_, z);
	double ans;
	for (int ii = 0; ii < p_; ++ii) {
		for (int jj = 0; jj <= ii; ++jj) {
			T2_->setAcf(dacfdt[jj]);
			T2_->multVec(tempVec2_, tempVec_);
			T2_->setAcf(dacfdt[ii]);
			T2_->multVec(tempVec1_, tempVec_);

			ans = vectorMult(d2zdt[ii][jj], tempVec_, N_);
			// temp0 = solve(Tz_, temp2)
			Tz_->solveVec(tempVec0_, tempVec2_);
			ans -= vectorMult(dzdt[ii], tempVec0_, N_);
			ans += vectorMult(tempVec1_, tempVec0_, N_);
			// temp0 = solve(Tz_, temp1)
			Tz_->solveVec(tempVec0_, tempVec1_);
			ans -= vectorMult(dzdt[jj], tempVec0_, N_);
			// temp0 = solve(Tz_, dZ[,jj])
			Tz_->solveVec(tempVec0_, dzdt[jj]);
			ans += vectorMult(dzdt[ii], tempVec0_, N_);
			ans *= 2;

			T2_->setAcf(d2acfdt[ii][jj]);
			T2_->multVec(tempVec0_, tempVec_);
			ans -= vectorMult(tempVec_, tempVec0_, N_);
			ans += Tz_->traceProd(d2acfdt[ii][jj]);
			ans -= Tz_->traceDeriv(dacfdt[ii], dacfdt[jj]);

			d2ldt[ii][jj] = -ans / 2;
		}
	}

	if (p_ > 1) {
		for (int ii = 0; ii < p_; ++ii) {
			for (int jj = ii + 1; jj < p_; ++jj) {
				d2ldt[ii][jj] = d2ldt[jj][ii];
			}
		}
	}

}

inline void NormalToeplitz::hess(double** d2ldt,
	double* z, double** dzdt, double*** d2zdt,
	Toeplitz* Tz, double** dacfdt, double*** d2acfdt) {
	// acf is already stored in Tz
	Tz->solveVec(tempVec_, z);
	double ans;
	for (int ii = 0; ii < p_; ++ii) {
		for (int jj = 0; jj <= ii; ++jj) {
			T2_->setAcf(dacfdt[jj]);
			T2_->multVec(tempVec2_, tempVec_);
			T2_->setAcf(dacfdt[ii]);
			T2_->multVec(tempVec1_, tempVec_);

			ans = vectorMult(d2zdt[ii][jj], tempVec_, N_);
			// temp0 = solve(Tz, temp2)
			Tz->solveVec(tempVec0_, tempVec2_);
			ans -= vectorMult(dzdt[ii], tempVec0_, N_);
			ans += vectorMult(tempVec1_, tempVec0_, N_);
			// temp0 = solve(Tz, temp1)
			Tz->solveVec(tempVec0_, tempVec1_);
			ans -= vectorMult(dzdt[jj], tempVec0_, N_);
			// temp0 = solve(Tz, dZ[,jj])
			Tz->solveVec(tempVec0_, dzdt[jj]);
			ans += vectorMult(dzdt[ii], tempVec0_, N_);
			ans *= 2;

			T2_->setAcf(d2acfdt[ii][jj]);
			T2_->multVec(tempVec0_, tempVec_);
			ans -= vectorMult(tempVec_, tempVec0_, N_);
			ans += Tz->traceProd(d2acfdt[ii][jj]);
			ans -= Tz->traceDeriv(dacfdt[ii], dacfdt[jj]);

			d2ldt[ii][jj] = -ans / 2;
		}
	}

	if (p_ > 1) {
		for (int ii = 0; ii < p_; ++ii) {
			for (int jj = ii + 1; jj < p_; ++jj) {
				d2ldt[ii][jj] = d2ldt[jj][ii];
			}
		}
	}

}

inline void NormalToeplitz::grad_full(double* dldz, double* dldacf,
	double* z, double* acf) {
	Tz_->setAcf(acf);

	// gradient with respect to z
	Tz_->solveVec(dldz, z);
	for (int ii = 0; ii < N_; ii++) {
		dldz[ii] = -dldz[ii];
	}

	// gradient with respect to acf
	Tz_->solveVec(tempVec_, z); // tempVec_ = Vz

	tempVec0_[0] = 1;
	std::fill(tempVec0_ + 1, tempVec0_ + N_, 0); // tempVec0_ = [1,0,0,...,0]
	Tz_->solveVec(tempVec1_, tempVec0_); // tempVec1_ = tau
	double tau1 = tempVec1_[0];

	T2_->setAcf(tempVec_);
	T2_->mult0Vec(dldacf, tempVec_); // dldacf = upper.toep(Vz) %*% Vz = ip

	tempVec2_[0] = 0;
	for (int ii = 1; ii < N_; ++ii) {
		tempVec2_[ii] = tempVec1_[N_ - ii];
	} // tempVec2_ = tau2

	for (int ii = 0; ii < N_; ++ii) {
		tempVec0_[ii] = (N_ - ii) * tempVec1_[ii];
	} // tempVec0_ = (N_:1 * tau)

	T2_->setAcf(tempVec1_);
	T2_->mult0Vec(tempVec_, tempVec0_); // tempVec_ = upper.toep(tau) %*% (N_:1 * tau) = tr

	for (int ii = 0; ii < N_; ++ii) {
		tempVec0_[ii] = (N_ - ii) * tempVec2_[ii];
	} // tempVec0_ = (N_:1 * tau2)
	T2_->setAcf(tempVec2_);
	T2_->mult0Vec(tempVec1_, tempVec0_); // tempVec1_ = upper.toep(tau2) %*% (N_:1 * tau2)

	for (int ii = 0; ii < N_; ++ii) {
		tempVec_[ii] -= tempVec1_[ii];
		tempVec_[ii] /= tau1;
		dldacf[ii] -= tempVec_[ii];
	} // tempVec_ = (tempVec_ - tempVec1_) / tau[1] = tr
	  // dldacf = ip - tr

	dldacf[0] /= 2;
}

inline void NormalToeplitz::grad_full(double* dldz, double* dldacf,
	double* z, Toeplitz* Tz) {
	// acf is already stored in Tz

	// gradient with respect to z
	Tz->solveVec(dldz, z);
	for (int ii = 0; ii < N_; ii++) {
		dldz[ii] = -dldz[ii];
	}

	// gradient with respect to acf
	Tz->solveVec(tempVec_, z); // tempVec_ = Vz

	tempVec0_[0] = 1;
	std::fill(tempVec0_ + 1, tempVec0_ + N_, 0); // tempVec0_ = [1,0,0,...,0]
	Tz->solveVec(tempVec1_, tempVec0_); // tempVec1_ = tau
	double tau1 = tempVec1_[0];

	T2_->setAcf(tempVec_);
	T2_->mult0Vec(dldacf, tempVec_); // dldacf = upper.toep(Vz) %*% Vz = ip

	tempVec2_[0] = 0;
	for (int ii = 1; ii < N_; ++ii) {
		tempVec2_[ii] = tempVec1_[N_ - ii];
	} // tempVec2_ = tau2

	for (int ii = 0; ii < N_; ++ii) {
		tempVec0_[ii] = (N_ - ii) * tempVec1_[ii];
	} // tempVec0_ = (N_:1 * tau)

	T2_->setAcf(tempVec1_);
	T2_->mult0Vec(tempVec_, tempVec0_); // tempVec_ = upper.toep(tau) %*% (N_:1 * tau) = tr

	for (int ii = 0; ii < N_; ++ii) {
		tempVec0_[ii] = (N_ - ii) * tempVec2_[ii];
	} // tempVec0_ = (N_:1 * tau2)
	T2_->setAcf(tempVec2_);
	T2_->mult0Vec(tempVec1_, tempVec0_); // tempVec1_ = upper.toep(tau2) %*% (N_:1 * tau2)

	for (int ii = 0; ii < N_; ++ii) {
		tempVec_[ii] -= tempVec1_[ii];
		tempVec_[ii] /= tau1;
		dldacf[ii] -= tempVec_[ii];
	} // tempVec_ = (tempVec_ - tempVec1_) / tau[1] = tr
	  // dldacf = ip - tr

	dldacf[0] /= 2;
}

#endif
