/// @file GSchur.h

///////////////////////////////////////////////
// Generalized Schur Algorithm
///////////////////////////////////////////////

#ifndef GSchur_h
#define GSchur_h 1

#include "VectorFFT.h"
#include <vector>

/// Convert integer to modulo-binary representation.
///
/// Given an integer `x`, returns an integer vector \f$s = {r, s_1,
/// s_2, ..., s_k}\f$, where r < b, \f$s_i = 2^m \times b\f$ equals power of 2 times b and `sum(s) = x`.
///
/// @param[in] x Integer of which the modulo-binary represention is computed.
/// @param[in] b Integer giving the binary modulus.
/// @return Integer vector containing the modulo-binary representation.
inline std::vector<int> int2bin(int x, int b = 1) {
	std::vector<int> s;
	int x1 = x / b;
	int dif = x - x1 * b;
	int n = b;
	do {
		if (x1 & 1) {
			s.push_back(n);
		}
		n <<= 1;
	} while ((x1 >>= 1) > 0);
	if (dif) {
		s.insert(s.begin(), dif);
	}
	return s;
}

/// Memory for one iteration of Generalized Schur Algorithm. 
/// Each iteration involves polynomials \f${\alpha_{0,2n}, \beta_{0,2n}, \eta_{0,n}, \xi_{0,n}, \tilde\eta_{0,n}, \tilde\xi_{0,n}, \eta_{n,n}, \xi_{n,n}, \eta_{0,2n}, \xi_{0,2n} \f$
///
/// Note: size of input is 2n for original part and (n+m) for the merge part.
class GSchur2K {
public:

	/// VectorFFT class members for Fast Fourier Transformations.
	VectorFFT* FFT;

	double* alpha0; ///< Polynomial \f$\alpha_{0,2n}\f$.
	std::complex<double>* alpha0_fft;
	double* alphan; ///< Polynomial \f$\alpha_{n,n}\f$.
	std::complex<double>* alphan_fft;
	double* beta0; ///< Polynomial \f$\beta_{0,2n}\f$.
	std::complex<double>* beta0_fft;
	double* betan; ///< Polynomial \f$\beta_{n,n}\f$.
	std::complex<double>* betan_fft;
	double* eta0; ///< Polynomial \f$\eta_{0,n}\f$.
	std::complex<double>* eta0_fft;
	double* etat; ///< Polynomial \f$\tilde\eta_{0,n}\f$.
	std::complex<double>* etat_fft;
	double* etan; ///< Polynomial \f$\eta_{n,n}\f$.
	std::complex<double>* etan_fft;
	double* eta2n; ///< Polynomial \f$\eta_{ 0,2n }\f$.
	std::complex<double>* eta2n_fft;
	double* xi0; ///< Polynomial \f$\xi_{0,n}\f$.
	std::complex<double>* xi0_fft;
	double* xit; ///< Polynomial \f$\tilde\xi_{0,n}\f$.
	std::complex<double>* xit_fft;
	double* xin; ///< Polynomial \f$\xi_{n,n}\f$.
	std::complex<double>* xin_fft;
	double* xi2n; ///< Polynomial \f$\xi_{0,2n}\f$.
	std::complex<double>* xi2n_fft;
	double* gamma;  ///< Polynomial \f$\gamma\f$ stored in length-n vector.

	/// Constructor.
	GSchur2K(int);
	/// Destructor.
	~GSchur2K();
};

/// @param[in] n Size of the input vector.
inline GSchur2K::GSchur2K(int n) {
	FFT = new VectorFFT(n);
	alpha0 = new double[n];
	std::fill(alpha0, alpha0 + n, 0);
	alpha0_fft = new std::complex<double>[n];
	alphan = new double[n];
	alphan_fft = new std::complex<double>[n];
	beta0 = new double[n];
	std::fill(beta0, beta0 + n, 0);
	beta0_fft = new std::complex<double>[n];
	betan = new double[n];
	betan_fft = new std::complex<double>[n];
	eta0 = new double[n];
	std::fill(eta0, eta0 + n, 0);
	eta0_fft = new std::complex<double>[n];
	etat = new double[n];
	std::fill(etat, etat + n, 0);
	etat_fft = new std::complex<double>[n];
	etan = new double[n];
	etan_fft = new std::complex<double>[n];
	eta2n = new double[n];
	eta2n_fft = new std::complex<double>[n];
	xi0 = new double[n];
	std::fill(xi0, xi0 + n, 0);
	xi0_fft = new std::complex<double>[n];
	xit = new double[n];
	std::fill(xit, xit + n, 0);
	xit_fft = new std::complex<double>[n];
	xin = new double[n];
	xin_fft = new std::complex<double>[n];
	xi2n = new double[n];
	xi2n_fft = new std::complex<double>[n];
	gamma = new double[n];
}

inline GSchur2K::~GSchur2K() {
	delete FFT;

	delete[] alpha0;
	delete[] alpha0_fft;
	delete[] alphan;
	delete[] alphan_fft;
	delete[] beta0;
	delete[] beta0_fft;
	delete[] betan;
	delete[] betan_fft;
	delete[] eta0;
	delete[] eta0_fft;
	delete[] etat;
	delete[] etat_fft;
	delete[] etan;
	delete[] etan_fft;
	delete[] eta2n;
	delete[] eta2n_fft;
	delete[] xi0;
	delete[] xi0_fft;
	delete[] xit;
	delete[] xit_fft;
	delete[] xin;
	delete[] xin_fft;
	delete[] xi2n;
	delete[] xi2n_fft;
	delete[] gamma;
}

/// Generalized Schur Algorithm for arbitrary size N.
///
/// The implementation of Generalized Schur algorithm for any `N` contains
/// following steps:
///
/// 1) Decompose integer `N` into the binary form w.r.t modulus `b` and obtain a
/// vector `s` s.t. \f$s = {r, s_1, s_2, ..., s_k}\f$. This step is done with function `int2bin`. As in GSchur
/// Algorithm, whole size-N transformation \f$T_{0,N}\f$ can be decomposed into
/// binary pieces 
///
/// 2) Compute the power-2-size transformations \f$T_{..., s_i}\f$ using function `GenStep`. 
///
/// 3) Merge the binary pieces in ascending order in following way with function `GenMerge`.
/// \f[
///   T_{0,N} = T_{0, r} \circ T_{r, s_1} \circ ... \circ T_{..., s_m}
/// \f]
class GSchurN {
	int N;  ///< Size of input vector
	int b;  ///< Integer giving binary modulus. Binary pieces of size smaller than
			///< `b` are computed using Progressive Schur Algorithm.
	double* alpha;       ///< Input polynomial in form of length-`N` vector.
	double* beta;        ///< Input polynomial in form of length-`N` vector.
	std::vector<int> s;  ///< Vector that records the binary decomposition of `N`.
						 ///< \f$s = {2^{k_0} \times b, 2^{k_1} \times b, ...,
						 ///< 2^{k_T} \times b, r}\f$.
	int k; ///< Length of s.
	GSchur2K**
		gs;  ///< Vector of `GSchur2K` for the binary pieces \f$ T_{0, s[1]}
			 ///< ,T_{s[1], s[2]},...\f$. Because of the doubling procedure, the
			 ///< computation of \f$ T_{, 2^k\times b}\f$ requires `GSchur2K` of
			 ///< size \f$(b, 2\times b, ..., 2^k\times b)\f$. As a result the
			 ///< length of vector `gs` is determined by the largest element (also
			 ///< the first element) of vector `s`.
	GSchur2K** gsM;  ///< Vector of GSchur2K for merging the binary pieces \f$
					 ///< T_{0, s[1]} ,T_{s[1], s[2]},...\f$. The length of vector
					 ///< `gsM` is determined by the length of vector `s`.

	/// Progress function that computes \f$ \{\alpha_{n,n} \beta_{n,n} \} \f$ from \f$ \alpha_{0,2n}, \beta_{0,2n},
	/// \eta_{0,n}, \xi_{0,n} \f$.
	void alpha2Beta(GSchur2K*, int, int);
	/// Progress function that computes \f$ \{\eta_{0,2n} \xi_{0,2n} \} \f$ from \f$ \eta_{0,n}, \xi_{0,n},
	/// \eta_{n,n}, \xi_{n,n} \f$.
	void eta2Xi(GSchur2K*, int);

	/// Progressive Schur Algorithm for arbitrary size.
	void ProgStep(double*, double*, int);

	/// Compute the binary pieces \f$T_{..,s_i}\f$.
	void GenStep(double*, double*, int, int);

	/// Merge the binary pieces \f$T_{..,s_i}\f$ to get \f$T_{0,N}\f$
	void GenMerge();

public:
	double* delta;  ///< Real output vector. The first column of inverse Toeplitz matrix.
	double ldV;   ///< Real number. Log-determinant of Toeplitz matrix.

	/// Constructor.
	GSchurN(int);
	/// Destructor.
	~GSchurN();
	/// Perform the Generalized Schur algorithm on the input data.
	void Compute(double*);
};

/// @param[in] N_ Size of Toeplitz matrix (Size of Generalized Schur Algorithm
/// is `N_-1`).
/// @param[in] b_ Integer for binary modulus.
inline GSchurN::GSchurN(int N_) {
	N = N_;
	b = 1;
	alpha = new double[N - 1];
	beta = new double[N - 1];
	delta = new double[N];

	/// s = [r, s_1, s_2, ..., s_k]
	s = int2bin(N - 1, b);
	k = s.size();
	int gs_size = b;
	int gs_layer = log2(ceil((double)s[k - 1] / b));
	gs = new GSchur2K * [gs_layer + 1];
	gs[0] = new GSchur2K(2 * gs_size);
	for (int ii = 0; ii < gs_layer; ++ii) {
		gs[ii + 1] = new GSchur2K(2 * gs_size);
		gs_size <<= 1;
	}

	// Length of gsM is max(k, 2) - 1
	gsM = new GSchur2K * [k > 2 ? k : 2 - 1];
	if (k == 1) {
		gsM[0] = new GSchur2K(s[0]);
	}
	else {
		int gsM_size = s[0];
		for (int ii = 0; ii < k - 1; ++ii) {
			gsM_size = gsM_size + s[ii + 1];
			gsM[ii] = new GSchur2K(gsM_size);
		}
	}
}

inline GSchurN::~GSchurN() {
	delete[] alpha;
	delete[] beta;
	delete[] delta;

	for (int ii = 0; ii <= log2(ceil((double)s[k - 1] / b)); ++ii) {
		delete gs[ii];
	}
	delete[] gs;

	if (k == 1) {
		delete gsM[0];
	}
	else {
		for (int ii = 0; ii < k - 1; ++ii) {
			delete gsM[ii];
		}
	}
	delete[] gsM;
}

/// For given length-N polynomials \f$\alpha_{0,2n}, \beta_{0,2n},
/// \eta_{0,n}(x), \xi_{0,n}(x)\f$, compute polynomials \f$
/// \alpha_{n,n}, \beta_{n,n} \f$ defined as \f[
///   \alpha_{n,n} = (\alpha_{0,2n} \eta_{0,n} - \beta_{0,2n} \xi_{0,n}) / x^n
///	  \beta_n = (\beta_{0,2n} \tilde\eta_{0,n} - \alpha_{0,2n} \tilde\xi_{0,n}) / x^n
/// \f]
/// where \f$ \tilde\eta_n = x^n \eta_n(1/x) \f$
///
/// Note: by negating the odd parts of the
/// conjugation of fft(\f$\eta_{0,n}, 0, 0, ...\f$), we can directly get
/// fft(\f$0, rev(\eta_{0,n}), 0, ...\f$). But it only works when length of zero-padding is also n.
///
/// @param[in] gsr GSchur2K class.
/// @param[in] n1 Length of \f$\eta_{0,n}, \xi_{0,n}\f$
/// @param[in] n2 Length of zero-padding. size of FFT is n1 + n2.
/// @return Output complex vectors for \f$\alpha_{n,n} and \beta_{n,n} \f$ stored in `gsr->betam_IFF->out` and `gsr->betam_IFF->out`
inline void GSchurN::alpha2Beta(GSchur2K* gsr, int n1, int n2) {

	gsr->FFT->fft(gsr->alpha0_fft, gsr->alpha0);
	gsr->FFT->fft(gsr->beta0_fft, gsr->beta0);
	gsr->FFT->fft(gsr->eta0_fft, gsr->eta0);
	gsr->FFT->fft(gsr->xi0_fft, gsr->xi0);

	if (n1 - n2) {
		for (int ii = 0; ii < n1; ++ii) {
			gsr->etat[n1 - ii] = gsr->eta0[ii];
			gsr->xit[n1 - ii] = gsr->xi0[ii];
		}

		gsr->FFT->fft(gsr->etat_fft, gsr->etat);
		gsr->FFT->fft(gsr->xit_fft, gsr->xit);

	}
	else {
		for (int ii = 0; ii < n1; ++ii) {

			gsr->etat_fft[2 * ii] = std::conj(gsr->eta0_fft[2 * ii]);
			gsr->etat_fft[2 * ii + 1] = -std::conj(gsr->eta0_fft[2 * ii + 1]);
			gsr->xit_fft[2 * ii] = std::conj(gsr->xi0_fft[2 * ii]);
			gsr->xit_fft[2 * ii + 1] = -std::conj(gsr->xi0_fft[2 * ii + 1]);

		}
	}

	int nn = (n1 + n2) / 2 + 1;

	Complex_Mult(gsr->alphan_fft, gsr->alpha0_fft, gsr->eta0_fft, nn);
	Complex_Mult_Minus(gsr->alphan_fft, gsr->xi0_fft, gsr->beta0_fft, nn);
	Complex_Mult(gsr->betan_fft, gsr->beta0_fft, gsr->etat_fft, nn);
	Complex_Mult_Minus(gsr->betan_fft, gsr->xit_fft, gsr->alpha0_fft, nn);

	gsr->FFT->ifft(gsr->alphan, gsr->alphan_fft);
	gsr->FFT->ifft(gsr->betan, gsr->betan_fft);
}

/// For given length-N polynomials \f$\eta_{n,n}, \xi_{n,n}, \eta_{0,n},
/// \xi_{0,n}\f$, compute polynomials \f$ \eta_{2n}, \xi_{2n} \f$
/// defined as
/// \f[
///   \xi_{2n} = \tilde\eta_{0,n} \xi_{n,n} + \xi_n \eta_{n,n}, 
///   \eta_{2n} = \tilde\xi_{0,n} \xi_{n,n} + \eta_n \eta_{n,n}
/// \f]
///
/// @param[in] gsr GSchur2K class.
/// @param[in] n Length of inputs (integer). Size of gsr members is n.
/// @return Output complex vectors for \f$ \eta_{2n}, \xi_{2n} \f$ stored in `gsr->xi2n` and `gsr->eta2n`.
inline void GSchurN::eta2Xi(GSchur2K* gsr, int n) {
	gsr->FFT->fft(gsr->xin_fft, gsr->xin);
	gsr->FFT->fft(gsr->etan_fft, gsr->etan);

	int n2 = n / 2 + 1;

	Complex_Mult(gsr->xi2n_fft, gsr->etat_fft, gsr->xin_fft, n2);
	Complex_Mult_Plus(gsr->xi2n_fft, gsr->xi0_fft, gsr->etan_fft, n2);

	Complex_Mult(gsr->eta2n_fft, gsr->xit_fft, gsr->xin_fft, n2);
	Complex_Mult_Plus(gsr->eta2n_fft, gsr->eta0_fft, gsr->etan_fft, n2);

	gsr->FFT->ifft(gsr->xi2n, gsr->xi2n_fft);
	gsr->FFT->ifft(gsr->eta2n, gsr->eta2n_fft);

}

/// It is an \f$O(n^2)\f$ algorithm . Numerical experiments shows that for
/// \f$N\leq 64\f$ PSchur Algorithm is faster than GSchur algorithm. Thus we
/// choose binary modulus to be `64` in implementation, i.e. computation of
/// pieces smaller than `64` is done with PSchur and the merging of such pieces
/// are done with GSchur.
///
/// @param[in] alpha0 Input real vector.
/// @param[in] beta0 Input real vector.
/// @param[in] n Size of inputs (integer).
inline void GSchurN::ProgStep(double* alpha0, double* beta0, int n) {
	double* tmpPtr, * xi1, * xi2, * eta1, * eta2, alpha_0, alpha_1;

	// Initialize the memory
	std::fill(gs[0]->xin, gs[0]->xin + 2 * n, 0);
	std::fill(gs[0]->etan, gs[0]->etan + 2 * n, 0);

	// Take advantage of the memory of existing `GSchur2K` class
	eta1 = gs[0]->etan;
	eta2 = gs[0]->etan + n;
	xi1 = gs[0]->xin;
	xi2 = gs[0]->xin + n;

	eta1[0] = 1.0;
	xi1[0] = alpha0[0] / beta0[0];
	gs[0]->gamma[0] = xi1[0];
	gs[0]->beta0[0] = beta0[0] * (1 - xi1[0] * xi1[0]);

	for (int kk = 1; kk < n; ++kk) {
		alpha_1 = alpha0[kk];
		gs[0]->beta0[kk] = beta0[kk];

		for (int jj = 1; jj <= kk; ++jj) {
			alpha_0 = alpha_1 - gs[0]->gamma[jj - 1] * gs[0]->beta0[kk - jj + 1];
			gs[0]->beta0[kk - jj + 1] -= gs[0]->gamma[jj - 1] * alpha_1;
			alpha_1 = alpha_0;
		}

		gs[0]->gamma[kk] = alpha_0 / gs[0]->beta0[0];
		gs[0]->beta0[0] *= 1 - gs[0]->gamma[kk] * gs[0]->gamma[kk];

		eta2[0] = 1.0;
		xi2[0] = alpha0[0] / beta0[0];
		for (int jj = 1; jj <= kk; ++jj) {
			xi2[jj] = xi1[jj] + gs[0]->gamma[kk] * eta1[kk - jj];
			eta2[jj] = eta1[jj] + gs[0]->gamma[kk] * xi1[kk - jj];
		}
		tmpPtr = xi1;
		xi1 = xi2;
		xi2 = tmpPtr;
		tmpPtr = eta1;
		eta1 = eta2;
		eta2 = tmpPtr;
	}
	
	std::copy(xi1, xi1 + n, gs[0]->xi2n);
	std::copy(eta1, eta1 + n, gs[0]->eta2n);
	return;
}

/// This function implements the recursive doubling procedure in Generalized
/// Schur Algorithm. For an input size `si`, function use memory gs[0], ..., gs[layer] to compute
/// polynomials \f$ \eta_{0,s[i]}, xi_{0,s[i]} and \gamma_{s[i]} \f$ from \f$
/// \alpha_{0,s[i]}, \beta_{0,s[i]}\f$. Size `si` is an element of vector `s` and is either
/// \f$\leq b\f$ (where layer = 0) or \f$= 2^{layer} \times b\f$.
///
/// @param[in] alpha0 First input real vector.
/// @param[in] beta0 Second input real vector.
/// @param[in] si Length of inputs (integer).
/// @param[in] layer Index of address for computation, i.e. computation happens in gs[0,1,..., layer]
/// @return polynomials \f$ \{ \eta_{0,s_i}, \xi_{0,s_i}, \gamma_{0,s_i}\} \f$
/// stored in `gs[layer]->eta2n`, `gs[layer]->xi2n`,
/// `gs[layer]->gamma`.
inline void GSchurN::GenStep(double* alpha0, double* beta0, int si, int layer) {
	// When s[i] equals r.

	if (si <= b) {
		ProgStep(alpha0, beta0, si);
		return;
	}

	// When s[i] equals 2^layer * b.
	ProgStep(alpha0, beta0, b);
	int n = b;
	for (int m = 0; m < layer; ++m) {
		// Prepare the inputs for function alpha2Beta, including alpha-0n, beta-0n,
		// xi-0n, eta-0n and gamma-0n

		std::copy(alpha0, alpha0 + 2 * n, gs[m + 1]->alpha0);
		std::copy(beta0, beta0 + 2 * n, gs[m + 1]->beta0);
		std::copy(gs[m]->eta2n, gs[m]->eta2n + n, gs[m + 1]->eta0);
		std::copy(gs[m]->xi2n, gs[m]->xi2n + n, gs[m + 1]->xi0);
		std::copy(gs[m]->gamma, gs[m]->gamma + n, gs[m + 1]->gamma);


		// Compute alpha-nn and beta-nn
		alpha2Beta(gs[m + 1], n, n);
		// Compute eta-nn, xi-nn and gamma-nn using function GenStep
		GenStep(gs[m + 1]->alphan + n, gs[m + 1]->betan + n, n, m);

		// Prepare the inputs for function eta2Xi, including xi-nn and eta-nn
		std::copy(gs[m]->xi2n, gs[m]->xi2n + n, gs[m + 1]->xin);
		std::fill(gs[m + 1]->xin + n, gs[m + 1]->xin + 2 * n, 0);

		std::copy(gs[m]->eta2n, gs[m]->eta2n + n, gs[m + 1]->etan);
		std::fill(gs[m + 1]->etan + n, gs[m + 1]->etan + 2 * n, 0);

		// Compute eta-02n and xi-02n
		eta2Xi(gs[m + 1], 2 * n);
		// Paste the computed gamma-nn after gamma-0n to get gamma-02n
		std::copy(gs[m]->gamma, gs[m]->gamma + n, gs[m + 1]->gamma + n);
		n <<= 1;
	}
	return;
}

/// @return polynomials \f$ \{ \eta_{0,N}, \xi_{0,N}, \gamma_{0,N}\} \f$ stored
/// in `gsM[0]->eta2n`, `gsM[0]->xi2n`, `gsM[0]->gamma`.
inline void GSchurN::GenMerge() {
	int layer = log2(ceil((double)s[0] / b));

	GenStep(alpha, beta, s[0], layer);

	// When vector `s` has only one element, directly pass the results to output.
	if (k == 1) {
		std::copy(gs[layer]->eta2n, gs[layer]->eta2n + s[0], gsM[0]->eta2n);
		std::copy(gs[layer]->xi2n, gs[layer]->xi2n + s[0], gsM[0]->xi2n);
		std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
		return;
	}

	int n = 0;

	for (int m = 0; m < k - 1; ++m) {
		n += s[m];

		if (m) {
			std::copy(gsM[m - 1]->eta2n, gsM[m - 1]->eta2n + n, gsM[m]->eta0);
			std::copy(gsM[m - 1]->xi2n, gsM[m - 1]->xi2n + n, gsM[m]->xi0);
		}
		else {
			
			std::copy(gs[layer]->eta2n, gs[layer]->eta2n + n,
				gsM[m]->eta0);
			std::copy(gs[layer]->xi2n, gs[layer]->xi2n + n,
				gsM[m]->xi0);
			std::copy(gs[layer]->gamma, gs[layer]->gamma + n, gsM[k - 2]->gamma);
		}

		// First paste alpha_{0,n} and beta_{0,n} into gsM[0]->alpha and
		// gsM[0]->beta, length is n
		
		std::copy(alpha, alpha + n + s[m + 1], gsM[m]->alpha0);
		std::copy(beta, beta + n + s[m + 1], gsM[m]->beta0);


		// alpha2Beta uses gsM[m]->alpha0, beta0, eta0 and xi0
		// to compute alphan and betan
		alpha2Beta(gsM[m], n, s[m + 1]);

		// Computation of GenStep happens inside gs[0,1,...,layer].
		// We just need alphan [s[m]+(1:s[m+1])]
		layer = log2(ceil((double)s[m + 1] / b));
		
		GenStep(gsM[m]->alphan + n, gsM[m]->betan + n,
			s[m + 1], layer);	  
		
		std::copy(gs[layer]->eta2n, gs[layer]->eta2n + s[m + 1], gsM[m]->etan);
		std::fill(gsM[m]->etan + s[m + 1], gsM[m]->etan + n + s[m + 1], 0);
		std::copy(gs[layer]->xi2n, gs[layer]->xi2n + s[m + 1], gsM[m]->xin);
		std::fill(gsM[m]->xin + s[m + 1], gsM[m]->xin + n + s[m + 1], 0);
		std::copy(gs[layer]->gamma, gs[layer]->gamma + s[m + 1], gsM[k - 2]->gamma + n);


		// uses xi_t_fft, xi_0_fft, xi_n_fft (eta), returns xi_0_ifft
		eta2Xi(gsM[m], s[m + 1] + n);
	}
}

/// This function generates the input polynomials \f$\alpha_0, \beta_0\f$ of
/// Generalized Schur algorithm from the first column of Toeplitz matrix, run
/// the Merged GSchur algorithm and finally convert the output of GSchur
/// algorithm into the log-determinant and first column of inverse Toeplitz
/// matrix
///
/// @param[in] acf Input real vector, the first column of Toeplitz matrix.
inline void GSchurN::Compute(double* acf) {

	// generate alpha_0 and beta_0 from acf.
	for (int ii = 0; ii < N - 1; ++ii) {
		alpha[ii] = -1 * acf[ii + 1];
		beta[ii] = acf[ii];
	}

	// run gschur algorithm
	GenMerge();

	// convert gamma into ldv.
	double sigma2 = log(acf[0]);
	ldV = sigma2;
	int k2 = k > 2 ? k : 2;
	for (int ii = 0; ii < N - 1; ++ii) {
		if (gsM[k2 - 2]->gamma[ii] < 1) {
			sigma2 += log(1 - gsM[k2 - 2]->gamma[ii] * gsM[k2 - 2]->gamma[ii]);
			ldV += sigma2;
		}
	}
	sigma2 = exp(sigma2);

	// convert eta and xi into delta, the first column of inverse Toeplitz matrix.
	std::copy(gsM[k2 - 2]->eta2n, gsM[k2 - 2]->eta2n + N - 1, delta);
	delta[N - 1] = 0;
	delta[0] /= sigma2;
	for (int ii = 1; ii < N; ++ii) {
		delta[ii] += gsM[k2 - 2]->xi2n[ii - 1];
		delta[ii] /= sigma2;
	}
	return;
}

#endif
