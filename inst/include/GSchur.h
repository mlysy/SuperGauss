/// @file GSchur.h

///////////////////////////////////////////////
// Generalized Schur Algorithm
///////////////////////////////////////////////

#ifndef GSchur_h
#define GSchur_h 1

#include "VectorFFT.h"

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
  VectorFFT* alpha_0_FFT;  ///< Polynomial \f$\alpha_{0,2n}\f$.
  VectorFFT* beta_0_FFT;   ///< Polynomial \f$\beta_{0,2n}\f$.
  VectorFFT* eta_0_FFT;  ///< Polynomial \f$\eta_{0,n}\f$.
  VectorFFT* eta_t_FFT;  ///< Polynomial \f$\tilde\eta_{0,n}\f$.
  VectorFFT* xi_0_FFT;  ///< Polynomial \f$\xi_{0,n}\f$.
  VectorFFT* xi_t_FFT;  ///< Polynomial \f$\tilde\xi_{0,n}\f$.
  VectorFFT* xi_n_FFT;  ///< Polynomial \f$\xi_{n,n}\f$.
  VectorFFT* eta_n_FFT;  ///< Polynomial \f$\eta_{n,n}\f$.
  VectorIFFT* alpha_n_IFFT;  ///< Polynomial \f$\alpha_{n,n}\f$.
  VectorIFFT* beta_n_IFFT;   ///< Polynomial \f$\beta_{n,n}\f$.
  VectorIFFT* xi_0_IFFT;  ///< Polynomial \f$\xi_{0,2n}\f$.
  VectorIFFT* eta_0_IFFT;  ///< Polynomial \f$\eta_{0,2n}\f$.
  double* gamma;  ///< Polynomial \f$\gamma\f$ stored in length-n vector.
  
  /// Constructor.
  GSchur2K(int);
  /// Destructor.
  ~GSchur2K();
};

/// @param[in] n Size of the input vector.
inline GSchur2K::GSchur2K(int n) {
  alpha_0_FFT = new VectorFFT(n);
  beta_0_FFT = new VectorFFT(n);
  eta_0_FFT = new VectorFFT(n);
  xi_0_FFT = new VectorFFT(n);
  eta_t_FFT = new VectorFFT(n);
  xi_t_FFT = new VectorFFT(n);
  xi_n_FFT = new VectorFFT(n);
  eta_n_FFT = new VectorFFT(n);
  alpha_n_IFFT = new VectorIFFT(n);
  beta_n_IFFT = new VectorIFFT(n);
  xi_0_IFFT = new VectorIFFT(n);
  eta_0_IFFT = new VectorIFFT(n);
  gamma = new double[n];
}

inline GSchur2K::~GSchur2K() {
  delete alpha_0_FFT;
  delete beta_0_FFT;
  delete eta_0_FFT;
  delete xi_0_FFT;
  delete eta_t_FFT;
  delete xi_t_FFT;
  delete xi_n_FFT;
  delete eta_n_FFT;
  delete alpha_n_IFFT;
  delete beta_n_IFFT;
  delete xi_0_IFFT;
  delete eta_0_IFFT;
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
  double* Phi;  ///< Real output vector. The first column of inverse Toeplitz matrix.
  double ldV;   ///< Real number. Log-determinant of Toeplitz matrix.

  /// Constructor.
  GSchurN(int, int);
  /// Destructor.
  ~GSchurN();
  /// Perform the Generalized Schur algorithm on the input data.
  void Compute(double*);
};

/// @param[in] N_ Size of Toeplitz matrix (Size of Generalized Schur Algorithm
/// is `N_-1`).
/// @param[in] b_ Integer for binary modulus.
inline GSchurN::GSchurN(int N_, int b_) {
  N = N_;
  b = b_;
  alpha = new double[N - 1];
  beta = new double[N - 1];
  Phi = new double[N];

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
  delete[] Phi;

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
/// @param[in] n2 Length of zero-padding.
/// @return Output complex vectors for \f$\alpha_{n,n} and \beta_{n,n} \f$ stored in `gsr->betam_IFF->out` and `gsr->betam_IFF->out`
inline void GSchurN::alpha2Beta(GSchur2K* gsr, int n1, int n2) {
  gsr->alpha_0_FFT->fft();
  gsr->beta_0_FFT->fft();
  gsr->eta_0_FFT->fft();
  gsr->xi_0_FFT->fft();

  if (n1 - n2) {
	  for (int ii = 0; ii < n1; ++ii) {
		  gsr->eta_t_FFT->in[n1 - ii] = gsr->eta_0_FFT->in[ii];
		  gsr->xi_t_FFT->in[n1 - ii] = gsr->xi_0_FFT->in[ii];
	  }

	  gsr->eta_t_FFT->fft();
	  gsr->xi_t_FFT->fft();
  }
  else {
	  for (int ii = 0; ii < n1; ++ii) {
		  gsr->eta_t_FFT->out[2 * ii][0] = gsr->eta_0_FFT->out[2 * ii][0];
		  gsr->eta_t_FFT->out[2 * ii][1] = -gsr->eta_0_FFT->out[2 * ii][1];
		  gsr->eta_t_FFT->out[2 * ii + 1][0] = -gsr->eta_0_FFT->out[2 * ii + 1][0];
		  gsr->eta_t_FFT->out[2 * ii + 1][1] = gsr->eta_0_FFT->out[2 * ii + 1][1];

		  gsr->xi_t_FFT->out[2 * ii][0] = gsr->xi_0_FFT->out[2 * ii][0];
		  gsr->xi_t_FFT->out[2 * ii][1] = -gsr->xi_0_FFT->out[2 * ii][1];
		  gsr->xi_t_FFT->out[2 * ii + 1][0] = -gsr->xi_0_FFT->out[2 * ii + 1][0];
		  gsr->xi_t_FFT->out[2 * ii + 1][1] = gsr->xi_0_FFT->out[2 * ii + 1][1];
	  }
  }

  int nn = (n1 + n2) / 2 + 1;

  vecConv(gsr->alpha_n_IFFT->in, gsr->alpha_0_FFT->out, gsr->eta_0_FFT->out, nn);
  vecConv_Sub(gsr->alpha_n_IFFT->in, gsr->xi_0_FFT->out, gsr->beta_0_FFT->out, nn);
  vecConv(gsr->beta_n_IFFT->in, gsr->beta_0_FFT->out, gsr->eta_t_FFT->out, nn);
  vecConv_Sub(gsr->beta_n_IFFT->in, gsr->xi_t_FFT->out, gsr->alpha_0_FFT->out, nn);

  gsr->alpha_n_IFFT->Ifft();
  gsr->beta_n_IFFT->Ifft();
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
/// @param[in] n Length of inputs (integer).
/// @return Output complex vectors for \f$ \eta_{2n}, \xi_{2n} \f$ stored in `gsr->xi_0_IFFT->out` and `gsr->eta_0_IFFT->out`.
inline void GSchurN::eta2Xi(GSchur2K* gsr, int n) {
  gsr->xi_n_FFT->fft();
  gsr->eta_n_FFT->fft();

  int nn = n / 2 + 1;

  vecConv(gsr->xi_0_IFFT->in, gsr->eta_t_FFT->out, gsr->xi_n_FFT->out, nn);
  vecConv_Add(gsr->xi_0_IFFT->in, gsr->xi_0_FFT->out, gsr->eta_n_FFT->out, nn);
  vecConv(gsr->eta_0_IFFT->in, gsr->xi_t_FFT->out, gsr->xi_n_FFT->out, nn);
  vecConv_Add(gsr->eta_0_IFFT->in, gsr->eta_0_FFT->out, gsr->eta_n_FFT->out, nn);

  gsr->xi_0_IFFT->Ifft();
  gsr->eta_0_IFFT->Ifft();
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
  double *tmpPtr, *xi1, *xi2, *eta1, *eta2, alpha_0, alpha_1;

  // Initialize the memory
  std::fill(gs[0]->xi_n_FFT->in, gs[0]->xi_n_FFT->in + 2 * n, 0);
  std::fill(gs[0]->eta_n_FFT->in, gs[0]->eta_n_FFT->in + 2 * n, 0);
  // Take advantage of the memory of existing `GSchur2K` class
  eta1 = gs[0]->eta_n_FFT->in;
  eta2 = gs[0]->eta_n_FFT->in + n;
  xi1 = gs[0]->xi_n_FFT->in;
  xi2 = gs[0]->xi_n_FFT->in + n;

  eta1[0] = 1.0;
  xi1[0] = alpha0[0] / beta0[0];
  gs[0]->gamma[0] = xi1[0];
  gs[0]->beta_0_FFT->in[0] = beta0[0] * (1 - xi1[0] * xi1[0]);

  for (int kk = 1; kk < n; ++kk) {
    alpha_1 = alpha0[kk];
    gs[0]->beta_0_FFT->in[kk] = beta0[kk];
    for (int jj = 1; jj <= kk; ++jj) {
      alpha_0 =
          alpha_1 - gs[0]->gamma[jj - 1] * gs[0]->beta_0_FFT->in[kk - jj + 1];
      gs[0]->beta_0_FFT->in[kk - jj + 1] -= gs[0]->gamma[jj - 1] * alpha_1;
      alpha_1 = alpha_0;
    }
    gs[0]->gamma[kk] = alpha_0 / gs[0]->beta_0_FFT->in[0];
    gs[0]->beta_0_FFT->in[0] *= 1 - gs[0]->gamma[kk] * gs[0]->gamma[kk];

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
  std::copy(xi1, xi1 + n, gs[0]->xi_0_IFFT->out);
  std::copy(eta1, eta1 + n, gs[0]->eta_0_IFFT->out);
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
/// stored in `gs[layer]->eta_0_IFFT->out`, `gs[layer]->xi_0_IFFT->out`,
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
    std::copy(alpha0, alpha0 + 2 * n, gs[m + 1]->alpha_0_FFT->in);
    std::copy(beta0, beta0 + 2 * n, gs[m + 1]->beta_0_FFT->in);
    std::copy(gs[m]->gamma, gs[m]->gamma + n, gs[m + 1]->gamma);
    std::copy(gs[m]->eta_0_IFFT->out, gs[m]->eta_0_IFFT->out + n,
              gs[m + 1]->eta_0_FFT->in);
    std::copy(gs[m]->xi_0_IFFT->out, gs[m]->xi_0_IFFT->out + n,
              gs[m + 1]->xi_0_FFT->in);
    // Compute alpha-nn and beta-nn
    alpha2Beta(gs[m + 1], n, n);
    // Compute eta-nn, xi-nn and gamma-nn using function GenStep
    GenStep(gs[m + 1]->alpha_n_IFFT->out + n, gs[m + 1]->beta_n_IFFT->out + n,
            n, m);
    // Prepare the inputs for function eta2Xi, including xi-nn and eta-nn
    std::copy(gs[m]->xi_0_IFFT->out, gs[m]->xi_0_IFFT->out + n,
              gs[m + 1]->xi_n_FFT->in);
    std::copy(gs[m]->eta_0_IFFT->out, gs[m]->eta_0_IFFT->out + n,
              gs[m + 1]->eta_n_FFT->in);
    // Compute eta-02n and xi-02n
    eta2Xi(gs[m + 1], 2 * n);
    // Paste the computed gamma-nn after gamma-0n to get gamma-02n
    std::copy(gs[m]->gamma, gs[m]->gamma + n, gs[m + 1]->gamma + n);
    n <<= 1;
  }
  return;
}

/// @return polynomials \f$ \{ \eta_{0,N}, \xi_{0,N}, \gamma_{0,N}\} \f$ stored
/// in `gsM[0]->eta_0_IFFT->out`, `gsM[0]->xi_0_IFFT->out`, `gsM[0]->gamma`.
inline void GSchurN::GenMerge() {
  int layer = log2(ceil((double)s[0] / b));

  GenStep(alpha, beta, s[0], layer);

  // When vector `s` has only one element, directly pass the results to output.
  if (k == 1) {
    std::copy(gs[layer]->eta_0_IFFT->out, gs[layer]->eta_0_IFFT->out + s[0],
              gsM[0]->eta_0_IFFT->out);
    std::copy(gs[layer]->xi_0_IFFT->out, gs[layer]->xi_0_IFFT->out + s[0],
              gsM[0]->xi_0_IFFT->out);
    std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
    return;
  }

  int n = 0;

  for (int m = 0; m < k - 1; ++m) {
	  n += s[m];

	  if (m) {
		  std::copy(gsM[m - 1]->eta_0_IFFT->out, gsM[m - 1]->eta_0_IFFT->out + n,
			  gsM[m]->eta_0_FFT->in);
		  std::copy(gsM[m - 1]->xi_0_IFFT->out, gsM[m - 1]->xi_0_IFFT->out + n,
			  gsM[m]->xi_0_FFT->in);
	  }
	  else {
		  std::copy(gs[layer]->eta_0_IFFT->out, gs[layer]->eta_0_IFFT->out + n,
			  gsM[m]->eta_0_FFT->in);
		  std::copy(gs[layer]->xi_0_IFFT->out, gs[layer]->xi_0_IFFT->out + n,
			  gsM[m]->xi_0_FFT->in);
		  std::copy(gs[layer]->gamma, gs[layer]->gamma + n, gsM[k - 2]->gamma);
	  }
	  
	  // First paste alpha_{0,n} and beta_{0,n} into gsM[0]->alpha and
	  // gsM[0]->beta, length is n
	  std::copy(alpha, alpha + n + s[m + 1], gsM[m]->alpha_0_FFT->in);
	  std::copy(beta, beta + n + s[m + 1], gsM[m]->beta_0_FFT->in);
	  
	  // alpha2Beta uses gsM[m]->alpha_0_FFT, beta_0_FFT, eta_0_FFT and xi_0_FFT
	  // to compute alpha_n_IFFT and beta_n_IFFT
	  alpha2Beta(gsM[m], n, s[m + 1]);
	  // produce alpha_n_IFFT->out [s[m]+(1:s[m])]
	  
	  // Computation of GenStep happens inside gs[0,1,...,layer].
	  // We just need alpha_n_IFFT->out [s[m]+(1:s[m+1])]
	  layer = log2(ceil((double)s[m + 1] / b));
	  GenStep(gsM[m]->alpha_n_IFFT->out + n, gsM[m]->beta_n_IFFT->out + n,
		  s[m + 1], layer);

	  // eta0m_t is gsM[m]->eta_t_FFT
	  // xi0m_t is gsM[m]->eta_t_FFT
	  // eta0m is gsM[m]->eta_0_FFT
	  // xi0m is gsM[m]->eta_0_FFT	  
	  std::copy(gs[layer]->eta_0_IFFT->out, gs[layer]->eta_0_IFFT->out + s[m + 1],
		  gsM[m]->eta_n_FFT->in);
	  std::copy(gs[layer]->xi_0_IFFT->out, gs[layer]->xi_0_IFFT->out + s[m + 1],
		  gsM[m]->xi_n_FFT->in);
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
  
  // convert eta and xi into Phi, the first column of inverse Toeplitz matrix.
  std::copy(gsM[k2 - 2]->eta_0_IFFT->out, gsM[k2 - 2]->eta_0_IFFT->out + N - 1, Phi);
  Phi[N - 1] = 0;
  Phi[0] /= sigma2;
  for (int ii = 1; ii < N; ++ii) {
    Phi[ii] += gsM[k2 - 2]->xi_0_IFFT->out[ii - 1];
    Phi[ii] /= sigma2;
  }
  return;
}

#endif
