/// @file GSchur.h

///////////////////////////////////////////////

// Generalized Schur Algorithm

///////////////////////////////////////////////

#ifndef GSchur_h
#define GSchur_h 1

// !!!FIXME: remove "using namespace" directives from VectorFFT.h!!!
#include "VectorFFT.h"

/// Convert integer to modulo-binary representation.
///
/// Given an integer `x`, returns an integer vector `s = {2^k0 * base, 2^k1 * base, ..., 2^kT * base, rem}`, where `k0 > ... > kT >= 0`, `rem < base`, and `sum(s) = x`.
///
/// @param[in] x Integer of which the modulo-binary represention is computed.
/// @param[in] base Integer giving the binary modulus.
/// @return Integer vector containing the modulo-binary representation.
inline vector<int> int2Bin(int x, int base = 1) {
  vector<int> s;
  int x1 = x / base;
  int dif = x - x1 * base;
  int M = base;
  do{
    if (x1 & 1){
      s.push_back(M);
    }
    M <<= 1;
  } while ((x1 >>= 1)>0);
  reverse(s.begin(), s.end());
  if (dif){
    s.push_back(dif);
  }
  return s;
}

// defining classes
//------------------------------------------------------

// 3, generalized schur algorithm
class GSchur2K {
 public:
  VectorFFT* alpha_FFT;
  VectorFFT* beta_FFT;
  VectorFFT* eta_FFT;
  VectorFFT* xi_FFT;
  VectorFFT* xi_m_FFT;
  VectorFFT* eta_m_FFT;
  VectorIFFT* alpham_IFFT;
  VectorIFFT* betam_IFFT;
  VectorIFFT* xi_IFFT;
  VectorIFFT* eta_IFFT;
  double* gamma;
  fftw_complex* eta_t;
  fftw_complex* xi_t;

  GSchur2K(int);
  ~GSchur2K();
};

//--------------------------------------------------------------------------------------------------
// constructor and destructor of class GSchur2K
inline GSchur2K::GSchur2K(int M) {
  alpha_FFT = new VectorFFT(M);
  beta_FFT = new VectorFFT(M);
  eta_FFT = new VectorFFT(M);
  xi_FFT = new VectorFFT(M);
  xi_m_FFT = new VectorFFT(M);
  eta_m_FFT = new VectorFFT(M);
  alpham_IFFT = new VectorIFFT(M);
  betam_IFFT = new VectorIFFT(M);
  xi_IFFT = new VectorIFFT(M);
  eta_IFFT = new VectorIFFT(M);
  gamma = new double[M];
  eta_t = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)  * M);
  xi_t = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)  * M);
}


inline GSchur2K::~GSchur2K() {
  delete alpha_FFT;
  delete beta_FFT;
  delete eta_FFT;
  delete xi_FFT;
  delete xi_m_FFT;
  delete eta_m_FFT;
  delete alpham_IFFT;
  delete betam_IFFT;
  delete xi_IFFT;
  delete eta_IFFT;
  delete[]gamma;
  fftw_free(eta_t);
  fftw_free(xi_t);
}
//------------------------------------------------------

// 4, GSchurN
class GSchurN {
  int n; // indicating the size of input
  int base; // ideal base is 64
  double* alpha;
  double* beta;
  vector<int> s;
  GSchur2K** gs;
  GSchur2K** gsM;

  void alpha2Beta(GSchur2K*, int);
  void GSchur_Merge(GSchur2K*, int);

  void ProgStep(double*, double*, int);
  void GenStep(double*, double*, int, int);
  void GenStep_Merge();
 public:
  GSchurN(int, int);
  ~GSchurN();
  double* Phi;
  double ldV;
  void Compute(double*);
};

//---------------------------------------------------------------------------------------

// given representation s(x) = 2^k_0 * base + ... + 2^k_T * base + rem
// gs is an array of length k_0 + 1,
// where element t=0...k_0 contains fft storage of size base * 2^t.
// gsM is an array of length T+1,
// each should contain fft storage for the merge at that step,
// starting at the end?
// FIXME: get storage exactly right (currently it's next power of 2).
// FIXME: doesn't work for n_ = 0.
// FIXME: remove size_t => int conversions
inline GSchurN::GSchurN(int n_, int base_) {
  n = n_;
  base = base_;
  alpha = new double[n - 1];
  beta = new double[n - 1];
  s = int2Bin(n - 1, base);
  Phi = new double[n];
  int jj = base;
  gs = new GSchur2K*[(int)log2(ceil((double)s[0]/base)) + 1];
  gs[0] = new GSchur2K(jj<<1);
  for (int ii = 0; ii < (int)log2(ceil((double)s[0]/base)); ++ii){
    gs[ii + 1] = new GSchur2K(jj<<1);
    jj <<= 1;
  }
  gsM = new GSchur2K*[s.size()];
  for(int ii = 0; ii < (int)s.size(); ++ii) {
    gsM[ii] = new GSchur2K(s[ii]<<1);
  }
}

inline GSchurN::~GSchurN() {
  delete[] alpha;
  delete[] beta;
  delete[] Phi;

  for (int ii = 0; ii < (int)s.size(); ++ii){
    delete gsM[ii];
  }
  delete[] gsM;

  for (int ii = 0; ii <= (int)log2(ceil((double)s[0]/base)); ++ii){
    delete gs[ii];
  }
  delete[] gs;
}

//---------------------------------------------------------------------------------------
// helper functions applied in GSchur2K and GSchurMerge
// generate the alpham_IFFT and betam_IFFT, the size is M
// input of this function is gsr->alpha/beta/eta/xi_FFT->in
inline void GSchurN::alpha2Beta(GSchur2K* gsr, int M) {
  gsr->alpha_FFT->fft();
  gsr->beta_FFT->fft();
  gsr->eta_FFT->fft();
  gsr->xi_FFT->fft();

  for (int ii = 0; ii < M; ++ii){
    gsr->eta_t[2 * ii][0] = gsr->eta_FFT->out[2 * ii][0];
    gsr->eta_t[2 * ii][1] = -gsr->eta_FFT->out[2 * ii][1];
    gsr->eta_t[2 * ii + 1][0] = -gsr->eta_FFT->out[2 * ii + 1][0];
    gsr->eta_t[2 * ii + 1][1] = gsr->eta_FFT->out[2 * ii + 1][1];

    gsr->xi_t[2 * ii][0] = gsr->xi_FFT->out[2 * ii][0];
    gsr->xi_t[2 * ii][1] = -gsr->xi_FFT->out[2 * ii][1];
    gsr->xi_t[2 * ii + 1][0] = -gsr->xi_FFT->out[2 * ii + 1][0];
    gsr->xi_t[2 * ii + 1][1] = gsr->xi_FFT->out[2 * ii + 1][1];
  }
  // ifft

  vecConv(gsr->alpham_IFFT->in, gsr->alpha_FFT->out, gsr->eta_FFT->out, 2 * (M / 2 + 1));
  vecConv_Sub(gsr->alpham_IFFT->in, gsr->xi_FFT->out, gsr->beta_FFT->out, 2 * (M / 2 + 1));

  vecConv(gsr->betam_IFFT->in, gsr->beta_FFT->out, gsr->eta_t, 2 * (M / 2 + 1));
  vecConv_Sub(gsr->betam_IFFT->in, gsr->xi_t, gsr->alpha_FFT->out, 2 * (M / 2 + 1));
  gsr->alpham_IFFT->Ifft();
  gsr->betam_IFFT->Ifft();
}

// merging the pervious xi/eta and current xi/eta
// input is gsr->eta/xi_t, gsr->xi/eta_FFT->out, gsr->xi/eta_m_FFT->in, output is gsr->xi/eta_IFFT->out
inline void GSchurN::GSchur_Merge(GSchur2K* gsr, int M) {
  gsr->xi_m_FFT->fft();
  gsr->eta_m_FFT->fft();
  vecConv(gsr->xi_IFFT->in, gsr->eta_t, gsr->xi_m_FFT->out, 2 * (M / 2 + 1));
  vecConv_Add(gsr->xi_IFFT->in, gsr->xi_FFT->out, gsr->eta_m_FFT->out, 2 * (M / 2 + 1));

  vecConv(gsr->eta_IFFT->in, gsr->xi_t, gsr->xi_m_FFT->out, 2 * (M / 2 + 1));
  vecConv_Add(gsr->eta_IFFT->in, gsr->eta_FFT->out, gsr->eta_m_FFT->out, 2 * (M / 2 + 1));

  gsr->xi_IFFT->Ifft();
  gsr->eta_IFFT->Ifft();
}

//--------------------------------------------------------------------------------------------------
// pschur function
// the final result is stored in gs->gamma, ->eta_IFFT->out, ->xi_IFFT->out.
inline void GSchurN::ProgStep(double* alpha0, double* beta0, int size) {
  // alpha_FFT->in, beta_FFT->in, xi_m_FFTPS->in, eta_m_FFT->in, xi_IFFT->out, eta_IFFT->out, gamma are n
  std::fill(gs[0]->xi_m_FFT->in,  gs[0]->xi_m_FFT->in + 2*size, 0);
  std::fill(gs[0]->eta_m_FFT->in, gs[0]->eta_m_FFT->in + 2*size, 0);
  double *tmpPtr, *xi1, *xi2, *eta1, *eta2, alpha_0, alpha_1;
  eta1 = gs[0]->eta_m_FFT->in;
  eta2 = gs[0]->eta_m_FFT->in + size;
  xi1 = gs[0]->xi_m_FFT->in;
  xi2 = gs[0]->xi_m_FFT->in + size;
  eta1[0] = 1.0;
  xi1[0] = alpha0[0]/beta0[0];
  gs[0]->gamma[0] = xi1[0];
  gs[0]->beta_FFT->in[0] = beta0[0] * (1 - xi1[0] * xi1[0]);

  for (int kk = 1; kk < size; ++kk){
    alpha_1 = alpha0[kk];
    gs[0]->beta_FFT->in[kk] = beta0[kk];
    for (int jj = 1; jj <= kk; ++jj){
      alpha_0 = alpha_1 - gs[0]->gamma[jj - 1] * gs[0]->beta_FFT->in[kk - jj + 1];
      gs[0]->beta_FFT->in[kk - jj + 1] -= gs[0]->gamma[jj - 1] * alpha_1;
      alpha_1 = alpha_0;
    }
    gs[0]->gamma[kk] = alpha_0 / gs[0]->beta_FFT->in[0];
    gs[0]->beta_FFT->in[0] *= 1 - gs[0]->gamma[kk] * gs[0]->gamma[kk];

    eta2[0] = 1.0;
    xi2[0] = alpha0[0]/beta0[0];
    for (int jj = 1; jj <= kk; ++jj){
      xi2[jj] = xi1[jj] + gs[0]->gamma[kk] * eta1[kk - jj];
      eta2[jj] = eta1[jj] + gs[0]->gamma[kk] *xi1[kk - jj];
    }
    tmpPtr = xi1;
    xi1 = xi2;
    xi2 = tmpPtr;
    tmpPtr = eta1;
    eta1 = eta2;
    eta2 = tmpPtr;
  }
  std::copy(xi1, xi1 + size, gs[0]->xi_IFFT->out);
  std::copy(eta1, eta1 + size, gs[0]->eta_IFFT->out);
  return;
}

//------------------------------------------------------------------------------------------------------------


//gschur function
// "size" is the size of input
inline void GSchurN::GenStep(double* alpha0, double* beta0,
			     int size, int layer) {
  if(size <= base){
    ProgStep(alpha0, beta0, size);
    return;
  }

  ProgStep(alpha0, beta0, base);

  int M = base;
  for (int m = 0; m < layer; ++m){
    //copying information from layer[m] to [m+1]
    std::copy(alpha0, alpha0 + 2 * M, gs[m + 1]->alpha_FFT->in);
    std::copy(beta0, beta0 + 2 * M, gs[m + 1]->beta_FFT->in);
    std::copy(gs[m]->gamma, gs[m]->gamma + M, gs[m + 1]->gamma);
    std::copy(gs[m]->eta_IFFT->out, gs[m]->eta_IFFT->out + M, gs[m + 1]->eta_FFT->in);
    std::copy(gs[m]->xi_IFFT->out, gs[m]->xi_IFFT->out + M, gs[m + 1]->xi_FFT->in);

    //generating alpham and betam
    alpha2Beta(gs[m + 1], M);
    //-----------------------------------------------------------------------------------------------
    // returns in layer[m]
    GenStep(gs[m + 1]->alpham_IFFT->out + M, gs[m + 1]->betam_IFFT->out + M, M, m);
    //-----------------------------------------------------------------------------------------------
    //copying information from layer[m] to [m+1]
    std::copy(gs[m]->gamma, gs[m]->gamma + M, gs[m + 1]->gamma + M);
    std::copy(gs[m]->xi_IFFT->out, gs[m]->xi_IFFT->out + M, gs[m + 1]->xi_m_FFT->in);
    std::copy(gs[m]->eta_IFFT->out, gs[m]->eta_IFFT->out + M, gs[m + 1]->eta_m_FFT->in);

    //merging xi and eta
    GSchur_Merge(gs[m + 1], M);

    M <<= 1;
  }
  return;
}

//------------------------------------------------------------------------------------------------------------

//gschur.merge function
inline void GSchurN::GenStep_Merge() {
  int layer = (int)log2(ceil((double)s[0]/base));
  GenStep(alpha, beta, s[0], layer);
  if(s.size() == 1){
    std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[0], gsM[0]->eta_IFFT->out);
    std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[0], gsM[0]->xi_IFFT->out);
    std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
    return;
  }
  std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[0], gsM[0]->eta_FFT->in);
  std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[0], gsM[0]->xi_FFT->in);
  std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
  int M = s[0];
  for(int m = 0; m < (int)s.size()-1; ++m){
    if(m == 0){
      std::copy(alpha, alpha + n - 1, gsM[m]->alpha_FFT->in);
      std::copy(beta, beta + n - 1, gsM[m]->beta_FFT->in);
    }
    else{
      std::copy(gsM[m - 1]->alpham_IFFT->out + s[m - 1], gsM[m - 1]->alpham_IFFT->out + s[m - 1] + 2 * s[m], gsM[m]->alpha_FFT->in);
      std::copy(gsM[m - 1]->betam_IFFT->out + s[m - 1], gsM[m - 1]->betam_IFFT->out + s[m - 1] + 2 * s[m], gsM[m]->beta_FFT->in);
    }
    alpha2Beta(gsM[m], s[m]);

    layer = (int)log2(ceil((double)s[m + 1]/base));
    GenStep(gsM[m]->alpham_IFFT->out + s[m], gsM[m]->betam_IFFT->out + s[m], s[m + 1], layer);

    std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[m + 1], gsM[m + 1]->eta_FFT->in);
    std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[m + 1], gsM[m + 1]->xi_FFT->in);
    std::copy(gs[layer]->gamma, gs[layer]->gamma + s[m + 1], gsM[0]->gamma + M);

    M += s[m + 1];

  }
  M = s[(int)s.size() - 1];
  for(int m = (int)s.size()-2; m >= 0; --m){
    if(m == (int)s.size()-2){
      std::copy(gsM[m + 1]->xi_FFT->in, gsM[m + 1]->xi_FFT->in + M, gsM[m]->xi_m_FFT->in);
      std::copy(gsM[m + 1]->eta_FFT->in, gsM[m + 1]->eta_FFT->in + M, gsM[m]->eta_m_FFT->in);
    }
    else{
      std::copy(gsM[m + 1]->xi_IFFT->out, gsM[m + 1]->xi_IFFT->out + M, gsM[m]->xi_m_FFT->in);
      std::copy(gsM[m + 1]->eta_IFFT->out, gsM[m + 1]->eta_IFFT->out + M, gsM[m]->eta_m_FFT->in);
    }
    GSchur_Merge(gsM[m], s[m]);
    
    M += s[m];
  }
}

//output function
inline void GSchurN::Compute(double* acf) {
  for(int ii = 0; ii < n - 1; ++ii) {
    alpha[ii] = -1 * acf[ii + 1];
    beta[ii] = acf[ii];
  }
  GenStep_Merge();

  double sigma2 = log(acf[0]);
  ldV = sigma2;
  for (int ii = 0; ii < n - 1; ++ii) {
    if(gsM[0]->gamma[ii] < 1){
      sigma2 += log(1 - gsM[0]->gamma[ii] * gsM[0]->gamma[ii]);
      ldV += sigma2;
    }
  }
  sigma2 = exp(sigma2);

  std::copy(gsM[0]->eta_IFFT->out, gsM[0]->eta_IFFT->out + n - 1, Phi);
  Phi[n-1] = 0;
  Phi[0] /= sigma2;
  for(int ii = 1; ii < n; ++ii){
    Phi[ii] += gsM[0]->xi_IFFT->out[ii-1];
    Phi[ii] /= sigma2;

  }
  return;
}

# endif
