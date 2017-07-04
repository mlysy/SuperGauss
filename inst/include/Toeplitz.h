///////////////////////////////////////////////

// Toeplitz
// In this version, trace part is not available
// can only returns T^{-1} * x and T * x and det(T)
///////////////////////////////////////////////

#ifndef Toeplitz_h
#define Toeplitz_h 1

// usual header

#include "GSchur.h"

// defining classes
//------------------------------------------------------

// InverseProduct
class Toeplitz {
  int n; // size of matrix
  // int d;
  // fft space
  VectorFFT* L1fft;
  VectorFFT* L11fft;
  VectorFFT* L2fft;
  VectorFFT* L22fft;
  VectorFFT* xfft;
  VectorFFT* Lxfft;
  VectorFFT* U1fft;
  VectorFFT* U2fft;
  VectorFFT* Toepfft;
  // inverse-fft space
  VectorIFFT* Invfft;
  double* phi2;
  double* temVec;

  // Flag controling the inner loop
  bool has_acf; // has acf been set yet
  bool has_mult; // have multiplication fft's been computed yet
  bool has_solve; // has Gohberg-Semencul decomposition been computed yet
  //bool acf_is_0;

  GSchurN* Gs; // Gohberg-Sememcul decomposition
  double* acf; // acf of the Toeplitz matrix

  // trace of L*U, triangular Toeplitz matrices
  double trace_LU(double* acf1, double* acf2, int n);
  bool all_zeros(double* acf, int n); // check if vector is all zeros
  void mult_setup(); // one-time calcs for multiplication
  void solve_setup(); // one-time calcs for solving linear systems
 public:
  Toeplitz(int n_); // constructor
  ~Toeplitz(); // destructor
  void setAcf(double* acfIn); // acf input
  void getAcf(double* acfOut); // acf output
  int size(); // size of the matrix
  bool hasAcf(); // whether or not the acf has been set
  void multVec(double* yOut, double* xIn); // Toeplitz * vector
  void solveVec(double* yOut, double* xIn); // Toeplitz^-1 * vector
  double logDet(); // log-determinant
  // trace(Toeplitz^-1 * Toeplitz_2)
  double traceProd(double* acf2);
  // trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
  double traceDeriv(double* acf2, double* acf3);
  void getPhi(double* phiOut); // first column of inverse matrix
};

//------------------------------------------------------

inline Toeplitz::Toeplitz(int n_) {
  n = n_;
  Gs = new GSchurN(n, 64); // default base 64
  L1fft = new VectorFFT(2 * n);
  L11fft = new VectorFFT(2 * n);
  L2fft = new VectorFFT(2 * n);
  L22fft = new VectorFFT(2 * n);
  xfft = new VectorFFT(2 * n);
  Lxfft = new VectorFFT(2 * n);
  U1fft = new VectorFFT(2 * n);
  U2fft = new VectorFFT(2 * n);
  Toepfft = new VectorFFT(2 * n);
  Invfft = new VectorIFFT(2 * n);
  acf = new double[n];
  phi2 = new double[n];
  temVec = new double[n];
  has_acf = false;
}

inline Toeplitz::~Toeplitz() {

  delete Gs;
  delete L1fft;
  delete L11fft;
  delete L2fft;
  delete L22fft;
  delete xfft;
  delete Lxfft;
  delete U1fft;
  delete U2fft;
  delete Toepfft;
  delete Invfft;
  delete[] acf;
  delete[] phi2;
  delete[] temVec;
}

// calculating the trace of L1 L2', lower trangular toeplitz composed of acf of size n
inline double Toeplitz::trace_LU(double* acf1, double* acf2, int n) {
  double trace = 0;
  for(int ii = 0; ii < n; ++ii){
    trace += (n - ii) * acf1[ii] * acf2[ii];
  }
  return trace;
}

inline bool Toeplitz::all_zeros(double* acf, int n) {
  bool acf_is_0 = true;
  for(int ii = 0; ii < n; ++ii){
    if(fabs(acf[ii]) > 0.0){
      acf_is_0 = false;
      break;
    }
  }
  return acf_is_0;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::setAcf(double* acfIn) {
  std::copy(acfIn, acfIn + n, acf);
  has_acf = true;
  has_mult = false;
  has_solve = false;
  return;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::getAcf(double *acfOut) {
  std::copy(acf, acf + n, acfOut);
  return;
}

//-------------------------------------------------------------------------------------------------
inline bool Toeplitz::hasAcf() {
  return has_acf;
}

//-------------------------------------------------------------------------------------------------
inline int Toeplitz::size() {
  return n;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::getPhi(double* phiOut) {
  std::copy(Gs->Phi, Gs->Phi + n, phiOut);
  return;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::mult_setup(){
  std::copy(acf, acf + n, Toepfft->in);
  std::copy(acf + 1, acf + n, Toepfft->in + n + 1);
  std::reverse(Toepfft->in + n + 1, Toepfft->in + 2 * n);
  Toepfft->fft();
  has_mult = true;
  return;
}

// Toeplitz matrix vector multiplication
inline void Toeplitz::multVec(double* yOut, double* xIn) {
  bool acf_is_0 = all_zeros(acf, n);

  if(acf_is_0) {
    std::fill(yOut, yOut + n, 0);
    return;
  }

  if(!has_mult){
    mult_setup();
  }

  std::copy(xIn, xIn + n, xfft->in);
  xfft->fft();
  vecConv(Invfft->in, Toepfft->out, xfft->out, 2*(n/2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, yOut);
  return;
}

//-------------------------------------------------------------------------------------------------
// compute part, generate required parts for TraceProd and InverseMult
inline void Toeplitz::solve_setup() {

  Gs->Compute(acf);

  L11fft->in[0] = Gs->Phi[0];
  std::copy(Gs->Phi + 1, Gs->Phi + n, L11fft->in + n + 1);
  std::reverse(L11fft->in + n + 1, L11fft->in + 2 * n);
  L11fft->fft(); // L1'

  std::copy(Gs->Phi, Gs->Phi + n, L1fft->in);
  L1fft->fft(); // L1

  std::copy(Gs->Phi + 1, Gs->Phi + n, L22fft->in + n + 1);
  L22fft->fft(); // L2'

  std::copy(Gs->Phi + 1, Gs->Phi + n, L2fft->in + 1);
  std::reverse(L2fft->in + 1, L2fft->in + n);
  L2fft->fft(); // L2

  has_solve = true;
  return;
}


inline void Toeplitz::solveVec(double* yOut, double* xIn) {
  if(!has_solve){
    solve_setup();
  }

  std::copy(xIn, xIn + n, xfft->in);
  xfft->fft();

  //L1 * L1' * x
  vecConv(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
  Lxfft->fft();
  vecConv(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  // stored in yOut
  std::copy(Invfft->out, Invfft->out + n, yOut);

  // L2 * L2' * x
  vecConv(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
  Lxfft->fft();
  vecConv(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  // 1/sigma2 (L1 L1'x - L2 L2'x)
  for(int ii = 0; ii < n; ++ii){
    yOut[ii] -= Invfft->out[ii];
    yOut[ii] /= Gs->Phi[0];
  }
  return;
}

inline double Toeplitz::logDet() {
  if(!has_solve){
    solve_setup();
  }
  return Gs->ldV;
}


// trace(Toeplitz^-1 * Toeplitz_2)
inline double Toeplitz::traceProd(double* acf2) {
  bool is_small_acf20;
  double t0;
  double trace;
  if(all_zeros(acf2, n)) {
    return 0.0;
  }

  if(!has_solve){
    solve_setup();
  }

  std::copy(acf2 + 1, acf2 + n, U1fft->in + 1);
  // check first term
  is_small_acf20 = fabs(acf2[0]) < 0.0001;
  t0 = is_small_acf20 ? acf2[0] + 1.0 : acf2[0];

  U1fft->in[0] = t0;
  U1fft->fft(); //U1
  std::copy(acf2 + 1, acf2 + n, U2fft->in + 1);
  U2fft->fft(); //U2

  // tr{U1'L1L1'U1}
  vecConv(Invfft->in, L1fft->out, U1fft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  trace = trace_LU(Invfft->out, Invfft->out, n);

  // tr{U1'L2L2'U1}
  vecConv(Invfft->in, L2fft->out, U1fft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  trace -= trace_LU(Invfft->out, Invfft->out, n);

  // tr{U2'L1L1'U2}
  vecConv(Invfft->in, L1fft->out, U2fft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  trace -= trace_LU(Invfft->out, Invfft->out, n);

  // tr{U2'L2L2'U2}
  vecConv(Invfft->in, L2fft->out, U2fft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  trace += trace_LU(Invfft->out, Invfft->out, n);

  // trace
  trace /= Gs->Phi[0];
  trace /= t0;

  if(is_small_acf20) {
    t0 = 0.0;
    for(int ii=1; ii<n; ii++) {
      t0 += (n-2*ii) * Gs->Phi[ii] * Gs->Phi[ii];
    }
    trace -= n * Gs->Phi[0] + t0/Gs->Phi[0];
  }
  return trace;
}


// trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
// TODO: this need to be corrected for when acf3[0] == 0.
inline double Toeplitz::traceDeriv(double* acf2, double* acf3) {
  double trace2;
  double t0;
  bool is_small_acf30;

  if(all_zeros(acf2, n) | all_zeros(acf3, n)) {
    return 0.0;
  }

  if(!has_solve) {
    solve_setup();
  }

  // check first term
  is_small_acf30 = fabs(acf3[0]) < 0.0001;
  t0 = is_small_acf30 ? acf3[0] + 1.0 : acf3[0];

  // phi2 = - Toeplitz^-1 * Toeplitz_2 * phi ------------------------------
  // phi2 = Toeplitz_2 * phi
  std::copy(acf2, acf2 + n, Lxfft->in);
  std::copy(acf2 + 1, acf2 + n, Lxfft->in + n + 1);
  std::reverse(Lxfft->in + n + 1, Lxfft->in + 2 * n);
  Lxfft->fft();
  std::fill(Lxfft->in + n, Lxfft->in + 2 * n, 0); // (readjust back to 0 assignment for latter half)
  std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
  xfft->fft();
  vecConv(Invfft->in, Lxfft->out, xfft->out, 2*(n/2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, phi2);

  // phi2 = - Toeplitz^-1 * phi2
  std::copy(phi2, phi2 + n, xfft->in);
  xfft->fft();

  //L1 * L1' * x
  vecConv(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
  Lxfft->fft();
  vecConv(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  std::copy(Invfft->out, Invfft->out + n, phi2);

  //L2 * L2' * x
  vecConv(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();
  std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
  Lxfft->fft();
  vecConv(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  // phi2 = -1/phi[1] * (L1 L1'x - L2 L2'x)
  for(int ii = 0; ii < n; ++ii) {
    phi2[ii] -= Invfft->out[ii];
    phi2[ii] /= Gs->Phi[0];
  }
  // phi2 is obtained ----------------------------------------------------

  // tr = phi2[1] * tr(Toeplitz^-1 * Toeplitz_3) + 2 * tr(L1(-phi2) * L1(phi)'
  //      * Toeplitz_3) - 2 * tr(_L2(-phi2) * _L2(phi)' * Toeplitz_3)


  // 1, phi2[1] * tr(Toeplitz^-1 * Toeplitz_3)
  trace2 = phi2[0] * traceProd(acf3);

  // 2, tr(L1(-phi2) * L1(phi)' * Toeplitz_3) = -tr(L1(phi2) * L1(phi)' * L1'(acf3)
  //    * L1(acf3)) + tr(L1(phi2) * L1(phi)' * L2'(acf3) * L2(acf3))

  // 2.1, tr(L1(phi2) * L1(phi)' * L1'(acf3) * L1(acf3)) = trace_LU(L1(acf3) * L1(phi2),
  //      L1(acf3) * L1(phi)) = temp
  std::copy(acf3, acf3 + n, U1fft->in);
  U1fft->in[0] = t0;
  U1fft->fft();

  std::copy(phi2, phi2 + n, Lxfft->in);
  Lxfft->fft();
  vecConv(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  std::copy(Invfft->out, Invfft->out + n, temVec);

  std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
  xfft->fft();
  vecConv(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  trace2 -= 2 * trace_LU(Invfft->out, temVec, n) / t0;

  // 2.2, tr(L1(phi2) * L1(phi)' * L2'(acf3) * L2(acf3)) = trace_LU(L2(acf3) * L1(phi2),
  //      L2(acf3) * L1(phi)) = temp
  std::copy(acf3 + 1, acf3 + n, U2fft->in + 1);
  U2fft->in[0] = 0;
  U2fft->fft();

  vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  std::copy(Invfft->out, Invfft->out + n, temVec);

  vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  trace2 += 2 * trace_LU(Invfft->out, temVec, n) / t0;

  // 3, tr(_L2(phi2) * _L2(phi)' * Toeplitz_3) = tr(_L2(phi2) * _L2'(phi) * L1'(acf3)
  //    * L1(acf3)) / acf3[1] - tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) / acf3[1]

  // 3.1, tr(L2(phi2) * L2'(phi) * L1'(acf3) * L1(acf3)) = trace_LU(L1(acf3) * _L2(phi2),
  //      L1(acf3) * _L2(phi)) = temp
  Lxfft->in[0] = 0;
  std::copy(phi2 + 1, phi2 + n, Lxfft->in + 1);
  std::reverse(Lxfft->in + 1, Lxfft->in + n);
  Lxfft->fft();
  vecConv(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  std::copy(Invfft->out, Invfft->out + n, temVec);

  xfft->in[0] = 0;
  std::copy(Gs->Phi + 1, Gs->Phi + n, xfft->in + 1);
  std::reverse(xfft->in + 1, xfft->in + n);
  xfft->fft();
  vecConv(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  trace2 += 2 * trace_LU(Invfft->out, temVec, n) / t0;

  // 3.2, tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) = trace_LU(L2(acf3) * _L2(phi2),
  //      L2(acf3) * _L2(phi)) = temp
  vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  std::copy(Invfft->out, Invfft->out + n, temVec);

  vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
  Invfft->Ifft();

  trace2 -= 2 * trace_LU(Invfft->out, temVec, n) / t0;

  // tr = tr / phi[1]
  trace2 /= Gs->Phi[0];
  trace2 *= -1.0;

  if(is_small_acf30) {
    t0 = 0.0;
    for(int ii=1; ii<n; ii++) {
      t0 += (n-2*ii) * Gs->Phi[ii] * phi2[ii];
    }
    trace2 -= 2 * n * phi2[0] + 2 * t0/Gs->Phi[0];
  }

  return trace2;
}

# endif
