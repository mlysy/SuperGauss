///////////////////////////////////////////////

// Generalized Schur Algorithm

///////////////////////////////////////////////

#include "VectorFFT.h"

// defining classes
//------------------------------------------------------

// 3, generalized schur algorithm
class GSchur
{
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
	   GSchur(int);
	   ~GSchur();
};
//------------------------------------------------------


// functions for convolution 
void CompMult(fftw_complex*, fftw_complex*, fftw_complex*, int);

void CompMult_Minus(fftw_complex*, fftw_complex*, fftw_complex*, int);

void CompMult_Plus(fftw_complex*, fftw_complex*, fftw_complex*, int);

vector<int> intoBit(int, int);

void TraceL(double*, int, double);

//------------------------------------------------------
// two helper functions applied in gschur
void alpha_beta(GSchur*, int);

void Merge(GSchur*, int);

void TraceComp(double*, int, double&);
//------------------------------------------------------
// 4, InverseToeplitz
class InverseToeplitz
{
	int n; // indicating the size of input
	int base; // ideal base is 64
    double* alpha;
	double* beta;
	vector<int> s;
	GSchur** gs;
	GSchur** gsM;
	void Pschur(double*, double*, int);
	void Gschur(double*, double*, int, int);
	void GschurMerge();
public:
	InverseToeplitz(int, int);
	~InverseToeplitz();
	double* Phi;
	double ldV;
	void Inverse(double*);
};