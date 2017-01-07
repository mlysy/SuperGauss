///////////////////////////////////////////////

// Toeplitz
// In this version, trace part is not available
// can only returns T^{-1} * x and T * x and det(T)
///////////////////////////////////////////////

// usual header

#include "GSchur.h"

// defining classes
//------------------------------------------------------

// InverseProduct
class Toep {
	int n;
	int d;
	
	// InverseProd part
	VectorFFT* L1fft; // used
	VectorFFT* L11fft; // used
	VectorFFT* L2fft; // used
	VectorFFT* L22fft; // used
	VectorFFT* xfft;
	VectorFFT* Lxfft;
	
	
	// Trace part
	VectorFFT* U1fft;
	VectorFFT* U2fft;
	
	// Toep mult part
    VectorFFT* Toepfft; // used

	// Trace derv part
	// VectorFFT* U1fft; // substituting L1Cfft
	// VectprFFT* U2fft; // substituting L2Cfft
    // VectorFFT* L1Dfft;
	// VectorFFT* L2Dfft;
	// VectorFFT* xfft; 

	// One InverseFFT is required
	VectorIFFT* Invfft;

	// double* phiD;

	// Flag controling the inner loop
	bool hasMult;
	bool hasAcf;
	bool hasInv;

	InverseToeplitz* Gs;

	double* acf;

	// void initialize(int n, int d);

public:		
	Toep(int, int); // second dimension default to be 1
	~Toep(); // destructor

	// input the acf
	void acfInput(double*);

	// Toeplitz times matrix
	void computeMult(); // controled by hasMult, if hasMult = F, run this ,otherwise skip
	void mult(double**); 

	// Inverse Toeplitz times matrix 
	void computeInv(); // controled by hasInv, if hasInv = F, run this ,otherwise skip
	void solve(double**);
	void detCheck();
	// tracr of Inverse Topelitz times Toeplitz
	void TraceProd(double*);
	
	// derivative of TraceProd
	// void TraceDerv(double*, double*); // not done yet

	// return of mult and solve, matrix n by d
	double** Mult;
	double det;

	// return of TraceProd and TraceDerv, scale
	double trace;

};
//------------------------------------------------------