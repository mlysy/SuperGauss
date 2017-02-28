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
	// dimension
	int n;
	// int d;
	
	// fft space
	VectorFFT* L1fft; // pre-assigned in .computeInv: fft(phi, 0)
	VectorFFT* L11fft; // pre-assigned in .computeInv: transpose of L1fft
	VectorFFT* L2fft; // pre-assigned in .computeInv: fft(0, rev(phi[-1]), 0)
	VectorFFT* L22fft; // pre-assigned in .computeInv: transpose of L2fft

	VectorFFT* xfft; // free space
	VectorFFT* Lxfft; // free space
	
	VectorFFT* U1fft; // free space
	VectorFFT* U2fft; // free space
	
    VectorFFT* Toepfft; // pre-assigned in .computeMult: fft(acf, 0, rev(acf[-1]))

	// inverse-fft space
	VectorIFFT* Invfft; // free space

	//
	double* phi2;
	double* temVec;

	// Flag controling the inner loop
	bool hasMult;
	bool hasAcf;
	bool hasInv;
	bool acf_is_0;
	
	// Superfast Algorithm Space
	InverseToeplitz* Gs;

	// acf input
	double* acf;

public:		
	Toep(int);
	~Toep();

	// input the acf
	void acfInput(double*);

	// Toeplitz * vector
	void computeMult(); // controled by hasMult, if hasMult = F, run this ,otherwise skip
	void mult(double*); 

	// Toeplitz^-1 * vector 
	void computeInv(); // controled by hasInv, if hasInv = F, run this ,otherwise skip
	void solve(double*);
	void detCheck();
	
	// trace(Toeplitz^-1 * Toeplitz_2)
	void traceProd(double*);
	
	// trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
	void traceDerv(double*, double*);

	// return of mult and solve
	double* Mult;

	// log.deternimant
	double det;

	// return of TraceProd
	double trace;

	// return of TraceDerv
	double trace2;
};
//------------------------------------------------------