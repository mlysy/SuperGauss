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
	// dimension
	int n;
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
	bool hasMult;
	bool hasSolve;
	bool acf_is_0;
	
	// Superfast Algorithm Space
	GSchurN* Gs;

	// acf input
	double* acf;

    double compTrace(double*, double*, int);
    bool zeroTest(double*, int);
    void flagMult();
    void flagSolve(); 
public:		
	Toeplitz(int);
	~Toeplitz();

	// input the acf
	void acfInput(double*);

	// Toeplitz * vector
	void multVec(double*, double*); 

	// Toeplitz^-1 * vector 
	void solveVec(double*, double*);
	
    // return the determinant
    double logDet();

	// trace(Toeplitz^-1 * Toeplitz_2)
	double traceProd(double*);
	
	// trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
	double traceDeriv(double*, double*);
};
//------------------------------------------------------

inline Toeplitz::Toeplitz(int n_){
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
}
 
inline Toeplitz::~Toeplitz(){
    
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
inline double Toeplitz::compTrace(double* acf1, double* acf2, int n){
    double trace = 0;
    for(int ii = 0; ii < n; ++ii){
        trace += (n - ii) * acf1[ii] * acf2[ii];
    }
    return trace;
}

inline bool Toeplitz::zeroTest(double* acf, int n){
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
inline void Toeplitz::acfInput(double* acfA){
    std::copy(acfA, acfA + n, acf);

    hasMult = FALSE;
    hasSolve = FALSE;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::flagMult(){
    std::copy(acf, acf + n, Toepfft->in);
    std::copy(acf + 1, acf + n, Toepfft->in + n + 1);
    std::reverse(Toepfft->in + n + 1, Toepfft->in + 2 * n);
    Toepfft->fft();
    hasMult = TRUE;
}

// Toeplitz matrix vector multiplication
inline void Toeplitz::multVec(double* yOut, double* xIn){
    bool acf_is_0 = zeroTest(acf, n);
    
    if(acf_is_0){
        std::fill(yOut, yOut + n, 0);
        return;
    }

    if(!hasMult){
        flagMult();
    }

    std::copy(xIn, xIn + n, xfft->in);
    xfft->fft();
    vecConv(Invfft->in, Toepfft->out, xfft->out, 2*(n/2 + 1));
    Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, yOut);
} 

//-------------------------------------------------------------------------------------------------
// compute part, generate required parts for TraceProd and InverseMult
inline void Toeplitz::flagSolve(){

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

    hasSolve = TRUE;
}


inline void Toeplitz::solveVec(double* yOut, double* xIn){
    if(!hasSolve){
        flagSolve();
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
}

inline double Toeplitz::logDet(){
    if(!hasSolve){
        flagSolve();
    }
    return Gs->ldV;
}


// trace(Toeplitz^-1 * Toeplitz_2)
inline double Toeplitz::traceProd(double* acf2){
    bool acf2_is_0 = zeroTest(acf2, n);
    if(acf2_is_0){
        return 0.0;
    }

    if(!hasSolve){
        flagSolve();
    }
    
    std::copy(acf2, acf2 + n, U1fft->in);
    U1fft->fft(); //U1
    std::copy(acf2 + 1, acf2 + n, U2fft->in + 1);
    U2fft->fft(); //U2

    double trace;

    // tr{U1'L1L1'U1}
    vecConv(Invfft->in, L1fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace = compTrace(Invfft->out, Invfft->out, n);

    // tr{U1'L2L2'U1}
    vecConv(Invfft->in, L2fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace -= compTrace(Invfft->out, Invfft->out, n);

    // tr{U2'L1L1'U2}
    vecConv(Invfft->in, L1fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace -= compTrace(Invfft->out, Invfft->out, n);

    // tr{U2'L2L2'U2}
    vecConv(Invfft->in, L2fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace += compTrace(Invfft->out, Invfft->out, n);

    // trace 
    trace /= Gs->Phi[0];
    trace /= acf2[0];
    return trace;
    //return value stored in trace;
}


// trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
inline double Toeplitz::traceDeriv(double* acf2, double* acf3){
    bool acf2_is_0 = zeroTest(acf2, n);    
    bool acf3_is_0 = zeroTest(acf3, n);

    if(acf2_is_0 | acf3_is_0){
        return 0.0;
    }

    if(!hasSolve){
        flagSolve();
    }

    double trace2;

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
    for(int ii = 0; ii < n; ++ii)
    {
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

    // 2.1, tr(L1(phi2) * L1(phi)' * L1'(acf3) * L1(acf3)) = compTrace(L1(acf3) * L1(phi2), 
    //      L1(acf3) * L1(phi)) = temp
    std::copy(acf3, acf3 + n, U1fft->in);
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

    trace2 -= 2 * compTrace(Invfft->out, temVec, n) / acf3[0];

    // 2.2, tr(L1(phi2) * L1(phi)' * L2'(acf3) * L2(acf3)) = compTrace(L2(acf3) * L1(phi2), 
    //      L2(acf3) * L1(phi)) = temp
    U2fft->in[0] = 0;
    std::copy(acf3 + 1, acf3 + n, U2fft->in + 1);
    U2fft->fft();

    vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    trace2 += 2 * compTrace(Invfft->out, temVec, n) / acf3[0];

    // 3, tr(_L2(phi2) * _L2(phi)' * Toeplitz_3) = tr(_L2(phi2) * _L2'(phi) * L1'(acf3)
    //    * L1(acf3)) / acf3[1] - tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) / acf3[1]

    // 3.1, tr(L2(phi2) * L2'(phi) * L1'(acf3) * L1(acf3)) = compTrace(L1(acf3) * _L2(phi2), 
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
    
    trace2 += 2 * compTrace(Invfft->out, temVec, n) / acf3[0];

    // 3.2, tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) = compTrace(L2(acf3) * _L2(phi2), 
    //      L2(acf3) * _L2(phi)) = temp
    vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    trace2 -= 2 * compTrace(Invfft->out, temVec, n) / acf3[0];

    // tr = tr / phi[1]
    trace2 /= Gs->Phi[0];
    trace2 *= -1;
    return trace2;
}

# endif