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

    // return of mult and solve
	double* Mult;

	// log.deternimant
	double det;

	// return of TraceProd
	double trace;

	// return of TraceDerv
	double trace2;

public:		
	Toeplitz(int);
	~Toeplitz();

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

    // functions to show acf, Mult, det, trace, trace2
    // TODO
};
//------------------------------------------------------

// calculating the trace of L1 L2', lower trangular toeplitz composed of acf of size n
inline void TraceComp(double* acf1, double* acf2, int n, double& trace){
    trace = 0;
    for(int ii = 0; ii < n; ++ii){
        trace += (n - ii) * acf1[ii] * acf2[ii];
    }
    return;
}

inline void TestZero(double* acf, int n, bool& acf_is_0){
    acf_is_0 = true;
    for(int ii = 0; ii < n; ++ii){
        if(acf[ii]){
            acf_is_0 = false;
            break;
        }
    }
    return;        
}

inline Toeplitz::Toeplitz(int n_){
    n = n_;
    
	Gs = new InverseToeplitz(n, 64); // default base is choosen to be 64 ,which is the most efficient, according to numerical test
	
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
    Mult = new double[n];

    hasAcf = FALSE;

    phi2 = new double[n];
    temVec = new double[n];
}
 
// requires edition
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

    delete[] Mult;
    
    delete[] phi2;
    delete[] temVec;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::acfInput(double* acfA){

    TestZero(acfA, n, acf_is_0);

    if(!acf_is_0){
        std::copy(acfA, acfA + n, acf);
    }

    hasMult = FALSE;
    hasInv = FALSE;
    hasAcf = TRUE;
}

//-------------------------------------------------------------------------------------------------
inline void Toeplitz::computeMult(){
    std::copy(acf, acf + n, Toepfft->in);
    std::copy(acf + 1, acf + n, Toepfft->in + n + 1);
    std::reverse(Toepfft->in + n + 1, Toepfft->in + 2 * n);
    Toepfft->fft();
    hasMult = TRUE;
}

// Toeplitz matrix vector multiplication
inline void Toeplitz::mult(double* x){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(acf_is_0){
        std::fill(Mult, Mult + n, 0);
        return;
    }

    if(!hasMult){
        computeMult();
    }

    std::copy(x, x + n, xfft->in);
    xfft->fft();
    CompMult(Invfft->in, Toepfft->out, xfft->out, 2*(n/2 + 1));
    Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, Mult);
} 

//-------------------------------------------------------------------------------------------------
// compute part, generate required parts for TraceProd and InverseMult
inline void Toeplitz::computeInv(){

    Gs->Inverse(acf);

    det = Gs->ldV;

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

    hasInv = TRUE;
}


inline void Toeplitz::solve(double* x){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(acf_is_0){
        cout << "illegal to inverse a zero matrix" << endl;
        return;
    }

    if(!hasInv){
        computeInv();
    }

    std::copy(x, x + n, xfft->in);
	xfft->fft();

    //L1 * L1' * x
	CompMult(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in); 
	Lxfft->fft();
	CompMult(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 

    // stored in Mult
    std::copy(Invfft->out, Invfft->out + n, Mult);

    // L2 * L2' * x
	CompMult(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
	Lxfft->fft();
	CompMult(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 

    // 1/sigma2 (L1 L1'x - L2 L2'x) 
    for(int ii = 0; ii < n; ++ii)
    {
        Mult[ii] -= Invfft->out[ii];
        Mult[ii] /= Gs->Phi[0];
    }
}

inline void Toeplitz::detCheck(){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(acf_is_0){
        cout << "determinant of a zero matrix is 0" << endl;
        return;
    }

    if(!hasInv){
        computeInv();
    }
    return;
}


// trace(Toeplitz^-1 * Toeplitz_2)
inline void Toeplitz::traceProd(double* acf2){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(acf_is_0){
        cout << "illegal to inverse a zero matrix" << endl;
        return;
    }

    bool acf2_is_0;
    TestZero(acf2, n, acf2_is_0);
    if(acf2_is_0){
        trace = 0;
        return;
    }

    if(!hasInv){
        computeInv();
    }
    
    double temp;
    std::copy(acf2, acf2 + n, U1fft->in);
    U1fft->fft(); //U1
    std::copy(acf2 + 1, acf2 + n, U2fft->in + 1);
    U2fft->fft(); //U2

    // tr{U1'L1L1'U1}
    CompMult(Invfft->in, L1fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    TraceComp(Invfft->out, Invfft->out, n, temp);
    trace = temp;

    // tr{U1'L2L2'U1}
    CompMult(Invfft->in, L2fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    TraceComp(Invfft->out, Invfft->out, n, temp);
    trace -= temp; 

    // tr{U2'L1L1'U2}
    CompMult(Invfft->in, L1fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    TraceComp(Invfft->out, Invfft->out, n, temp);
    trace -= temp; 

    // tr{U2'L2L2'U2}
    CompMult(Invfft->in, L2fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    TraceComp(Invfft->out, Invfft->out, n, temp);
    trace += temp;
    // trace 
    trace /= Gs->Phi[0];
    trace /= acf2[0];
    //return value stored in trace;
}


// trace(Toeplitz^-1 * Toeplitz_2 * Toeplitz^-1 * Toeplitz_3)
inline void Toeplitz::traceDerv(double* acf2, double* acf3){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(acf_is_0){
        cout << "illegal to inverse a zero matrix" << endl;
        return;
    }

    bool acf2_is_0, acf3_is_0;
    TestZero(acf2, n, acf2_is_0);    
    TestZero(acf3, n, acf3_is_0);
    if(acf2_is_0 | acf3_is_0){
        trace2 = 0;
        return;
    }

    if(!hasInv){
        computeInv();
    }
    trace2 = 0;
    double temp;

    // phi2 = - Toeplitz^-1 * Toeplitz_2 * phi ------------------------------
    // phi2 = Toeplitz_2 * phi
    std::copy(acf2, acf2 + n, Lxfft->in);
    std::copy(acf2 + 1, acf2 + n, Lxfft->in + n + 1);
    std::reverse(Lxfft->in + n + 1, Lxfft->in + 2 * n);
    Lxfft->fft();
    std::fill(Lxfft->in + n, Lxfft->in + 2 * n, 0); // (readjust back to 0 assignment for latter half)
    std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
    xfft->fft();
    CompMult(Invfft->in, Lxfft->out, xfft->out, 2*(n/2 + 1));
    Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, phi2);
    
    // phi2 = - Toeplitz^-1 * phi2
    std::copy(phi2, phi2 + n, xfft->in);
	xfft->fft();

    //L1 * L1' * x
	CompMult(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in); 
	Lxfft->fft();
	CompMult(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 

    std::copy(Invfft->out, Invfft->out + n, phi2);

    //L2 * L2' * x
	CompMult(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
	Invfft->Ifft(); 
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
	Lxfft->fft();
	CompMult(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
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
    traceProd(acf3);
    trace2 = phi2[0] * trace;

    // 2, tr(L1(-phi2) * L1(phi)' * Toeplitz_3) = -tr(L1(phi2) * L1(phi)' * L1'(acf3)
    //    * L1(acf3)) + tr(L1(phi2) * L1(phi)' * L2'(acf3) * L2(acf3))

    // 2.1, tr(L1(phi2) * L1(phi)' * L1'(acf3) * L1(acf3)) = TraceComp(L1(acf3) * L1(phi2), 
    //      L1(acf3) * L1(phi)) = temp
    std::copy(acf3, acf3 + n, U1fft->in);
    U1fft->fft();

    std::copy(phi2, phi2 + n, Lxfft->in);
    Lxfft->fft();
    CompMult(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
    xfft->fft();
    CompMult(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    TraceComp(Invfft->out, temVec, n, temp);

    trace2 -= 2 * temp / acf3[0];

    // 2.2, tr(L1(phi2) * L1(phi)' * L2'(acf3) * L2(acf3)) = TraceComp(L2(acf3) * L1(phi2), 
    //      L2(acf3) * L1(phi)) = temp
    U2fft->in[0] = 0;
    std::copy(acf3 + 1, acf3 + n, U2fft->in + 1);
    U2fft->fft();

    CompMult(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    CompMult(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    TraceComp(Invfft->out, temVec, n, temp);

    trace2 += 2 * temp / acf3[0];

    // 3, tr(_L2(phi2) * _L2(phi)' * Toeplitz_3) = tr(_L2(phi2) * _L2'(phi) * L1'(acf3)
    //    * L1(acf3)) / acf3[1] - tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) / acf3[1]

    // 3.1, tr(L2(phi2) * L2'(phi) * L1'(acf3) * L1(acf3)) = TraceComp(L1(acf3) * _L2(phi2), 
    //      L1(acf3) * _L2(phi)) = temp
    Lxfft->in[0] = 0;
    std::copy(phi2 + 1, phi2 + n, Lxfft->in + 1);
    std::reverse(Lxfft->in + 1, Lxfft->in + n);
    Lxfft->fft();
    CompMult(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    xfft->in[0] = 0;
    std::copy(Gs->Phi + 1, Gs->Phi + n, xfft->in + 1);
    std::reverse(xfft->in + 1, xfft->in + n);
    xfft->fft();
    CompMult(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    
    TraceComp(Invfft->out, temVec, n, temp);

    trace2 += 2 * temp / acf3[0];

    // 3.2, tr(_L2(phi2) * _L2(phi)' * L2'(acf3) * L2(acf3)) = TraceComp(L2(acf3) * _L2(phi2), 
    //      L2(acf3) * _L2(phi)) = temp
    CompMult(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    std::copy(Invfft->out, Invfft->out + n, temVec);

    CompMult(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    TraceComp(Invfft->out, temVec, n, temp);

    trace2 -= 2 * temp / acf3[0];

    // tr = tr / phi[1]
    trace2 /= Gs->Phi[0];
    trace2 *= -1;
}

# endif