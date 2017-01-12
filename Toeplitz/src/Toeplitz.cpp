// implementing Inverse Product
// for input vector acf and matrix x
// computing T^{-1} * x, called as InverseProd

// InnerProduction
// basically doing x * T^-1 * x' and returns log|T|, the log determinant

#include "Toeplitz.h"

//----------------------------------------------------------------------------------------


// calculating the trace of L1 L2', lower trangular toeplitz composed of acf of size n
void TraceComp(double* acf1, double* acf2, int n, double& trace){
    trace = 0;
    for(int ii = 0; ii < n; ++ii){
        trace += (n - ii) * acf1[ii] * acf2[ii];
    }
    return;
}


Toep::Toep(int n_, int d_ = 1){
    n = n_;
    d = d_;

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
    Mult = new double*[d];

    for(int ii = 0; ii < d; ++ii){
        Mult[ii] = new double[n];
    }

    hasAcf = FALSE;

    phi2 = new double[n];
    temVec = new double[n];
}
 
// requires edition
Toep::~Toep(){
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

    for(int ii = 0; ii < d; ++ii){
        delete Mult[ii];
    }
    delete[] Mult;
    
    delete[] phi2;
    delete[] temVec;
}

//-------------------------------------------------------------------------------------------------
void Toep::acfInput(double* acfA){
    // copy the acfA into Toep::acf
    std::copy(acfA, acfA + n, acf);
    // turn the flag into False
    hasMult = FALSE;
    hasInv = FALSE;
    hasAcf = TRUE;
}

//-------------------------------------------------------------------------------------------------
void Toep::computeMult(){
    // preparation for Toeplitz times Matrix
    // fft(c(acf, acf[n-1:2]))
    std::copy(acf, acf + n, Toepfft->in);
    std::copy(acf + 1, acf + n, Toepfft->in + n + 1);
    std::reverse(Toepfft->in + n + 1, Toepfft->in + 2 * n);
    Toepfft->fft();
    hasMult = TRUE;
}

// Toeplitz matrix vector multiplication
void Toep::mult(double** x){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }
    
    if(!hasMult){
        computeMult();
    }

    for(int ii = 0; ii < d; ++ii){
        std::copy(x[ii], x[ii] + n, xfft->in);
        xfft->fft();
        CompMult(Invfft->in, Toepfft->out, xfft->out, 2*(n/2 + 1));
        Invfft->Ifft(); 
        std::copy(Invfft->out, Invfft->out + n, Mult[ii]);
    }
    // return value stored in Mult.
} 

//-------------------------------------------------------------------------------------------------
// compute part, generate required parts for TraceProd and InverseMult
void Toep::computeInv(){

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


void Toep::solve(double** x){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }

    if(!hasInv){
        computeInv();
    }

	for(int ii = 0; ii < d; ++ii){
        std::copy(x[ii], x[ii] + n, xfft->in);
		xfft->fft();

        //L1 * L1' * x
		CompMult(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 
        std::copy(Invfft->out, Invfft->out + n, Lxfft->in); 
		Lxfft->fft();
		CompMult(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 

        // stored in Mult[ii]
        std::copy(Invfft->out, Invfft->out + n, Mult[ii]);

        // L2 * L2' * x
		CompMult(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 
        std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
		Lxfft->fft();
		CompMult(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 

        // 1/sigma2 (L1 L1'x - L2 L2'x) 
        for(int jj = 0; jj < n; ++jj)
        {
            Mult[ii][jj] -= Invfft->out[jj];
            Mult[ii][jj] /= Gs->Phi[0];
        }
    }
    // return value stored in Mult
}

void Toep::detCheck(){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }
    if(!hasInv){
        computeInv();
    }
    return;
}


// trace(Toeplitz^-1 * Toeplitz_2)
void Toep::traceProd(double* acf2){
    if(!hasAcf){
        cout << "Please input acf" << endl;
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
void Toep::traceDerv(double* acf2, double* acf3){
    if(!hasAcf){
        cout << "Please input acf" << endl;
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
        // phi2[ii] *= -1;
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
}

