// implementing Inverse Product
// for input vector acf and matrix x
// computing T^{-1} * x, called as InverseProd

// InnerProduction
// basically doing x * T^-1 * x' and returns log|T|, the log determinant

#include "Toeplitz.h"

//----------------------------------------------------------------------------------------


// calculating the trace of L1 L1', lower trangular toeplitz composed of acf ofsize n
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
    xfft = new VectorFFT(2 * n);

    Invfft = new VectorIFFT(2 * n);

    acf = new double[n];
    Mult = new double*[d];
    for(int ii = 0; ii < d; ++ii){
        Mult[ii] = new double[n];
    }

    hasAcf = FALSE;
    // phiD = new double[n];
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
	delete Invfft;

    // delete U1fft;
    // delete U2fft;
    
    delete Toepfft;
    // delete xfft;

    delete[] acf;
    for(int ii = 0; ii < d; ++ii){
        delete Mult[ii];
    }
    delete[] Mult;
    // delete[] phiD;
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

        //L1' * x
		CompMult(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 
        std::copy(Invfft->out, Invfft->out + n, Lxfft->in); 
        
        //L1 * L1' * x
		Lxfft->fft();
		CompMult(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 

        // results stored in Mult[ii]
        std::copy(Invfft->out, Invfft->out + n, Mult[ii]);

        //here we finish L2' * x
		CompMult(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 
        std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
        //here we finish L2 * L2' * x
		Lxfft->fft();
		CompMult(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
		Invfft->Ifft(); 

        // omputing 1/sigma2 (L1 L1'x - L2 L2'x) 
        for(int jj = 0; jj < n; ++jj)
        {
            Mult[ii][jj] -= Invfft->out[jj];
            Mult[ii][jj] /= Gs->sigma2;
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


// trace part
void Toep::TraceProd(double* acf2){
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
    trace /= Gs->sigma2;
    trace /= acf2[0];
    //return value stored in trace;
}


/*
// New trace part, the TraceDerv function
void Toep::TraceDerv(double* acf2, double* acf3){
    if(!hasAcf){
        cout << "Please input acf" << endl;
        return;
    }
    if(!hasInv){
        computeInv();
    }
    trace = 0;
    // phiD = acf2 %*% phi * sigma2
    // phiD = acf^-1 %*% phiD
    // avaliable spaces are 
}
*/