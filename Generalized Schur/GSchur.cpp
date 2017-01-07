// implementing Generalized Schur Algorithm
// for input vector acf
// computing Toeplitz(acf)^-1, expressed in the form of xi, eta and gamma

#include "GSchur.h"

//defining fundamential functions

void CompMult(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n)
{
	for (int ii = 0; ii < n; ++ii)
	{
		y[ii][0] = alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] = alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

void CompMult_Minus(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n)
{
	for (int ii = 0; ii < n; ++ii)
	{
		y[ii][0] -= alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] -= alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

void CompMult_Plus(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n)
{
	for (int ii = 0; ii < n; ++ii)
	{
		y[ii][0] += alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] += alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

vector<int> intoBit(int x, int base = 1)
{
	vector<int> s;
	int x1 = x / base;
	int dif = x - x1 * base;
	int M = base;
	do
	{
		if (x1 & 1)
		{
			s.push_back(M);
		}
		M <<= 1;
	} while ((x1 >>= 1)>0);
	reverse(s.begin(), s.end());
	if (dif) 
	{
		s.push_back(dif);
	}
	return s;
}


//--------------------------------------------------------------------------------------------------
// constructor and destructor of class GSchur
GSchur::GSchur(int M)
{
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


GSchur::~GSchur()
{
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

//---------------------------------------------------------------------------------------
InverseToeplitz::InverseToeplitz(int n_, int base_)
{
	n = n_;
	base = base_;
	alpha = new double[n - 1];
	beta = new double[n - 1];
	s = intoBit(n - 1, base);
	Phi = new double[n];
	int jj = base;
	gs = new GSchur*[(int)log2(ceil((double)s[0]/base)) + 1];
	gs[0] = new GSchur(jj<<1);
	for (int ii = 0; ii < (int)log2(ceil((double)s[0]/base)); ++ii)
	{
		gs[ii + 1] = new GSchur(jj<<1);
		jj <<= 1;
	}
	gsM = new GSchur*[s.size()];
	for(int ii = 0; ii < (int)s.size(); ++ii)
	{
		gsM[ii] = new GSchur(s[ii]<<1);
	}
}

InverseToeplitz::~InverseToeplitz()
{
    delete[] alpha;
	delete[] beta;
	delete[] Phi;

	for (int ii = 0; ii < (int)s.size(); ++ii)
	{
		delete gsM[ii];
	}
	delete[] gsM;

	for (int ii = 0; ii <= (int)log2(ceil((double)s[0]/base)); ++ii)
	{
		delete gs[ii];
	}
	delete[] gs;	
}

//---------------------------------------------------------------------------------------
// helper functions applied in GSchur and GSchurMerge

// generate the alpham_IFFT and betam_IFFT used next, the size is M
// input of this function is gsr->alpha/beta/eta/xi_FFT->in
void alpha_beta(GSchur* gsr, int M)
{
	gsr->alpha_FFT->fft();
	gsr->beta_FFT->fft();
	gsr->eta_FFT->fft();
	gsr->xi_FFT->fft();

	for (int ii = 0; ii < M; ++ii)
	{
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

	CompMult(gsr->alpham_IFFT->in, gsr->alpha_FFT->out, gsr->eta_FFT->out, 2 * (M / 2 + 1));
	CompMult_Minus(gsr->alpham_IFFT->in, gsr->xi_FFT->out, gsr->beta_FFT->out, 2 * (M / 2 + 1));

	CompMult(gsr->betam_IFFT->in, gsr->beta_FFT->out, gsr->eta_t, 2 * (M / 2 + 1));
	CompMult_Minus(gsr->betam_IFFT->in, gsr->xi_t, gsr->alpha_FFT->out, 2 * (M / 2 + 1));
	gsr->alpham_IFFT->Ifft();
	gsr->betam_IFFT->Ifft();
}

// merging the pervious xi/eta and current xi/eta
// input is gsr->eta/xi_t, gsr->xi/eta_FFT->out, gsr->xi/eta_m_FFT->in, output is gsr->xi/eta_IFFT->out
void Merge(GSchur* gsr, int M)
{
	gsr->xi_m_FFT->fft();
	gsr->eta_m_FFT->fft();
	CompMult(gsr->xi_IFFT->in, gsr->eta_t, gsr->xi_m_FFT->out, 2 * (M / 2 + 1));
	CompMult_Plus(gsr->xi_IFFT->in, gsr->xi_FFT->out, gsr->eta_m_FFT->out, 2 * (M / 2 + 1));
		
	CompMult(gsr->eta_IFFT->in, gsr->xi_t, gsr->xi_m_FFT->out, 2 * (M / 2 + 1));
	CompMult_Plus(gsr->eta_IFFT->in, gsr->eta_FFT->out, gsr->eta_m_FFT->out, 2 * (M / 2 + 1));

	gsr->xi_IFFT->Ifft();
	gsr->eta_IFFT->Ifft();
}

//--------------------------------------------------------------------------------------------------
// pschur function
// the final result is stored in gs->gamma, ->eta_IFFT->out, ->xi_IFFT->out.
void InverseToeplitz::Pschur(double* alpha0, double* beta0, int size)
{
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

	for (int kk = 1; kk < size; ++kk)
	{
		alpha_1 = alpha0[kk];
		gs[0]->beta_FFT->in[kk] = beta0[kk];
		for (int jj = 1; jj <= kk; ++jj)
		{
			alpha_0 = alpha_1 - gs[0]->gamma[jj - 1] * gs[0]->beta_FFT->in[kk - jj + 1];
			gs[0]->beta_FFT->in[kk - jj + 1] -= gs[0]->gamma[jj - 1] * alpha_1;
			alpha_1 = alpha_0;
		}
		gs[0]->gamma[kk] = alpha_0 / gs[0]->beta_FFT->in[0];
		gs[0]->beta_FFT->in[0] *= 1 - gs[0]->gamma[kk] * gs[0]->gamma[kk];
		
		eta2[0] = 1.0;
        xi2[0] = alpha0[0]/beta0[0];		
		for (int jj = 1; jj <= kk; ++jj)
		{
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
void InverseToeplitz::Gschur(double* alpha0, double* beta0, int size, int layer)
{
	if(size <= base)
	{
		Pschur(alpha0, beta0, size);
		return;
	}
	
    Pschur(alpha0, beta0, base);
	
    int M = base;
	for (int m = 0; m < layer; ++m)
	{
		//copying information from layer[m] to [m+1]
		std::copy(alpha0, alpha0 + 2 * M, gs[m + 1]->alpha_FFT->in);
		std::copy(beta0, beta0 + 2 * M, gs[m + 1]->beta_FFT->in);
		std::copy(gs[m]->gamma, gs[m]->gamma + M, gs[m + 1]->gamma);
		std::copy(gs[m]->eta_IFFT->out, gs[m]->eta_IFFT->out + M, gs[m + 1]->eta_FFT->in);
		std::copy(gs[m]->xi_IFFT->out, gs[m]->xi_IFFT->out + M, gs[m + 1]->xi_FFT->in);

		//generating alpham and betam
		alpha_beta(gs[m + 1], M);
		//-----------------------------------------------------------------------------------------------
		// returns in layer[m]
		Gschur(gs[m + 1]->alpham_IFFT->out + M, gs[m + 1]->betam_IFFT->out + M, M, m);
		//-----------------------------------------------------------------------------------------------
		//copying information from layer[m] to [m+1]
		std::copy(gs[m]->gamma, gs[m]->gamma + M, gs[m + 1]->gamma + M);
		std::copy(gs[m]->xi_IFFT->out, gs[m]->xi_IFFT->out + M, gs[m + 1]->xi_m_FFT->in);
		std::copy(gs[m]->eta_IFFT->out, gs[m]->eta_IFFT->out + M, gs[m + 1]->eta_m_FFT->in);

		//merging xi and eta		
		Merge(gs[m + 1], M);

		M <<= 1;
	}
	return;
}

//------------------------------------------------------------------------------------------------------------

//gschur.merge function
void InverseToeplitz::GschurMerge()
{
	int layer = (int)log2(ceil((double)s[0]/base));
	Gschur(alpha, beta, s[0], layer);
	if(s.size() == 1)
	{
		std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[0], gsM[0]->eta_IFFT->out);
		std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[0], gsM[0]->xi_IFFT->out);
		std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
		return;
	}
	std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[0], gsM[0]->eta_FFT->in);
	std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[0], gsM[0]->xi_FFT->in);
	std::copy(gs[layer]->gamma, gs[layer]->gamma + s[0], gsM[0]->gamma);
	int M = s[0];
	for(int m = 0; m < (int)s.size()-1; ++m)
	{
		if(m == 0)
		{
			std::copy(alpha, alpha + n - 1, gsM[m]->alpha_FFT->in);
			std::copy(beta, beta + n - 1, gsM[m]->beta_FFT->in);
		}
		else
		{
			std::copy(gsM[m - 1]->alpham_IFFT->out + s[m - 1], gsM[m - 1]->alpham_IFFT->out + s[m - 1] + 2 * s[m], gsM[m]->alpha_FFT->in);
			std::copy(gsM[m - 1]->betam_IFFT->out + s[m - 1], gsM[m - 1]->betam_IFFT->out + s[m - 1] + 2 * s[m], gsM[m]->beta_FFT->in);
		}
		alpha_beta(gsM[m], s[m]);
		
		layer = (int)log2(ceil((double)s[m + 1]/base));
		Gschur(gsM[m]->alpham_IFFT->out + s[m], gsM[m]->betam_IFFT->out + s[m], s[m + 1], layer);
		
		std::copy(gs[layer]->eta_IFFT->out, gs[layer]->eta_IFFT->out + s[m + 1], gsM[m + 1]->eta_FFT->in);
		std::copy(gs[layer]->xi_IFFT->out, gs[layer]->xi_IFFT->out + s[m + 1], gsM[m + 1]->xi_FFT->in);
		std::copy(gs[layer]->gamma, gs[layer]->gamma + s[m + 1], gsM[0]->gamma + M);
		
		M += s[m + 1];
		
	}
	M = s[(int)s.size() - 1];
	for(int m = (int)s.size()-2; m >= 0; --m)
	{
		if(m == (int)s.size()-2)
		{
			std::copy(gsM[m + 1]->xi_FFT->in, gsM[m + 1]->xi_FFT->in + M, gsM[m]->xi_m_FFT->in);
			std::copy(gsM[m + 1]->eta_FFT->in, gsM[m + 1]->eta_FFT->in + M, gsM[m]->eta_m_FFT->in);
		}
		else
		{
			std::copy(gsM[m + 1]->xi_IFFT->out, gsM[m + 1]->xi_IFFT->out + M, gsM[m]->xi_m_FFT->in);
			std::copy(gsM[m + 1]->eta_IFFT->out, gsM[m + 1]->eta_IFFT->out + M, gsM[m]->eta_m_FFT->in);
		}
		Merge(gsM[m], s[m]);
		
		M += s[m];
	}
}

//output function
void InverseToeplitz::Inverse(double* acf)
{
	for(int ii = 0; ii < n - 1; ++ii)
	{
		alpha[ii] = -1 * acf[ii + 1];
        beta[ii] = acf[ii];
	}
	GschurMerge();

	double sigma2 = log(acf[0]);
	ldV = sigma2;
	for (int ii = 0; ii < n - 1; ++ii)
	{
		if(gsM[0]->gamma[ii] < 1){
			sigma2 += log(1 - gsM[0]->gamma[ii] * gsM[0]->gamma[ii]);
			ldV += sigma2;
		}
	}
	sigma2 = exp(sigma2);

	std::copy(gsM[0]->eta_IFFT->out, gsM[0]->eta_IFFT->out + n - 1, Phi);
	Phi[n-1] = 0;
	Phi[0] /= sigma2;
	for(int ii = 1; ii < n; ++ii)
	{
		Phi[ii] += gsM[0]->xi_IFFT->out[ii-1];
		Phi[ii] /= sigma2;

	}
	return;
}