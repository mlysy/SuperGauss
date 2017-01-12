#include "Toeplitz.h"

// R wapper class

class Toeplitz : private Toep 
{
private:
    int n_R;
    int d_R;

    double* acf;
    double* acf2;
    double* acf3;
    double** x;
    double* vec;
public:
	Toeplitz(int, int);
    Toeplitz(int);
	~Toeplitz();
    Rcpp::List dimCheck_R(); // wrapper for dimension check
	void acfInput_R(NumericVector); // wapper for acfInput
    Rcpp::NumericMatrix mult_R(NumericMatrix); // wrapper for mult
    Rcpp::NumericMatrix solve_R(NumericMatrix); // wrapper for solve
    double det_R(); // wrapper returns the determinant of Toeplitz matrix
	double traceProd_R(NumericVector); // traceProd
    double traceDerv_R(NumericVector, NumericVector); // traceDerv
};

Toeplitz::Toeplitz(int n_, int d_): Toep(n_, d_)
{
    n_R = n_;
    d_R = d_;
    // Phi_R = NumericVector(n_R);
    acf = new double[n_R];
    acf2 = new double[n_R];
    acf3 = new double[n_R];
    x = new double*[d_R];
    for(int ii = 0; ii < d_R; ++ii) {
        x[ii] = new double[n_R];
    }
}

Toeplitz::Toeplitz(int n_): Toep(n_, 1)
{
    n_R = n_;
    d_R = 1;
    // Phi_R = NumericVector(n_R);
    acf = new double[n_R];
    acf2 = new double[n_R];
    acf3 = new double[n_R];
    x = new double*[d_R];
    for(int ii = 0; ii < d_R; ++ii) {
        x[ii] = new double[n_R];
    }
}

Toeplitz::~Toeplitz()
{
    delete[] acf;
    delete[] acf2;
    delete[] acf3;
    for(int ii = 0; ii < d_R; ++ii) {
        delete x[ii];
    }
    delete[] x;
}

Rcpp::List Toeplitz::dimCheck_R(){
    return List::create(_["N"] = n_R, _["d"] = d_R);
}

void Toeplitz::acfInput_R(NumericVector acf_R){
    if(acf_R.size() != n_R){
        cout << "non-conformable arguments" << endl;
        return;
    }
    std::copy(acf_R.begin(), acf_R.end(), acf);
    acfInput(acf);
    return;
}

// we have a special case for these 2 functions
NumericMatrix Toeplitz::mult_R(NumericMatrix x_R){
    if((x_R.nrow() != n_R) | (x_R.ncol() != d_R)) {
        if(d_R != 1){
            cout << "non-conformable arguments" << endl;
            return NumericMatrix();
        } 
        else{
            int d_xR = x_R.ncol();
            NumericMatrix Mult_R(n_R, d_xR);
            for(int ii = 0; ii < d_xR; ii++){
                std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x[0]);
                mult(x);
                std::copy(Mult[0], Mult[0] + n_R, Mult_R.column(ii).begin());
            }
            return Mult_R;
        }
    }
    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ++ii) {
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x[ii]);
    }
    mult(x);
    for(int ii = 0; ii < d_R; ++ii) {
        std::copy(Mult[ii], Mult[ii] + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

NumericMatrix Toeplitz::solve_R(NumericMatrix x_R){
    if((x_R.nrow() != n_R) | (x_R.ncol() != d_R)) {
        if(d_R != 1){
            cout << "non-conformable arguments" << endl;
            return NumericMatrix();
        } 
        else{
            int d_xR = x_R.ncol();
            NumericMatrix Mult_R(n_R, d_xR);
            for(int ii = 0; ii < d_xR; ii++){
                std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x[0]);
                solve(x);
                std::copy(Mult[0], Mult[0] + n_R, Mult_R.column(ii).begin());
            }
            return Mult_R;
        }
    }
    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ++ii) {
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x[ii]);
    }
    solve(x);
    for(int ii = 0; ii < d_R; ++ii) {
        std::copy(Mult[ii], Mult[ii] + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

// return the determinant
double Toeplitz::det_R() {
    double Det_R;
    detCheck();
    Det_R = det;
    return Det_R;
}


// return the TraceProd result
double Toeplitz::traceProd_R(NumericVector acf2_R)
{
    if(acf2_R.size() != n_R)
    {
        cout << "non-conformable acf2" << endl;
        return 0.0;
    }
    std::copy(acf2_R.begin(), acf2_R.end(), acf2);
    traceProd(acf2);
    double Trace_R;
    Trace_R = trace;
    return Trace_R;
}

double Toeplitz::traceDerv_R(NumericVector acf2_R, NumericVector acf3_R){
    if(!(acf2_R.size() == n_R & acf3_R.size() == n_R))
    {
        cout << "non-conformable acf2" << endl;
        return 0.0;
    }
    std::copy(acf2_R.begin(), acf2_R.end(), acf2);
    std::copy(acf3_R.begin(), acf3_R.end(), acf3);
    traceDerv(acf2, acf3);
    double Trace2_R;
    Trace2_R = trace2;
    return Trace2_R;
}


// wrapper

RCPP_MODULE(Toeplitz)
{
    using namespace Rcpp;
    class_<Toeplitz>("Toeplitz")
    .constructor<int, int>()
    .constructor<int>()
    .method("dimCheck", &Toeplitz::dimCheck_R)
    .method("acfInput", &Toeplitz::acfInput_R)
    .method("det", &Toeplitz::det_R)
    .method("solve", &Toeplitz::solve_R)
    .method("mult", &Toeplitz::mult_R)
    .method("traceprod", &Toeplitz::traceProd_R)
    .method("tracederv", &Toeplitz::traceDerv_R)
    ;
}


