#include "Toeplitz.h"

// R wapper class

class Toeplitz : private Toep 
{
private:
    int n_R;

    double* acf;
    double* acf2;
    double* acf3;
    double* x;
    double* vec;
public:
    Toeplitz(int);
	~Toeplitz();
    int dimCheck_R(); // wrapper for dimension check
	void acfInput_R(NumericVector); // wapper for acfInput
    Rcpp::NumericMatrix mult_R(NumericMatrix); // wrapper for mult
    Rcpp::NumericMatrix mult_Vec(NumericVector);
    Rcpp::NumericMatrix solve_R(NumericMatrix); // wrapper for solve
    Rcpp::NumericMatrix solve_Vec(NumericVector);
    double det_R(); // wrapper returns the determinant of Toeplitz matrix
	double traceProd_R(NumericVector); // traceProd
    double traceDerv_R(NumericVector, NumericVector); // traceDerv
};

Toeplitz::Toeplitz(int n_): Toep(n_)
{
    n_R = n_;
    acf = new double[n_R];
    acf2 = new double[n_R];
    acf3 = new double[n_R];
    x = new double[n_R];
}

Toeplitz::~Toeplitz()
{
    cout <<"rcpp destructor called" << endl;
    delete[] acf;
    delete[] acf2;
    delete[] acf3;
    delete[] x;
}

int Toeplitz::dimCheck_R(){
    return n_R;
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

NumericMatrix Toeplitz::mult_R(NumericMatrix x_R){
    if(x_R.nrow() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }
    int d_R = x_R.ncol();
    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ii++){
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x);
        mult(x);
        std::copy(Mult, Mult + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

NumericMatrix Toeplitz::mult_Vec(NumericVector x_R){
    if(x_R.size() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }
    NumericMatrix Mult_R(n_R, 1);
    std::copy(x_R.begin(), x_R.end(), x);
    mult(x);
    std::copy(Mult, Mult + n_R, Mult_R.column(0).begin());
    return Mult_R;
}

NumericMatrix Toeplitz::solve_R(NumericMatrix x_R){
    if(x_R.nrow() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    } 
    int d_R = x_R.ncol();
    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ii++){
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), x);
        solve(x);
        std::copy(Mult, Mult + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

NumericMatrix Toeplitz::solve_Vec(NumericVector x_R){
    if(x_R.size() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    } 
    NumericMatrix Mult_R(n_R, 1);
    std::copy(x_R.begin(), x_R.end(), x);
    solve(x);
    std::copy(Mult, Mult + n_R, Mult_R.column(0).begin());
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
    class_<Toeplitz>("Toeplitz")
    .constructor<int>()
    .method("DimCheck", &Toeplitz::dimCheck_R)
    .method("AcfInput", &Toeplitz::acfInput_R)
    .method("Det", &Toeplitz::det_R)
    .method("Solve", &Toeplitz::solve_R)
    .method("SolveVec", &Toeplitz::solve_Vec)
    .method("Mult", &Toeplitz::mult_R)
    .method("MultVec", &Toeplitz::mult_Vec)
    .method("TraceProd", &Toeplitz::traceProd_R)
    .method("TraceDeriv", &Toeplitz::traceDerv_R)
    ;
}


