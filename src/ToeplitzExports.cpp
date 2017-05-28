#include "Toeplitz.h"

// R wapper class

class Toeplitz_Cpp : private Toeplitz
{
private:
    int n_R;

    double* acf;
    double* acf2;
    double* acf3;
    double* xIn;
    double* yOut;
public:
  Toeplitz_Cpp(int);
  ~Toeplitz_Cpp();
  int dimCheck_R(); // wrapper for dimension check
  void acfInput_R(NumericVector); // wapper for acfInput
  Rcpp::NumericVector acfOutput_R(); // wrapper for acfOutput
  Rcpp::NumericMatrix mult_R(NumericMatrix); // wrapper for mult
  Rcpp::NumericMatrix mult_Vec(NumericVector);
  Rcpp::NumericMatrix solve_R(NumericMatrix); // wrapper for solve
  Rcpp::NumericMatrix solve_Vec(NumericVector);
  double det_R(); // wrapper returns the determinant of Toeplitz matrix
  double traceProd_R(NumericVector); // traceProd
  double traceDeriv_R(NumericVector, NumericVector); // traceDeriv
};

Toeplitz_Cpp::Toeplitz_Cpp(int n_): Toeplitz(n_)
{
    n_R = n_;
    acf = new double[n_R];
    acf2 = new double[n_R];
    acf3 = new double[n_R];
    xIn = new double[n_R];
    yOut = new double[n_R];
}

Toeplitz_Cpp::~Toeplitz_Cpp()
{
    delete[] acf;
    delete[] acf2;
    delete[] acf3;
    delete[] xIn;
    delete[] yOut;
}


int Toeplitz_Cpp::dimCheck_R(){
    return n_R;
}

void Toeplitz_Cpp::acfInput_R(NumericVector acf_R){
    if(acf_R.size() != n_R){
        cout << "non-conformable arguments" << endl;
        return;
    }
    std::copy(acf_R.begin(), acf_R.end(), acf);
    acfInput(acf);
    return;
}

NumericVector Toeplitz_Cpp::acfOutput_R(){
  NumericVector acf_R(n_R);
  std::copy(acf, acf + n_R, acf_R.begin());
  return acf_R;
}

NumericMatrix Toeplitz_Cpp::mult_R(NumericMatrix x_R){
    if(x_R.nrow() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }
    int d_R = x_R.ncol();

    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ii++){
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), xIn);
        multVec(yOut, xIn);
        std::copy(yOut, yOut + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

NumericMatrix Toeplitz_Cpp::mult_Vec(NumericVector x_R){
    if(x_R.size() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }

    NumericMatrix Mult_R(n_R, 1);
    std::copy(x_R.begin(), x_R.end(), xIn);
    multVec(yOut, xIn);
    std::copy(yOut, yOut + n_R, Mult_R.column(0).begin());
    return Mult_R;
}

NumericMatrix Toeplitz_Cpp::solve_R(NumericMatrix x_R){
    if(x_R.nrow() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }
    int d_R = x_R.ncol();
    NumericMatrix Mult_R(n_R, d_R);
    for(int ii = 0; ii < d_R; ii++){
        std::copy(x_R.column(ii).begin(), x_R.column(ii).end(), xIn);
        solveVec(yOut, xIn);
        std::copy(yOut, yOut + n_R, Mult_R.column(ii).begin());
    }
    return Mult_R;
}

NumericMatrix Toeplitz_Cpp::solve_Vec(NumericVector x_R){
    if(x_R.size() != n_R) {
        cout << "non-conformable arguments" << endl;
        return NumericMatrix();
    }
    NumericMatrix Mult_R(n_R, 1);
    std::copy(x_R.begin(), x_R.end(), xIn);
    solveVec(yOut, xIn);
    std::copy(yOut, yOut + n_R, Mult_R.column(0).begin());
    return Mult_R;
}

// return the determinant
double Toeplitz_Cpp::det_R() {
    double Det_R = logDet();
    return Det_R;
}

// return the TraceProd result
double Toeplitz_Cpp::traceProd_R(NumericVector acf2_R)
{
    if(acf2_R.size() != n_R)
    {
        cout << "non-conformable acf2" << endl;
        return NA_REAL;
    }
    std::copy(acf2_R.begin(), acf2_R.end(), acf2);
    double Trace_R = traceProd(acf2);
    return Trace_R;
}

double Toeplitz_Cpp::traceDeriv_R(NumericVector acf2_R, NumericVector acf3_R){
    if(!(acf2_R.size() == n_R & acf3_R.size() == n_R))
    {
        cout << "non-conformable acf2" << endl;
        return NA_REAL;
    }

    std::copy(acf2_R.begin(), acf2_R.end(), acf2);
    std::copy(acf3_R.begin(), acf3_R.end(), acf3);
    double Trace2_R = traceDeriv(acf2, acf3);
    return Trace2_R;
}

// wrapper
RCPP_MODULE(Toeplitz_Class)
{
    class_<Toeplitz_Cpp>("Toeplitz_Cpp")
    .constructor<int>()
      //.field("acf", &Toeplitz_Cpp::acf_Num)
    .method("DimCheck", &Toeplitz_Cpp::dimCheck_R)
    .method("setAcf", &Toeplitz_Cpp::acfInput_R)
    .method("getAcf", &Toeplitz_Cpp::acfOutput_R)
    .method("Det", &Toeplitz_Cpp::det_R)
    .method("Solve", &Toeplitz_Cpp::solve_R)
    .method("SolveVec", &Toeplitz_Cpp::solve_Vec)
    .method("Mult", &Toeplitz_Cpp::mult_R)
    .method("MultVec", &Toeplitz_Cpp::mult_Vec)
    .method("traceT2", &Toeplitz_Cpp::traceProd_R)
    .method("traceT4", &Toeplitz_Cpp::traceDeriv_R)
    ;
}


