#include "mex.h"
#include <cmath>
#include <complex>

/*
Inputs:

t - time vector (Nsteps x 1)
xy - trajectory (Nsteps x 2)
wv - wave-vectors -- each column is a pair of wavenumbers (K x 2)
domain - width and height of the domain ( 2 x 1 )
*/
#define p_t (prhs[0])
#define p_xy (prhs[1])
#define p_wv (prhs[2])
#define p_domain (prhs[3])

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {

  if (nrhs < 4) 
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "5 arguments required.");
  
  size_t Nsteps = mxGetNumberOfElements( p_t );
  double *t = mxGetPr( p_t );

  if ( mxGetM(p_xy) != Nsteps )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "Number of rows in xy and elements in t have to be the same");

  if ( mxGetN(p_xy) != 2 )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "xy has to have 2 columns.");
  
  if ( mxGetM(p_wv) != 2 )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "wv has to have 2 rows.");
  

  double *xy = mxGetPr( p_xy );
  double *wv = mxGetPr( p_wv );
  size_t K = mxGetN(p_wv);

  double *domain = mxGetPr(p_domain);

  if ( mxGetNumberOfElements(p_domain) != 2 || 
       domain[0] <= 0 || domain[1] <= 0 )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:domain",
                      "'domain' argument has to contain two positive numbers.");

  plhs[0] = mxCreateDoubleMatrix(K, 1, mxCOMPLEX );
  double *avg_r = mxGetPr( plhs[0] );
  double *avg_i = mxGetPi( plhs[0] );

  // loop through observables
  for ( size_t k = 0; k < K; k++ ) {

    std::complex<double> val = 0.;
    double wv_x = wv[ 2*k ] / domain[0];
    double wv_y = wv[ 2*k + 1 ] / domain[1];

    // loop through integration intervals
    double tstep;
    for ( size_t n = 0; n < (Nsteps-1); n++ ) {

      tstep = (t[n+1] - t[n]);
      
      // trapezoidal integral
      val += tstep*exp( std::complex<double>(0., 2*M_PI*(wv_x*xy[ n ] + 
                                                         wv_y*xy[ n + Nsteps]) ) )/2.;
      val += tstep*exp( std::complex<double>(0., 2*M_PI*(wv_x*xy[(n+1)]+ 
                                                         wv_y*xy[(n+1)+Nsteps]) ) )/2.;

    }

    // average by the length of time
    val /= (t[Nsteps-1]-t[0]); 

    avg_r[k] = real(val);
    avg_i[k] = imag(val);

  }

}
