#include "mex.h"
#include <cmath>
#include <complex>

/*
Inputs:
t - time vector of Nsteps length
xy - Nsteps x D trajectory
wv - wavevectors - D x K matrix
domain - width of the domain in each state variable ( 1 x D )
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

  const size_t D = mxGetM(p_wv);
  const size_t K = mxGetN(p_wv);
  const size_t Nsteps = mxGetNumberOfElements( p_t );

  if ( mxGetM(p_xy) != Nsteps )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "Number of rows in xy and elements in t"
                      " have to be the same");

  if ( mxGetN(p_xy) != D )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "Number of states in xy and rows of wv"
                      " have to be the same");

  if ( mxGetNumberOfElements(p_domain) != D )
    mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:input",
                      "Number of dimensions of state needs "
                      "to be consistent between xy, and wv.");

  double *t = mxGetPr( p_t );
  double *xy = mxGetPr( p_xy );
  double *wv = mxGetPr( p_wv );
  double *domain = mxGetPr(p_domain);

  for ( size_t d = 0; d < D; d++ )
    if (domain[d] <= 0)
      mexErrMsgIdAndTxt("DIFFEQ:computeAverages_helper:domain",
                        "'domain' argument has to contain "
                        "two positive numbers.");

  plhs[0] = mxCreateDoubleMatrix(K, 1, mxCOMPLEX );
  double *avg_r = mxGetPr( plhs[0] );
  double *avg_i = mxGetPi( plhs[0] );

  // Loop through observables and use trapezoidal rule to compute
  // their averages

  // loop through observables
  for ( size_t k = 0; k < K; k++ ) {

    std::complex<double> val = 0.;

    // loop through integration intervals
    double tstep;
    for ( size_t n = 0; n < (Nsteps-1); n++ ) {

      tstep = (t[n+1] - t[n]);
      
      // trapezoidal integral
      double argleft = 0.;
      double argright = 0.;
      for (size_t d = 0; d < D; d++) {
        argleft += wv[D*k + d] * xy[n + d*Nsteps] / domain[d];
        argright += wv[D*k + d] * xy[(n+1) + d*Nsteps] / domain[d];
      }
      val += tstep*exp( std::complex<double>(0., 2*M_PI*argleft ) )/2.;
      val += tstep*exp( std::complex<double>(0., 2*M_PI*argright ) )/2.;
    }

    // average by the length of time
    val /= (t[Nsteps-1]-t[0]); 

    avg_r[k] = real(val);
    avg_i[k] = imag(val);

  }

}
