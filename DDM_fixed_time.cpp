// in r: install rcpp library and ZigguratRCPP library, then source this script and the function becomes available (source('DDM_KSfitting'))
#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;
// [[Rcpp::export]]
NumericMatrix DDM_fixed_time(double v, int ntrials, double s, double dt, double time) {

  // initialize output
  NumericMatrix DATA(ntrials,2);

  // loop over trials
  for (int i = 0; i < ntrials; i++) {

    // initalize variables
    double evidence = 0;
    double t = 0;

    
    //Post-decisional processing
    for (int j = 0; j < time/dt; j++){
      t = t + dt;
      evidence = evidence + v * dt + s * sqrt(dt) * zigg.norm();
    }
    DATA(i,0) = evidence;
    DATA(i,1) = time;
  }

  return DATA; //evidence2, rt2
}
