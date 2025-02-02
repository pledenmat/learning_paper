// in r: install rcpp library, then source this script and the function becomes available
#include <Rcpp.h>
//#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double bin_data_double(double cj, int nbin) {
  return floor(cj * nbin)/nbin;
}

// [[Rcpp::export]]
double Mode(NumericVector x, bool narm = false) {
  if (narm) x = x[!is_na(x)];
  NumericVector ux = unique(x);
  double y = ux[which_max(table(match(x, ux)))];
  return y;
}

// [[Rcpp::export]]
NumericVector bin_data_vec(NumericVector cj, int nbin) {
  int n = cj.size();
  NumericVector temp_cj = clone(cj);
  for (int j = 0; j < n; j++){
    temp_cj[j] = floor(temp_cj[j] * nbin)/nbin;
  }
  return temp_cj;
}

// define a function to calculate the error based on the specified error_type
// [[Rcpp::export]]
double error(NumericVector y, NumericVector y_pred, String error_type){
  double err = 0;
  int n = y.length();
  if (error_type == "cross-entropy") {
    for (int i = 0; i < n; i++) {
      err -= y[i] * log(y_pred[i]);
    }
    err /= n;
  } else if (error_type == "mse") {
    for (int i = 0; i < n; i++) {
      err += pow(y[i] - y_pred[i], 2) / 2;
    }
    err /= n;
  } else {
    stop("Invalid error type");
  }
  return err;
}

// define a function to calculate the gradient of the error
// [[Rcpp::export]]
NumericVector grad(NumericVector x_j, double y_j, double pred_j, String error_type){
  NumericVector gradient(3);
  if (error_type == "cross-entropy") {
    gradient[0] = (y_j - pred_j) * x_j[0] / sqrt(x_j[2]);
    gradient[1] = (y_j - pred_j) * x_j[1] / sqrt(x_j[2]);
    gradient[2] = 0;
  } else if (error_type == "mse") {
    gradient[0] = pred_j * (1 - pred_j) * (y_j - pred_j) * x_j[0] / sqrt(x_j[2]);
    gradient[1] = pred_j * (1 - pred_j) * (y_j - pred_j) * x_j[1] / sqrt(x_j[2]);
    gradient[2] = 0;
  } else {
    stop("Invalid error type");
  }
  return gradient;
}

// Define helper functions
// [[Rcpp::export]]
double combine_input_1trial(NumericVector x, NumericVector w, bool binning = false,
                            int nbin = 6) {
  double z = w[2] * x[2] * (w[0] * x[0] + w[1] * x[1]);
  z =  1 / (1 + exp(-z)); // sigmoid activation function
  if (binning) {
    z = bin_data_double(z, nbin);
  }
  return z;
}

// [[Rcpp::export]]
NumericVector combine_input(NumericMatrix x, NumericVector w, bool binning = false,
                            int nbin = 6) {
  NumericVector z = w[2] * x(_,2) * (w[0] * x(_,0) + w[1] * x(_,1));
  z = 1 / (1 + exp(-z)); // sigmoid activation function
  z[z == 1] = 0.9999999; // Avoid log(0)
  z[z == 0] = 0.0000001;
  if (binning) {
    z = bin_data_vec(z,nbin);
  }
  return z;
}

// Define main function to train the model
// [[Rcpp::export]]
List train_model(NumericMatrix x, NumericVector w0, NumericVector y, NumericVector y_err, double eta = 0.1,
                 std::string error_type1 = "cross-entropy",
                 bool verbose = false, bool trace = false, bool binning = false, 
                 int nbin = 6.0, std::string cost = "per_trial",
                 NumericMatrix x_err = NumericMatrix(), int Nsim_error = 1000,
                 std::string error_type2 = "mse") {
  int ntrain = y.size();
  double err = 0;
  NumericVector w = clone(w0);
  NumericMatrix dat_iter(ntrain, 3);
  
  for (int j = 0; j < ntrain; j++) {
    if (!trace){
      NumericVector y_pred = combine_input(x_err(Range(j*Nsim_error,(j+1)*Nsim_error-1), _), w, binning = binning, nbin = nbin);
      double y_pred_mean = mean(y_pred);
      err += error(wrap(y_err[j]), wrap(y_pred_mean), error_type2);
    } 
    double pred_j = combine_input_1trial(x(j, _), w, binning = binning, nbin = nbin);
    // Update weights (time gradient is fixed to 0)
    NumericVector g = grad(x(j, _), y[j], pred_j, error_type1);
    w[0] += eta * g[0];
    w[1] += eta * g[1];
    w[2] += eta * g[2];
    
    if (verbose) {
      Rcout << "  -> updating data point " << j + 1 << " :\n";
      Rcout << "     -> alpha: " << w[0] << "\n";
      Rcout << "     -> beta: " << w[1] << "\n";
    }
    if (trace) {
      dat_iter(j, 0) = w[0];
      dat_iter(j, 1) = w[1];
    }
  }
  err/= ntrain;
  if (trace) {
    List results = List::create(Named("w") = w , _["trace"] = dat_iter);
    return results;
  } else {
    List results = List::create(Named("w") = w , _["err"] = err);
    return results;
  }
}

// Define main function to train the model
// [[Rcpp::export]]
List train_model_eta_sep(NumericMatrix x, NumericVector w0, NumericVector y, NumericVector y_err, double eta_a = 0.1,
                 double eta_b = 0.1, std::string error_type1 = "cross-entropy", 
                 bool verbose = false, bool trace = false, 
                 bool binning = false,int nbin = 6.0, std::string cost = "per_trial",
                 NumericMatrix x_err = NumericMatrix(), int Nsim_error = 1000,
                 std::string error_type2 = "mse") {
  int ntrain = y.size();
  double err = 0;
  NumericVector w = clone(w0);
  NumericMatrix dat_iter(ntrain, 3);
  
  for (int j = 0; j < ntrain; j++) {
    if (trace) {
      dat_iter(j, 0) = w[0];
      dat_iter(j, 1) = w[1];
    } else {
      NumericVector y_pred = combine_input(x_err(Range(j*Nsim_error,(j+1)*Nsim_error-1), _), w, binning = binning, nbin = nbin);
      double y_pred_mean = mean(y_pred);
      err += error(wrap(y_err[j]), wrap(y_pred_mean), error_type2);
    } 
    // Compute confidence prediction
    double pred_j = combine_input_1trial(x(j, _), w, binning = binning, nbin = nbin);
    // Update weights (time gradient is fixed to 0)
    NumericVector g = grad(x(j, _), y[j], pred_j, error_type1);
    w[0] += eta_a * g[0];
    w[1] += eta_b * g[1];
    w[2] += 0 * g[2];
    
    if (verbose) {
      Rcout << "  -> updating data point " << j + 1 << " :\n";
      Rcout << "     -> alpha: " << w[0] << "\n";
      Rcout << "     -> beta: " << w[1] << "\n";
    }
  }
  err/= ntrain;
  if (trace) {
    List results = List::create(Named("w") = w, _["trace"] = dat_iter);
    return results;
  } else {
    List results = List::create(Named("w") = w , _["err"] = err);
    return results;
  }
}