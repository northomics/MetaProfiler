#include "deoptim.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

template<class VectorType1, class VectorType2>
VectorType2 three_exponential_function(const VectorType1 & x, const VectorType2 & xr) {
  double u = ((x[2] + x[0] + x[3]) - sqrt(std::pow(x[2] + x[0] + x[3], 2) - (4 * x[0]*x[3]))) / 2;
  double v = ((x[2] + x[0] + x[3]) + sqrt(std::pow(x[2] + x[0] + x[3], 2) - (4 * x[0]*x[3]))) / 2;
  double yu = x[0] * x[1]*(u - x[3]) / ((u - v)*(u - x[1])*u);
  double yv = x[0] * x[1]*(v - x[3]) / ((v - u)*(v - x[1])*v);
  double ykbi = x[0] * (x[1] - x[3]) / ((u - x[1])*(v - x[1]));
  VectorType2 value(xr.size());
  for(size_t i = 0; i < xr.size(); ++i) {
    value[i] = (1 + yu * exp(-u * xr[i]) + yv * exp(-v * xr[i]) + ykbi * exp(-x[1] * xr[i]))*x[4];
  }
  return value;
}

template<class VectorType1, class VectorType2>
VectorType2 two_exponential_function(const VectorType1 & x, const VectorType2 & xr) {
  double v = x[0];
  double yv = x[1] / (x[0] - x[1]);
  double ykbi = x[0] / (x[1] - x[0]);
  VectorType2 value(xr.size());
  for(size_t i = 0; i < xr.size(); ++i) {
    value[i] = (1 + yv * exp(-v * xr[i]) + ykbi * exp(-x[1] * xr[i]))*x[2];
  }
  return value;
}

template<class VectorType1, class VectorType2>
VectorType2 one_exponential_function(const VectorType1 & x, const VectorType2 & xr) {
  VectorType2 value(xr.size());
  for(size_t i = 0; i < xr.size(); ++i) {
    value[i] = x[1]*(1 - exp(-x[0] * xr[i]));
  }
  return value;
}


template<int equation>
double compute_score(const Rcpp::NumericVector & x, const arma::colvec & xr,
                     const arma::colvec & yr) {
  arma::colvec predicted;
  switch(equation) {
  case 1 :
    predicted = one_exponential_function(x, xr);
    break;
  case 2 :
    predicted = two_exponential_function(x, xr);
    break;
  case 3 :
    predicted = three_exponential_function(x, xr);
    break;
  default :
    Rcpp::stop("Not a valid equation number");
  }
  arma::colvec ri = arma::abs(predicted - yr);
  double score = 0;
  for(auto r : ri) {
    if(r > 5) {
      score += 5. * (r - 0.5 * 5.);
    } else {
      score +=  0.5 * std::pow(r, 2);
    }
  }
  return score;
}

// [[Rcpp::export]]
arma::mat curve_fitting_c(const std::vector<arma::colvec> & xr,
                          const std::vector<arma::colvec> & yr,
                          const arma::colvec & minbound,                  // user-defined lower bounds
                          const arma::colvec & maxbound,                  // user-defined upper bounds
                          const int equation,
                          const Rcpp::List & control,                     // parameters
                          bool verbose = true)
{
  arma::mat consts(minbound.size() + 1, xr.size());
  for(size_t i = 0; i < yr.size(); i++)
  {
    std::vector<bool> ok(yr[i].size());
    int count = 0;
    for(size_t j = 0; j < yr[i].size(); ++j) {
      if(!Rcpp::NumericVector::is_na(yr[i][j])) {
        ok[j] = true;
        ++count;
      } else {
        ok[j] = false;
      }
    }
    arma::colvec xr_new(count);
    arma::colvec yr_new(count);
    int k = 0;
    for(size_t j = 0; j < xr[i].size(); ++j) {
      if(ok[j]) {
        xr_new[k] = xr[i][j];
        yr_new[k] = yr[i][j];
        ++k;
      }
    }
    
    arma::colvec vec(minbound.size() + 1);
    switch(equation) {
    case 1 :
      DEoptim_impl(vec, minbound, maxbound, compute_score<1>, control, xr_new, yr_new);
      break;
    case 2 :
      DEoptim_impl(vec, minbound, maxbound, compute_score<2>, control, xr_new, yr_new);
      break;
    case 3 :
      DEoptim_impl(vec, minbound, maxbound, compute_score<3>, control, xr_new, yr_new);
      break;
    default :
      Rcpp::stop("Not a valid equation number");
    }
    vec[minbound.size()] = std::sqrt(vec[minbound.size()]/(xr_new.size() - 2));
    consts.col(i) = vec;
    if(verbose) {
      Rcpp::Rcout << i << "/" << xr.size() << "              \r";
    }
  }
  if(verbose) {
    Rcpp::Rcout << xr.size() << "/" << xr.size();
  }
  return consts;
}

// [[Rcpp::export]]
arma::mat curve_fitting_test_c(const std::vector<arma::colvec> & xr,
                               const std::vector<arma::colvec> & yr,
                               const arma::colvec & minbound,                  // user-defined lower bounds
                               const arma::colvec & maxbound,                  // user-defined upper bounds
                               const int equation,
                               const Rcpp::List & control,                     // parameters
                               bool verbose = true)
{
  std::vector<int> timepoints = Rcpp::as<std::vector<int>>(control["timepoints"]);            // value to reach
  arma::mat RMSE(timepoints.size(), xr.size());
  for(size_t i = 0; i < yr.size(); i++)
  {
    std::vector<bool> ok(yr[i].size());
    int count = 0;
    for(size_t j = 0; j < yr[i].size(); ++j) {
      if(!Rcpp::NumericVector::is_na(yr[i][j])) {
        ok[j] = true;
        ++count;
      } else {
        ok[j] = false;
      }
    }
    arma::colvec xr_no_NA(count);
    arma::colvec yr_no_NA(count);
    int k = 0;
    for(size_t j = 0; j < yr[i].size(); ++j) {
      if(ok[j]) {
        xr_no_NA[k] = xr[i][j];
        yr_no_NA[k] = yr[i][j];
        ++k;
      }
    }
    auto xr_it1 = xr_no_NA.begin();
    auto yr_it1 = yr_no_NA.begin();
    
    arma::colvec predicted(xr_no_NA.size());
    auto p_it = predicted.begin();
    
    
    while((xr_it1 != xr_no_NA.end()) & (yr_it1 != yr_no_NA.end())) {
      arma::colvec xr_new(xr_no_NA.size() - 1);
      arma::colvec yr_new(yr_no_NA.size() - 1);
      
      auto xr_it2 = xr_no_NA.begin();
      auto yr_it2 = yr_no_NA.begin();
      
      auto xr_new_it = xr_new.begin();
      auto yr_new_it = yr_new.begin();
      
      std::vector<double> xr_left_out(1);
      std::vector<double> yr_left_out(1);
      
      while((xr_it2 != xr_no_NA.end()) & (yr_it2 != yr_no_NA.end()) & (xr_new_it != xr_new.end()) & (yr_new_it != xr_new.end())) {
        if((xr_it1 != xr_it2) & (yr_it1 != yr_it2)) {
          *xr_new_it = *xr_it2;
          *yr_new_it = *yr_it2;
          xr_new_it++;
          yr_new_it++;
          xr_it2++;
          yr_it2++;
        } else {
          xr_left_out[0] = *xr_it2;
          yr_left_out[0] = *yr_it2;
          xr_it2++;
          yr_it2++;
        }
      }
      arma::colvec vec(minbound.size() + 1);
      switch(equation) {
      case 1 :
        DEoptim_impl(vec, minbound, maxbound, compute_score<1>, control, xr_new, yr_new);
        break;
      case 2 :
        DEoptim_impl(vec, minbound, maxbound, compute_score<2>, control, xr_new, yr_new);
        break;
      case 3 :
        DEoptim_impl(vec, minbound, maxbound, compute_score<3>, control, xr_new, yr_new);
        break;
      default :
        Rcpp::stop("Not a valid equation number");
      }
      arma::colvec value;
      switch(equation) {
      case 1 :
        value = one_exponential_function(vec, xr_left_out);
        break;
      case 2 :
        value = two_exponential_function(vec, xr_left_out);
        break;
      case 3 :
        value = three_exponential_function(vec, xr_left_out);
        break;
      default :
        Rcpp::stop("Not a valid equation number");
      }
      *p_it = std::pow(value[0] - yr_left_out[0], 2);
      p_it++;
      xr_it1++;
      yr_it1++;
    }
    for(size_t j = 0, l = 0; j < timepoints.size(); ++j) {
      if(timepoints[j] == xr_no_NA[l]) {
        RMSE.col(i).row(j) = predicted[l];
        ++l;
      }
    }
    if(verbose) {
      Rcpp::Rcout << i << "/" << xr.size() << "              \r";
    }
  }
  if(verbose) {
    Rcpp::Rcout << xr.size() << "/" << xr.size() << "\n";
  }
  return RMSE;
}




