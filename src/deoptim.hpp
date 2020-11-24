// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Port of DEoptim (2.0.7) to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2015  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004  
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#pragma once
#include "devol.hpp"      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// #include <google/profiler.h>

void DEoptim_impl(arma::vec & consts,
                  const arma::colvec & minbound,                  // user-defined lower bounds
                  const arma::colvec & maxbound,                  // user-defined upper bounds
                  std::function<double(const Rcpp::NumericVector &,const arma::colvec &,const arma::colvec &)> fnS,                                       // function to be optimized, either R or C++
                  const Rcpp::List & control,                     // parameters 
                  const arma::colvec & xr,
                  const arma::colvec & yr) {                                    // optional environment
  
#if RCPP_DEV_VERSION >= RcppDevVersion(0,12,17,1)
  Rcpp::SuspendRNGSynchronizationScope rngScope;
#endif
  
  double VTR           = Rcpp::as<double>(control["VTR"]);            // value to reach
  int i_strategy       = Rcpp::as<int>(control["strategy"]);          // chooses DE-strategy
  int i_itermax        = Rcpp::as<int>(control["itermax"]);           // Maximum number of generations
  long l_nfeval        = 0;                                           // nb of function evals (NOT passed in)
  int i_D              = Rcpp::as<int>(control["npar"]);              // Dimension of parameter vector
  int i_NP             = Rcpp::as<int>(control["NP"]);                // Number of population members
  int i_storepopfrom   = Rcpp::as<int>(control["storepopfrom"]) - 1;  // When to start storing populations 
  int i_storepopfreq   = Rcpp::as<int>(control["storepopfreq"]);      // How often to store populations 
  int i_specinitialpop = Rcpp::as<int>(control["specinitialpop"]);    // User-defined inital population 
  double f_weight      = Rcpp::as<double>(control["F"]);              // stepsize 
  double f_cross       = Rcpp::as<double>(control["CR"]);             // crossover probability 
  int i_bs_flag        = Rcpp::as<int>(control["bs"]);                // Best of parent and child 
  int i_trace          = Rcpp::as<int>(control["trace"]);             // Print progress? 
  double i_pPct        = Rcpp::as<double>(control["p"]);              // p to define the top 100p% best solutions 
  double d_c           = Rcpp::as<double>(control["c"]);              // c as a trigger of the JADE algorithm
  double d_reltol      = Rcpp::as<double>(control["reltol"]);         // tolerance for relative convergence test, default to be sqrt(.Machine$double.eps)
  int i_steptol        = Rcpp::as<double>(control["steptol"]);        // maximum of iteration after relative convergence test is passed, default to be itermax
  
  // as above, doing it in two steps is faster
  Rcpp::NumericMatrix initialpopm = Rcpp::as<Rcpp::NumericMatrix>(control["initialpop"]);
  arma::mat initpopm(initialpopm.begin(), initialpopm.rows(), initialpopm.cols(), false);
  
  arma::mat ta_popP(i_D, i_NP*2);                                     // Data structures for parameter vectors 
  arma::mat ta_oldP(i_D, i_NP);
  arma::mat ta_newP(i_D, i_NP);
  arma::colvec t_bestP(i_D); 
  
  arma::colvec ta_popC(i_NP*2);                                       // Data structures for obj. fun. values 
  arma::colvec ta_oldC(i_NP);
  arma::colvec ta_newC(i_NP);
  double t_bestC; 
  
  arma::colvec t_bestitP(i_D);
  arma::colvec t_tmpP(i_D); 
  
  int i_nstorepop = static_cast<int>(ceil(static_cast<double>((i_itermax - i_storepopfrom) / i_storepopfreq)));
  arma::mat d_pop(i_D, i_NP); 
  Rcpp::List d_storepop(i_nstorepop);
  arma::mat d_bestmemit(i_D, i_itermax);       
  arma::colvec d_bestvalit(i_itermax);     
  int i_iter = 0;
  
  // call actual Differential Evolution optimization given the parameters
  devol(VTR, f_weight, f_cross, i_bs_flag, minbound, maxbound, fnS, xr, yr, i_trace, i_strategy, i_D, i_NP, 
        i_itermax, initpopm, i_storepopfrom, i_storepopfreq, i_specinitialpop,
        ta_popP, ta_oldP, ta_newP, t_bestP, ta_popC, ta_oldC, ta_newC, t_bestC, t_bestitP, t_tmpP,
        d_pop, d_storepop, d_bestmemit, d_bestvalit, i_iter, i_pPct, d_c, l_nfeval,
        d_reltol, i_steptol);
  consts.rows(0,i_D - 1) = t_bestP;
  consts.row(i_D) = t_bestC;  
}