#ifndef OVERLAP_UTILS_H
#define OVERLAP_UTILS_H

#include "iset.h"

void check_indices (const Rcpp::IntegerVector&, const Rcpp::IntegerVector&, const Rcpp::IntegerVector&, int);

void set_mode_values(Rcpp::IntegerVector, int&, int&);

#endif
