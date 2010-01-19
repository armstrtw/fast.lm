///////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009  Whit Armstrong                                    //
//                                                                       //
// This program is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published by  //
// the Free Software Foundation, either version 3 of the License, or     //
// (at your option) any later version.                                   //
//                                                                       //
// This program is distributed in the hope that it will be useful,       //
// but WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
// GNU General Public License for more details.                          //
//                                                                       //
// You should have received a copy of the GNU General Public License     //
// along with this program.  If not, see <http://www.gnu.org/licenses/>. //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <vector>

#include <armadillo>
#include <Rinternals.h>
#include "interface.hpp"
#include "index.ge.hpp"

using namespace arma;
using std::vector;
using std::cout;
using std::endl;

SEXP fast_lm(SEXP A_sexp, SEXP b_sexp) {
  SEXP ans;
  if(TYPEOF(A_sexp) != REALSXP || TYPEOF(b_sexp) != REALSXP) {
    cerr << "A and b must be vector<double>." << endl;
    return R_NilValue;
  }
  int NR = nrows(A_sexp);
  int NC = ncols(A_sexp);
  PROTECT(ans = allocVector(REALSXP,NC));
  mat A(REAL(A_sexp), NR, NC, false);
  vec b(REAL(b_sexp), NR, false);
  vec x;
  solve( x, A, b);
  for(int i = 0; i < NC; i++) {
    REAL(ans)[i] = x[i];
  }
  UNPROTECT(1);
  return ans;
}

double* getColFromName(SEXP x, const char* str) {
  double* ans = NULL;
  SEXP nms = getAttrib(x,R_NamesSymbol);
  for(R_len_t i = 0; i < length(nms); i++) {
    if(strcmp(CHAR(STRING_ELT(nms,i)),str) == 0) {
      ans = REAL(VECTOR_ELT(x,i));
      break;
    }
  }
  return ans;
}

SEXP fast_lm_dataframe(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp) {
  SEXP ans;
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  vector<double*> theData(NC);
  for(int i = 0; i < NC; i++) {
    theData[i] = getColFromName(panel_sexp,CHAR(STRING_ELT(left_hand_sides_sexp,i)));
  }
  PROTECT(ans = allocVector(REALSXP,NC));
  mat A(NR, NC);
  for(int i = 0; i < NC; i++) {
    std::copy(theData[i], theData[i] + NR, &(A.col(i)[0]));
  }
  vec b(getColFromName(panel_sexp,CHAR(STRING_ELT(right_hand_side_sexp,0))), NR, false);
  vec x;
  solve( x, A, b);
  for(int i = 0; i < NC; i++) {
    REAL(ans)[i] = x[i];
  }
  UNPROTECT(1);
  return ans;
}

void single_panel_lm(vec& x, vector<double*>& lhs_p, double* rhs, const int NR, double* scratch) {
  const int NC = lhs_p.size();
  vec b(rhs, NR, false);
  mat A(scratch, NR, NC, false);
  for(int i = 0; i < NC; i++) {
    std::copy(lhs_p[i], lhs_p[i] + NR, &(A.col(i)[0]));
  }
  solve( x, A, b);
}

// assumes that panel_sexp is already sorted by timestamp
SEXP expanding_lm_dataframe(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP min_rows_sexp) {
  SEXP ans;
  if(TYPEOF(min_rows_sexp) != INTSXP) {
    cerr << "min.rows must be an integer." << endl;
    return R_NilValue;
  }
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  int min_rows = INTEGER(min_rows_sexp)[0];
  if(min_rows > NR) {
    cerr << "nrow(panel) must be > min.rows." << endl;
    return R_NilValue;
  }
  double* rhs = getColFromName(panel_sexp,CHAR(STRING_ELT(right_hand_side_sexp,0)));
  vector<double*> theData(NC);
  for(int i = 0; i < NC; i++) {
    theData[i] = getColFromName(panel_sexp,CHAR(STRING_ELT(left_hand_sides_sexp,i)));
  }
  // this is the biggest regression we will do
  // during this call, smaller subsets can use the same scratch space
  double* scratch = new double[NR * NC];
  vec x;
  int ans_NR = NR - min_rows + 1;
  PROTECT(ans = allocMatrix(REALSXP, ans_NR, NC));
  double* ans_ptr = REAL(ans);
  for(int i = min_rows - 1; i < NR; i++) {
    single_panel_lm( x, theData, rhs, i + 1, scratch);
    for(int j = 0; j < NC; j++) {
      ans_ptr[j*ans_NR] = x[j];
    }
    ++ans_ptr;
  }
  delete []scratch;
  UNPROTECT(1);
  return ans;
}

// assumes that panel_sexp is already sorted by asofdate
SEXP expanding_panel_dataframe(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP asofdate_column_sexp, SEXP min_dates_sexp) {
  SEXP ans;
  if(TYPEOF(min_dates_sexp) != INTSXP) {
    cerr << "min.dates must be an integer." << endl;
    return R_NilValue;
  }
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  double* dts = getColFromName(panel_sexp,CHAR(STRING_ELT(asofdate_column_sexp,0)));
  vector<double> unique_dates;
  std::unique_copy(dts,dts+NR,back_inserter(unique_dates));
  int min_dates = INTEGER(min_dates_sexp)[0];
  if(min_dates > unique_dates.size()) {
    cerr << "size of unique dates must be > min.dates." << endl;
    return R_NilValue;
  }
  double* rhs = getColFromName(panel_sexp,CHAR(STRING_ELT(right_hand_side_sexp,0)));
  vector<double*> theData(NC);
  for(int i = 0; i < NC; i++) {
    theData[i] = getColFromName(panel_sexp,CHAR(STRING_ELT(left_hand_sides_sexp,i)));
  }
  // this is the biggest regression we will do
  // during this call, smaller subsets can use the same scratch space
  double* scratch = new double[NR * NC];
  vec x;
  int ans_NR = unique_dates.size()  - min_dates + 1;
  PROTECT(ans = allocMatrix(REALSXP, ans_NR, NC));
  double* ans_ptr = REAL(ans);

  // data cutoff for regression (start at row 0)
  double* subset = dts;
  for(size_t i = min_dates - 1; i < unique_dates.size(); i++) {
    subset = index_ge(subset, dts + NR, unique_dates[i]);
    single_panel_lm( x, theData, rhs, std::distance(dts,subset), scratch);
    for(int j = 0; j < NC; j++) {
      ans_ptr[j*ans_NR] = x[j];
    }
    ++ans_ptr;
  }
  delete []scratch;
  UNPROTECT(1);
  return ans;
}
