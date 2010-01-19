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

// assumes that panel_sexp is already sorted by asofdate
SEXP panel_lm(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP asofdate_column_sexp) {
  SEXP ans;
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  SEXP nms = getAttrib(panel_sexp, R_NamesSymbol);
  double* dts = getColFromName(panel_sexp,"asofdate");
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

void single_panel_lm(vec& x, vector<double*>& lhs_p, double* rhs, const int NR) {
  const int NC = lhs_p.size();
  vec b(rhs, NR, false);
  mat A(NR, NC);
  for(int i = 0; i < NC; i++) {
    std::copy(lhs_p[i], lhs_p[i] + NR, &(A.col(i)[0]));
  }
  solve( x, A, b);
}

// assumes that panel_sexp is already sorted by asofdate
SEXP expanding_panel(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP asofdate_column_sexp, SEXP min_dates_sexp) {
  SEXP ans;
  if(TYPEOF(min_dates_sexp) != INTSXP) {
    cerr << "min.dates_sexp must be an integer." << endl;
    return R_NilValue;
  }
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  int min_dates = INTEGER(min_dates_sexp)[0];
  cout << "NR:" << NR << endl;
  cout << "min_dates:" << min_dates << endl;
  if(min_dates > NR) {
    cerr << "rnow(panel) must be > min.dates." << endl;
    return R_NilValue;
  }

  SEXP nms = getAttrib(panel_sexp, R_NamesSymbol);
  double* dts = getColFromName(panel_sexp,"asofdate");
  double* rhs = getColFromName(panel_sexp,CHAR(STRING_ELT(right_hand_side_sexp,0)));
  vector<double*> theData(NC);
  for(int i = 0; i < NC; i++) {
    theData[i] = getColFromName(panel_sexp,CHAR(STRING_ELT(left_hand_sides_sexp,i)));
  }
  vec x;
  PROTECT(ans = allocMatrix(REALSXP, NR - min_dates + 1, NC));
  double* ans_ptr = REAL(ans);
  for(int i = min_dates; i < (NR-1); i++) {
    cout << "i:" << i << endl;
    single_panel_lm( x, theData, rhs, i);
    for(int j = 0; j < NC; j++) {
      cout << "j:" << j << endl;
      ans_ptr[j*NR] = x[j];
    }
    ++ans_ptr;
  }
  UNPROTECT(1);
  return ans;
}
