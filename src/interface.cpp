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
SEXP expanding_panel(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP asofdate_column_sexp, SEXP min_dates) {
  SEXP ans;
  if(TYPEOF(min_dates) != INTSXP) {
    cerr << "min_dates must be an integer." << endl;
  }
  R_len_t NR = length(VECTOR_ELT(panel_sexp,0));
  R_len_t NC = length(left_hand_sides_sexp);
  SEXP nms = getAttrib(panel_sexp, R_NamesSymbol);
  double* dts = getColFromName(panel_sexp,"asofdate");
  vector<double*> theData;
  for(int i = 0; i < NC; i++) {
    theData.push_back(getColFromName(panel_sexp,CHAR(STRING_ELT(left_hand_sides_sexp,i))));
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
