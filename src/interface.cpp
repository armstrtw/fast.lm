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

#include <armadillo>
#include <Rinternals.h>
#include "interface.hpp"

using namespace arma;
using std::cout;
using std::endl;

SEXP fast_lm(SEXP A_sexp, SEXP b_sexp) {
  SEXP ans;
  if(TYPEOF(A_sexp) != REALSXP || TYPEOF(b_sexp) != REALSXP) {
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
