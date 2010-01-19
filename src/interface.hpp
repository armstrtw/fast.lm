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

#ifndef INTERFACE_HPP
#define INTERFACE_HPP

extern "C" {
  SEXP fast_lm(SEXP A_sexp, SEXP b_sexp);
  SEXP expanding_panel(SEXP panel_sexp, SEXP right_hand_side_sexp, SEXP left_hand_sides_sexp, SEXP asofdate_column_sexp, SEXP min_dates);
}

#endif // INTERFACE_HPP
