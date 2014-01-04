/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* 
 * PAlgebra.cpp - Implementation of the class PAlgebra
 *
 * The class PAlgebra is the base class containing the structure of (Z/mZ)^*,
 * which is isomorphic to the Galois group over A = Z[X]/Phi_m(X)).
 */
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_p.h>
#include <NTL/GF2XFactoring.h>

NTL_CLIENT

#include <climits>     // defines INT_MAX
#include <cstring>
#include <algorithm>   // defines count(...), min(...)
#include <iostream>

#include "NumbTh.h"    // defines argmax(...)
#include "PAlgebra.h"


// Generate the representation of Z_m^* for a given odd integer m
void PAlgebra::init(unsigned mm, unsigned g) {
  if (m==mm) return; // nothing to do

  if (mm>NTL_SP_BOUND) return; //(mm&1)==0 ||

  m = mm;
  this->g = g;
  zz_p::init(mm);

  long idx;

  zmsIdx.assign(m,-1);  // allocate m slots, initialize them to -1
  for (unsigned i=idx=0; i<m; i++) if (GCD(i,m)==1) zmsIdx[i] = idx++;
  this->phim = idx;

  Phi_mX = Cyclotomic(m); // compute and store Phi_m(X)
}
