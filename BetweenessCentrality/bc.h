/*
 * Graph Kernel Collection
 *
 * Copyright 2020 Carnegie Mellon University.
 *
 * NO WARRANTY. THIS CARNEGIE MELLON UNIVERSITY AND SOFTWARE ENGINEERING
 * INSTITUTE MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. CARNEGIE MELLON 
 * UNIVERSITY MAKES NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, 
 * AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR 
 * PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF 
 * THE MATERIAL. CARNEGIE MELLON UNIVERSITY DOES NOT MAKE ANY WARRANTY OF ANY
 * KIND WITH RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT 
 * INFRINGEMENT.
 *
 * Released under a BSD (SEI)-style license, please see license.txt or
 * contact permission@sei.cmu.edu for full terms.
 *
 * [DISTRIBUTION STATEMENT A] This material has been approved for public
 * release and unlimited distribution.  Please see Copyright notice for 
 * non-US Government use and distribution.
 *
 * This Software includes and/or makes use of the following Third-Party
 * Software subject to its own license:
 *
 * 1. Matrix Market Loader code (https://math.nist.gov/MatrixMarket/mmio-c.html).
 *
 *      The code made publicly available at nist.gov is not marked with a 
 *      copyright notice and is therefore believed pursuant to section 105 of 
 *      the Copyright Act, to not be entitled to domestic copyright protection 
 *      under U.S. law and is therefore in the public domain.  Accordingly, it 
 *      is believed that no license is required for its use.
 *
 * This Software may include certain portions of copyrighted code that is 
 * initially being released only in binary form for validation and evaluation
 * purposes. It is expected that source code will be released as open source at
 * a future date. 
 *
 * DM20-0375
 */
 
#ifndef __BETWEENNESS_CENTRALITY_H__
#define __BETWEENNESS_CENTRALITY_H__
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <unistd.h>
#include <memory>
#include <fstream>
#include <queue>
#include <algorithm>
#include <immintrin.h>
#include <omp.h>
#include "bc.h"
#include "graph.h"
#include "utils.h"

typedef uint32_t VTYPE;
typedef double PATH_TYPE;
typedef double CENT_T;
typedef CENT_T DELTA_T;

inline void min_CAS(uint32_t * loc, uint32_t cmp, uint32_t swp){
 while (swp < cmp){
  uint32_t tmp = __sync_val_compare_and_swap(loc, cmp, swp);
  if (tmp == cmp) {
   break;
  }
  cmp = tmp;
 }
}

inline void max_CAS(uint32_t * loc, uint32_t cmp, uint32_t swp){
 while (swp > cmp){
  uint32_t tmp = __sync_val_compare_and_swap(loc, cmp, swp);
  if (tmp == cmp) {
   break;
  }
  cmp = tmp;
 }
}

/*
 * Main Brandes Wrapper
 */
void brandes_centralities(
  uint32_t * IA, uint32_t * JA, 
  uint32_t N, CENT_T *& centralities,
  uint32_t * sources, uint32_t num_srcs);

/* 
 * Parallel BFS (used in Brandes)
 */
VTYPE par_bc(
		VTYPE source_id, 
		VTYPE * depths, 
		PATH_TYPE * paths,
		VTYPE * frontiers_IA, 
		VTYPE * frontiers_JA, 
		uint8_t * JA_copy, 
		VTYPE * IA, 
		VTYPE * JA, 
		VTYPE NUM_VERTICES);
#endif
