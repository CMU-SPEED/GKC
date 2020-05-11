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
 
/*This file contains utitlity functions to operate on
adjaceny matrices and frontier vectors*/
#ifndef GRAPH_HEADER
#define GRAPH_HEADER
#include <stdlib.h>
#include <stdbool.h>
#include <memory.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include "utils.h"

// returns number of lines in <filename>
uint32_t tell_size(const char * filename);

// need not really return anything
// reads 0-indexed file
uint32_t read_binary( const char * filename, uint32_t *array);

// Read buffered. Faster!
uint32_t read_binary_buffers( const char * filename, uint32_t *array);

 //create adjacency matrix
uint32_t init( char ** argv, uint32_t ** IA, uint32_t ** JA);

void print_v( uint32_t * v, uint32_t length );

void print_v_f(float * v, uint32_t length );

void init_vector(uint32_t * v, uint32_t length,  uint32_t val);

void init_vector_f(float * v, uint32_t length,  float val);

uint32_t nz_v(uint32_t *v, uint32_t val, uint32_t length);

// convert to CSC
bool csr_to_csc(uint32_t *IAr, uint32_t * JAr, uint32_t ** IAc, uint32_t ** JAc, uint32_t length);
bool csr_to_csc_parallel(uint32_t *IAr, uint32_t * JAr, uint32_t ** IAc, uint32_t ** JAc, uint32_t length);

// Full symmetric matrix to lower tri:
void csr_to_lower(uint32_t * IAf, uint32_t * JAf, 
  uint32_t * IAl, uint32_t * JAl, uint32_t N);

void sort_neighborhoods(uint32_t * IA, uint32_t * JA, uint32_t N);

// Convert IA and JA to center-offset IA and JA, such that the values in 
// IA and JA can be signed integers and the centered pointers point into 
// the middle of each array.
// NOTE: IA_cent and JA_cent are allocated within this function! Don't pre-allocate!
void csr_to_center_csr(uint32_t * IA, uint32_t * JA, 
  int32_t ** IA_cent, int32_t ** JA_cent, uint32_t N);

#endif
