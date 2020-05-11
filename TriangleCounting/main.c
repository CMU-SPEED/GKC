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
 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "utils.h"
#include "graph.h"
#include "tc.h"
#include <omp.h>
#include <math.h>
#include "tri_count_checker.h"

#ifndef ITERS
#define ITERS 3
#endif


void usage(char * pname){
 fprintf(stderr, "USAGE: %s <IA fname> <JA fname>\n", pname);
 exit(EXIT_FAILURE);
}

int main(int argc, char** argv){
 uint32_t * IA;
 uint32_t * JA;
 uint32_t N;
 uint32_t M;
 double st, nd;

 if (argc != 3) {
  usage(argv[0]);
 }

 // File read in:
 N = tell_size(argv[1])-1;
 M = tell_size(argv[2]); 

 IA = (uint32_t *)malloc((N+1)*sizeof(uint32_t));
 JA = (uint32_t *)malloc(M*sizeof(uint32_t));

 if (!IA || !JA ) {
  fprintf(stderr, "COULD NOT ALLOCATE MEMORY\n");
  exit(EXIT_FAILURE);
 }

 read_binary_buffers(argv[1], IA);
 read_binary_buffers(argv[2], JA);

 char * trunc_fname = truncate_fname(argv[1]);


 // *************** Begin Processing ****************************
 uint64_t delta;
 double avg_time = 0;
 for (uint32_t p_iter = 0; p_iter < ITERS; p_iter++){
  st = omp_get_wtime();
  delta = tri_count(IA, JA, N);
  nd = omp_get_wtime();
  printf("TRIAL: %lu triangles in %f seconds\n", delta, nd-st);
  avg_time += nd-st;
#ifdef VERIFY
  if (check_tri_count(delta, trunc_fname)){
   printf("PASSED CHECK\n");
  }
  else {
   printf("FAILED CHECK\n");
  }
#endif
 }
 avg_time /= ITERS;
 printf("Average time: %f seconds\n", avg_time);
 free(IA);
 free(JA);
 free(trunc_fname);

 return 0;
}

