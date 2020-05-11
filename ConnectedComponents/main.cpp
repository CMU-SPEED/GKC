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
 
#include "graph.h"
#include "utils.h"
#include <omp.h>
#include <math.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef ITERS
#define ITERS 16
#endif

bool check_CC(uint32_t * IAr, uint32_t * JAr,
 	      uint32_t * IAc, uint32_t * JAc, 
	      uint32_t N, uint32_t * labels_to_check);

uint32_t CC(uint32_t * IA, uint32_t * JA, 
		       uint32_t * IAc, uint32_t * JAc,
		       uint32_t N, uint32_t * parents);

void usage(char * pname){
	fprintf(stderr, "USAGE: %s <IA fname> <JA fname>\n", pname);
	exit(EXIT_FAILURE);
}

int main(int argc, char ** argv){
  uint32_t *IAc;
  uint32_t *JAc;
  uint32_t *IA;
  uint32_t *JA;
  uint32_t iter=1, last_size;

  uint32_t NUM_EDGES;

  if (argc < 3)
    {
      usage(argv[0]);
      return 1;
    }


  uint32_t N = tell_size(argv[1])-1;
  uint32_t M = tell_size(argv[2]);

  IA = (uint32_t *)malloc((N+1)*sizeof(uint32_t));
  JA = (uint32_t *)malloc(M*sizeof(uint32_t));


  if (!IA || !JA ) {
    fprintf(stderr, "COULD NOT ALLOCATE MEMORY\n");
    exit(EXIT_FAILURE);
  }

  read_binary_buffers(argv[1], IA);
  read_binary_buffers(argv[2], JA);

  IAc = IA;
  JAc = JA;
  printf(" %s %u nodes %u edges\n", argv[1], N, IAc[N]);

  uint32_t * parents;

  double st, nd;
  char * trunc_fname = truncate_fname(argv[1]);
  double tot_time = 0.0;

  uint32_t num_threads = omp_get_max_threads();

  printf("Start Connected Component\n");
  printf("round, name, num of cc, time(s), threads\n");

  for (uint32_t iter=0; iter < ITERS; iter++){

    st = omp_get_wtime();
    parents = (uint32_t *)malloc(N * sizeof(uint32_t));
    uint32_t num_comps = CC(IA,JA,IAc,JAc,N,parents);
    nd = omp_get_wtime();

    printf("Round %u, %s, %u, %f sec, %u\n", iter, trunc_fname, num_comps, nd-st, num_threads);

#ifdef VERIFY	
    if (iter == ITERS - 1)
      {
	if (check_CC(IA, JA, IAc, JAc, N, parents))
	  printf("Passed\n");
	else
	  printf("Failed\n");
	
      }
#endif

    free(parents);
    tot_time += nd - st;

  }
  printf("Average time: %lf seconds.\n\n", tot_time/ITERS);

  free(trunc_fname);
  free(IA);
  free(JA);  

  return 0;
}
