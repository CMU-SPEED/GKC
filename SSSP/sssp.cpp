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
#include <limits>
#include "utils.h"
#include "graph.h"
#include "sssp_checker.h"
#include <omp.h>

#ifndef ITERS
#define ITERS 64
#endif


#ifndef VALIDATE
#define VALIDATE 0
#endif

typedef uint32_t VTYPE;

extern uint32_t sssp( uint32_t * IA, uint32_t * JA, uint32_t *A,
  uint32_t N, uint32_t * lens, uint32_t src, uint32_t delta);

void usage(char * pname){
 fprintf(stderr, "USAGE: %s <IA fname> <JA fname> <delta> [<sources> <A fname>]\n", pname);
 exit(EXIT_FAILURE);
}


int main(int argc, char** argv){
 uint32_t * IA;
 uint32_t * JA;
 uint32_t * A;

 uint32_t * srcs;

 uint32_t N;
 uint32_t M;

 if (argc < 3) {
  usage(argv[0]);
 }

 // File read in:
 N = tell_size(argv[1])-1;
 M = tell_size(argv[2]);

 IA = (uint32_t *)malloc((N+1)*sizeof(uint32_t));
 JA = (uint32_t *)malloc(M*sizeof(uint32_t));
 A = (uint32_t *)malloc(M*sizeof(uint32_t));

 srcs = (uint32_t *)malloc(64*sizeof(uint32_t));

 if (!IA || !JA || !A ) {
  fprintf(stderr, "COULD NOT ALLOCATE MEMORY\n");
  exit(EXIT_FAILURE);
 }

 read_binary_buffers(argv[1], IA);
 read_binary_buffers(argv[2], JA);

 uint32_t delta = atoi(argv[3]);
 printf("DELTA = %u\n", delta); fflush(NULL);

 FILE *source_file;
 uint32_t i = 0;
 if (argc >= 5){
  printf("Reading Source Files\n");
  if ((source_file = fopen(argv[4], "r")) != NULL){
   uint32_t tmp_src;
   while ( fscanf(source_file, "%u\n", &tmp_src) != EOF){
    if ((tmp_src) >= N){
     printf("WARNING: sources file contained too-large source: %u\n", tmp_src);
     printf("Skipping it!\n");
    }
    else
    {
     srcs[i] = tmp_src;
     i += 1;
    }
   }
   fclose(source_file);
   printf("Sucessfully read %u sources.\n", i);
  }
 }

 if (argc < 6){
  for (uint32_t i = 0; i != M; ++i){
   A[i] = rand() % 254 + 1;
  }
 }
 else{
  read_binary_buffers(argv[5], A);
 }
 printf("Number of vertices: %u\n", N); fflush(NULL);

 uint32_t * lens;

 double st, nd;
 char * trunc_fname = truncate_fname(argv[1]);
 double tot_time = 0.0;

 uint32_t num_threads = omp_get_max_threads();

 printf("Start SSSP\n");
 printf("round, name, sourceID, time(s), threads\n");

 for (uint32_t iter=0; iter < ITERS; iter++){

  st = omp_get_wtime();
  lens = (uint32_t *)calloc(N, sizeof(uint32_t));
  uint32_t tmp = sssp(IA, JA, A, N, lens, srcs[iter], delta);
  nd = omp_get_wtime();

  printf("Round %u, %s, %u, %f sec, %u\n", iter, trunc_fname, srcs[iter], nd-st, num_threads);

#ifdef VALIDATE
  if ( check_dists(IA, JA, A, N, srcs[iter], lens) )
   printf("Passed\n");
  else
   printf("Failed\n");
#endif
  free(lens);
  tot_time += nd - st;

 }
 printf("Average time: %f seconds.\n\n", tot_time/ITERS);

 free(trunc_fname);
 free(IA);
 free(JA);
 free(A);
 free(srcs);
 return 0;
}
