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
#include "bc.h"
#include "bc_checking.h"

#define NUM_SRCS 64
#define SRCS_PER_TRIAL 4

void usage(char * pname){
 fprintf(stderr, "USAGE: %s <IA fname> <JA fname> [optional: source or sources file] \n", pname);
 exit(EXIT_FAILURE);
}

int main(int argc, char** argv){
 uint32_t N;
 uint32_t M;
 uint32_t * IA;
 uint32_t * JA;
 CENT_T * centralities = NULL;
 FILE *source_file;
 std::vector<uint32_t> source_ids;

 /********************/
 /**** Argparsing ****/
 /********************/
 if (argc < 3)
 {
  usage(argv[0]);
  return 1;
 }

 // Read file sizes since the previous check means we at least have their names
 N=tell_size(argv[1])-1;
 M = tell_size(argv[2]);	
 if (argc >= 4){
  // Try to read file with source IDs. otherwise assume the arg is a single
  // source.
  if ((source_file = fopen(argv[3], "r")) != NULL){
   uint32_t tmp_src;
   while ( fscanf(source_file, "%u\n", &tmp_src) != EOF){
    if ((tmp_src) >= N){
     printf("WARNING: sources file contained too-large source: %u\n", tmp_src);
     printf("Skipping it!\n");
    }
    else
    {
     source_ids.push_back(tmp_src); 
    }
   }
   fclose(source_file);
   printf("Sucessfully read %lu sources.\n", source_ids.size());
  } else {
   uint32_t tmp_src = (uint32_t)atol(argv[3]);
   if (tmp_src >= N){
    printf("WARNING: you selected a source vertex id that does not exist. Defaulting to 0\n");
    tmp_src = 0;
   }
   source_ids.push_back(tmp_src);
  }
  printf("Ready to process %lu sources.\n", source_ids.size());
 }
 if (source_ids.size() == 0){
  source_ids.push_back(0);
  printf("WARNING: no sources provided. Processing source 0.\n");
 }

 // File read in:
 IA = (uint32_t *)malloc((N+1)*sizeof(uint32_t));
 JA = (uint32_t *)malloc(  M  *sizeof(uint32_t));
 uint32_t * sources = source_ids.data();
 uint32_t num_srcs = MIN(NUM_SRCS, source_ids.size());

 if (!IA || !JA ) {
  fprintf(stderr, "COULD NOT ALLOCATE MEMORY\n");
  exit(EXIT_FAILURE);
 }

 read_binary_buffers(argv[1], IA);
 read_binary_buffers(argv[2], JA);

 double st, nd;
 char * trunc_fname = truncate_fname(argv[1]);
 double avg_time = 0;	

 PATH_TYPE * paths_out = NULL;
 printf("name,time(s),threads\n");
 uint32_t num_threads = omp_get_max_threads();
 printf("input,source1,source2,source3,source4,time,threads\n");
 for (uint32_t src_idx = 0; src_idx < NUM_SRCS; src_idx+=SRCS_PER_TRIAL){
  st = omp_get_wtime();
  brandes_centralities(IA,JA,N,centralities,sources+src_idx, SRCS_PER_TRIAL);
  nd = omp_get_wtime();
  avg_time += nd-st;
  printf("%s,%u,%u,%u,%u,%f,%u\n", 
    trunc_fname, 
    sources[src_idx], sources[src_idx+1], sources[src_idx+2], sources[src_idx+3],
    nd-st, num_threads
    );
#ifdef VERIFY
  if (BC_checker(IA, JA, paths_out, centralities, sources+src_idx, SRCS_PER_TRIAL, N)){
   printf("++++++SUCCESS: PASSED CHECK!+++++\n");
  } else {
   printf("~~~~~~~~~FAILED CHECK!!~~~~~~~~~~\n");
  }
#endif
	free(centralities);
 }


 printf("Average time: %f seconds.\n\n", avg_time/(NUM_SRCS/SRCS_PER_TRIAL));
 free(trunc_fname);


 free(IA);
 free(JA);

 return 0;
}

