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
 
#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <memory>
#include <cstdint>
#include <unistd.h>
#include "graph.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>
#include <vector>
#include "bfs_core.h"
#include "utils.h"
#include <omp.h>

#ifndef ITERS
#define ITERS 1
#endif
typedef uint32_t VTYPE;
typedef int64_t  PTYPE;

void usage(const char * exec_name)
{
 printf("USAGE: %s IA_FILE JA_FILE [optional:source_id(int)]\n", exec_name);
}

int main(int argc, char **argv){
 PTYPE * parent;
 double t0,t1;
 uint32_t NUM_EDGES,NUM_VERTICES;
 VTYPE * IAr, * JAr, * IAc, * JAc;
 FILE *source_file;
 std::vector<uint32_t> source_ids;
 uint32_t NUM_THREADS = omp_get_max_threads();

 /********************/
 /**** Argparsing ****/
 /********************/
 if (argc < 3)
 {
  usage(argv[0]);
  return 1;
 }

 // Read file sizes since the previous check means we at least have their names
 NUM_VERTICES=tell_size(argv[1])-1;
 NUM_EDGES=tell_size(argv[2]);
 if (argc >= 4){
  // Try to read file with source IDs. otherwise assume the arg is a single
  // source.
  if ((source_file = fopen(argv[3], "r")) != NULL){
   uint32_t tmp_src;
   while ( fscanf(source_file, "%u\n", &tmp_src) != EOF){
    if ((tmp_src) >= NUM_VERTICES){
     printf("WARNING: sources file contained too-large source: %u\n", tmp_src);
     printf("Skipping it!\n");
    }
    else
    {
     source_ids.push_back(tmp_src); 
    }
   }
   fclose(source_file);
   printf("Successfully read %ld sources.\n", source_ids.size());
  } else {
   uint32_t tmp_src = (uint32_t)atol(argv[3]);
   if (tmp_src >= NUM_VERTICES){
    printf("WARNING: you selected a source vertex id that does not exist. Defaulting to 0\n");
    tmp_src = 0;
   }
   source_ids.push_back(tmp_src);
  }
  printf("Ready to process %ld sources.\n", source_ids.size());
 }
 if (source_ids.size() == 0){
  source_ids.push_back(0);
  printf("WARNING: no sources provided. Processing source 0.\n");
 }

 /***************************************/
 /**** File Read in and Memory Alloc ****/
 /***************************************/
 printf("Read IA returns %u\n", NUM_VERTICES); 
 printf("Read JA returns %u\n", NUM_EDGES);
 // Compute rounded up number of entries for data
 // (such that each segment will be cache-aligned):
 uint64_t I_bytes = (NUM_VERTICES+1) * sizeof(uint32_t);
 I_bytes = ((I_bytes >> 6) << 6)  + 64 * ((I_bytes & 0x3F) != 0);
 uint64_t J_bytes = NUM_EDGES * sizeof(uint32_t);
 J_bytes  = ((J_bytes >> 6) << 6) + 64 * ((J_bytes & 0x3F) != 0);


 // Set pointers:
 IAr = (uint32_t*) aligned_alloc(64, I_bytes);
 JAr = (uint32_t*) aligned_alloc(64, J_bytes);
 if (!IAr || !JAr){
  fprintf(stderr, "ERROR: could not allocate memory for I and J arrays.\n");
  exit(EXIT_FAILURE);
 }

 read_binary_buffers(argv[1],IAr);
 read_binary_buffers(argv[2],JAr); 

 // Convert to CSC:
 IAc=NULL; JAc=NULL;
 if (!csr_to_csc_parallel(IAr, JAr, &IAc, &JAc, NUM_VERTICES)){
  fprintf(stderr, "ERROR: failed to transpose matrix!");
  exit(EXIT_FAILURE);
 }
 printf("Completed transpose. Moving to BFS.\n");

 std::cout << "Going to run with " << omp_get_max_threads() 
  << " threads." << std::endl;
 // ********** End of setup *********


 double total_time = 0;
 parent = (PTYPE * )malloc(NUM_VERTICES* sizeof(PTYPE));
 for (auto srcs_itr = source_ids.begin(); srcs_itr != source_ids.end(); srcs_itr++){
  uint32_t source_id = *srcs_itr;
  uint32_t depth;	
  double trial_time = 0;
  for (int i = 0; i < ITERS; i++){
   t0 = omp_get_wtime();
   depth = par_bfs(source_id,parent,IAr,JAr,IAc,JAc,NUM_VERTICES);
   t1 = omp_get_wtime();
   trial_time += (t1-t0);
  }
  double avg_time = (double)trial_time / (double)ITERS;
  total_time += avg_time;
  uint32_t num_unvisited = 0;
  for (uint32_t idx = 0; idx < NUM_VERTICES; idx++){
   if (parent[idx] == NUM_VERTICES) {
    num_unvisited++;
   }
  }
  printf("name,source,time_avg,unreached,depth,threads\n");
  printf("%s,%u,%f,%d,%u,%d\n", argv[1], source_id, avg_time, num_unvisited, depth, NUM_THREADS);

#ifdef CHECK_DEPTHS
  // Check depths:
  uint32_t * depth_table = (uint32_t * )malloc(NUM_VERTICES* sizeof(uint32_t ));
  init_vector(depth_table, NUM_VERTICES, NUM_VERTICES);
  make_depth_table( source_id, depth_table, IAr, JAr, NUM_VERTICES);
  if (!check_parents_vs_depths(parent, depth_table, source_id, NUM_VERTICES)){
   std::cerr << "FAILED PARENT VS DEPTH CHECK" << std::endl;
  }
  else {
   std::cerr << "PASSED PARENT VS DEPTH CHECK." << std::endl;
  }
  free(depth_table);
#endif
 }
 printf("Average time for all sources: %f\n", total_time/source_ids.size());
 
 free(parent);
 free(IAr);
 free(JAr);
 free(IAc);
 free(JAc);
}
