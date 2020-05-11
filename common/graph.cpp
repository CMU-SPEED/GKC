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
 
/*Implementation of Graph utitlty functions as
described in "graph.h"*/
#include "graph.h"


uint32_t tell_size(const char * filename) {
  FILE *fptr= fopen(filename,"rb");
  uint32_t n;
  if ( fread(&n, sizeof(n), 1, fptr) != 1) {
    fprintf(stderr, "ERROR reading size from bin file!\n");
    exit(EXIT_FAILURE);
  }
  fclose(fptr);
  return n;
}
uint32_t read_binary( const char * filename, uint32_t *array) {
  FILE *fptr= fopen(filename,"rb");
  fseek(fptr,0,SEEK_END);
  uint64_t end = ftell(fptr);
  fseek(fptr,0,SEEK_SET);
  printf("File has %lu entries.\n", end/4);
  uint64_t counter =0;
  uint32_t n;

  if ( fread(&n, sizeof(n), 1, fptr) != 1) {
    fprintf(stderr, "ERROR reading size from bin file!\n");
    exit(EXIT_FAILURE);
  }
  printf("First 4 bytes indicate %u entries.\n", n);
  counter+=sizeof(uint32_t);
  fseek(fptr,counter,SEEK_SET);
  uint32_t el;
  while(counter<end) {
    if ( fread(&el,sizeof(el),1,fptr) != 1){
      fprintf(stderr, "ERROR reading data from bin file!\n");
      exit(EXIT_FAILURE);
    }
    counter+=sizeof(el);
    fseek(fptr,counter,SEEK_SET);
    array[counter/4-2]=el;
  }

  fclose(fptr);
  return n;
}

uint32_t read_binary_buffers( const char * filename, uint32_t *array) {
  FILE *fptr= fopen(filename,"rb");
  uint32_t n;
  // reads first element from file (number of elements to follow)
  size_t s = fread(&n, sizeof(uint32_t), 1, fptr);
  printf(" %u elements to read \n", n);

  uint32_t stride = 500000000;
  uint32_t index = 0;
  uint32_t loop_iters = n/stride + (n % stride > 0);
  printf(" %u reads\n ", loop_iters);
  for(uint32_t i = 0; i < loop_iters ; i++){
    uint32_t read_amt = MIN(stride, n-index);
    s = fread(array+index, sizeof(uint32_t), read_amt, fptr);
    index += read_amt;
    printf("%lu bytes read\n", (uint64_t)s*sizeof(uint32_t));
  }
  // after reading the number of elements, the file stream points to the
  // second element in the file.
  // read the rest of the elements into the array.


  fclose(fptr);
  return n;
}

uint32_t init( char ** argv, uint32_t ** IA, uint32_t ** JA){

  uint32_t NUM_VERTICES=tell_size(argv[1])-1;
  uint32_t NUM_EDGES=tell_size(argv[2]);

  // Allocate contiguous sections of memory for each thread in the test:

  // for non-triangularized data:
  char *I_J;

  // Compute rounded up number of entries for data
  // (such that each segment will be cache-aligned):
  uint64_t I_bytes = (NUM_VERTICES+1) * sizeof(uint32_t);
  I_bytes = ((I_bytes >> 6) << 6)  + 64 * ((I_bytes & 0x3F) != 0);
  uint64_t J_bytes = NUM_EDGES * sizeof(uint32_t);
  J_bytes  = ((J_bytes >> 6) << 6) + 64 * ((J_bytes & 0x3F) != 0);
  // Num bytes for original data:
  uint64_t bytes_aligned = I_bytes + J_bytes;

  // Allocate the contiguous segments of memory:
  I_J = (char*)aligned_alloc(64, bytes_aligned);
  // Set pointers:
  *IA = (uint32_t*) (I_J + 0);
  *JA = (uint32_t*) (I_J + I_bytes);

  read_binary_buffers(argv[1],*IA);
  read_binary_buffers(argv[2],*JA);

  return NUM_VERTICES;
}


void print_v( uint32_t * v, uint32_t length ){
  for(uint32_t i = 0; i< length /*&& v[i]<length*/ ; i++){
    printf(" %d ", v[i]);
  }
  printf("\n" );
}

void print_v_f(float * v, uint32_t length ){
  for(uint32_t i = 0; i< length /*&& v[i]<length*/ ; i++){
    printf(" %f \n", v[i]);
  }
  printf("\n" );
}

void init_vector(uint32_t * v, uint32_t length,  uint32_t val){
  for(uint32_t i = 0; i< length; i++){
    v[i]= val;
  }
}
void init_vector_f(float * v, uint32_t length,  float val){
  for(uint32_t i = 0; i< length; i++){
    v[i]= val;
  }
}

uint32_t nz_v(uint32_t *v, uint32_t val, uint32_t length){
  uint32_t nz = 0;
  for(uint32_t i = 0; i < length; i++ ){
    if(v[i] == val)
      nz++;
  }
  return nz;
}

bool csr_to_csc(uint32_t *IAr, uint32_t * JAr, uint32_t ** IAc, uint32_t ** JAc, uint32_t length){
  // Memory alloc:
  if ((*IAc)==NULL || (*JAc)==NULL){
		printf("one or both null\n");
    if (*IAc!=NULL || *JAc!=NULL) return false; // Either both or none should be allocated.
    // Allocate cache-aligned data
    printf("Allocating memory for transpose.\n");
		uint32_t edges = IAr[length];
    unsigned long long I_bytes = (length+1) * sizeof(uint32_t);
    I_bytes = ((I_bytes >> 6) << 6)  + 64 * ((I_bytes & 0x3F) != 0);
    unsigned long long J_bytes = edges * sizeof(uint32_t);
    J_bytes  = ((J_bytes >> 6) << 6) + 64 * ((J_bytes & 0x3F) != 0);

    *IAc = (uint32_t*)aligned_alloc(64, I_bytes);
    *JAc = (uint32_t*)aligned_alloc(64, J_bytes);
    if (!*IAc || !*JAc) return false;
  }
  printf("Transposing.\n");
  // Transpose:
  std::vector<std::vector<uint32_t>> remap(length);
  for (uint32_t idx = 0; idx < length; idx++){
    uint32_t st = IAr[idx];
    uint32_t nd = IAr[idx+1];
    for (uint32_t edx = st; edx < nd; edx++) {
      uint32_t jdx = JAr[edx];
      remap[jdx].push_back(idx);
    }  
  }
  // Sort and copy (note swapped jdx and idx):
  *(IAc)[0] = 0;
  for (uint32_t jdx = 0; jdx < length; jdx++){
    auto st = remap[jdx].begin();
    auto nd = remap[jdx].end();
    std::sort(st, nd);
    // Refresh start and end iterators (is that even necessary?)
    st = remap[jdx].begin();
    nd = remap[jdx].end();
    // Cumulative sum:
    (*IAc)[jdx+1] = (*IAc)[jdx] + remap[jdx].size();
    // Copy to CSC:
    for (uint32_t edx = 0; edx < remap[jdx].size(); edx++) {
      (*JAc)[(*IAc)[jdx] + edx] = remap[jdx][edx];
    }
  }
  return true;
}

bool csr_to_csc_parallel(uint32_t *IAr, uint32_t * JAr, uint32_t ** IAc, uint32_t ** JAc, uint32_t length){
 // Memory alloc:
  uint32_t edges = IAr[length];
  unsigned long long I_bytes = (length+1) * sizeof(uint32_t);
  I_bytes = ((I_bytes >> 6) << 6)  + 64 * ((I_bytes & 0x3F) != 0);
  unsigned long long J_bytes = edges * sizeof(uint32_t);
  J_bytes  = ((J_bytes >> 6) << 6) + 64 * ((J_bytes & 0x3F) != 0);
 if ((*IAc)==NULL || (*JAc)==NULL){
  printf("one or both null\n");
  if ((*IAc) || (*JAc)) return false; // Either both or none should be allocated.
  // Allocate cache-aligned data
  printf("Allocating memory for transpose.\n");

  *IAc = (uint32_t*)aligned_alloc(64, I_bytes);
  *JAc = (uint32_t*)aligned_alloc(64, J_bytes);
  if (!(*IAc) || !(*JAc)) return false;
 }
 uint32_t * degrees = (uint32_t *)calloc(length, sizeof(uint32_t));
 uint32_t * IAc_main = *IAc;
 if (!degrees) return false;
 
 printf("Transposing.\n");
 // First get size of incoming neighborhoods:
#pragma omp parallel for schedule(dynamic, 64)
 for (uint32_t edx = 0; edx < edges; edx++) {
  uint32_t jdx = JAr[edx];
#pragma omp atomic
  degrees[jdx]++;
 }

 uint32_t num_elems = 1024*1024 / sizeof(uint32_t); // About the size of L2 cache on SKX
	uint32_t num_blocks= (length+num_elems-1) / num_elems;

	// Coarse parallel sum
	IAc_main[0] = 0;
#pragma omp parallel for
	for (uint32_t bidx = 0; bidx < num_blocks; bidx++){
		uint32_t b_st = bidx * num_elems;
		uint32_t b_nd = MIN(length, (bidx+1) * num_elems);
		uint32_t t_sum = 0;
		for (uint32_t idx = b_st; idx < b_nd; idx++){
			// TODO: SIMD reduction
			t_sum += degrees[idx];
			IAc_main[idx+1] = t_sum;
		}
	}
	// At this point each b_nd location of IAc_main
	// has a prefix sum for a block
	
	// Coarse prefix summation of sums at each b_nd:
	for (uint32_t bidx = 1; bidx < num_blocks; bidx++){
		uint32_t b_st = bidx * num_elems;
		uint32_t b_nd = MIN(length, (bidx+1) * num_elems);
		IAc_main[b_nd] += IAc_main[b_st];
	}

	// Now each block has a local sum, except that the entry at the 
	// block start (1 + last block end) has a global prefix sum value

	// Parallel local summation
#pragma omp parallel for
	for (uint32_t bidx = 1; bidx < num_blocks; bidx++){
		uint32_t b_st = bidx * num_elems;
		uint32_t b_nd = MIN(length, (bidx+1) * num_elems);
		uint32_t t_sum = IAc_main[b_st];
		// TODO: SIMD
		for (uint32_t idx = b_st+1; idx < b_nd; idx++){
			IAc_main[idx] += t_sum;
		}
	}

 // IAc_shadow can be used as a temporary for size of unwritten neighborhoods:
 uint32_t * IAc_shadow = (uint32_t*)aligned_alloc(64, I_bytes);
#pragma omp parallel for
 for (uint32_t idx = 0; idx < length; idx++){
  IAc_shadow[idx] = 0;
 }

 // Then copy source over (expensive :( )
#pragma omp parallel for schedule(dynamic, 64)
 for (uint32_t idx = 0; idx < length; idx++){
  uint32_t st = IAr[idx];
  uint32_t nd = IAr[idx+1];
  for (uint32_t edx = st; edx < nd; edx++) {
   uint32_t jdx = JAr[edx];
   uint32_t nbhd_st = IAc_main[jdx];
   uint32_t curr_nbhd_sz = __sync_fetch_and_add(&(IAc_shadow[jdx]), 1);
   (*JAc)[nbhd_st + curr_nbhd_sz] = idx;
  }  
 }
 // Sort new neighborhoods:
#pragma omp parallel for schedule(dynamic, 64)
 for (uint32_t idx = 0; idx < length; idx++){
  uint32_t j_st = IAc_main[idx];
  uint32_t len = IAc_main[idx+1] - IAc_main[idx];
  qsort(*JAc + j_st, len, sizeof(uint32_t), cmp_u32);
 }

 free(degrees);
 free(IAc_shadow);
 return true;
}

// Convert symmetric, full matrix from CSR/CSC to lower triangular CSR.
void csr_to_lower(uint32_t * IAf, uint32_t * JAf, 
  uint32_t * IAl, uint32_t * JAl, uint32_t N){
  uint32_t Ml = 0;
  IAl[0] = 0;
  for (uint32_t idx = 0; idx < N; idx++){
    uint32_t st = IAf[idx];
    uint32_t nd = IAf[idx+1];
    //uint32_t old_M = Ml;
    for (auto j = st; j != nd; j++){
      if ( JAf[j] < idx) {
        JAl[Ml] = JAf[j];
        Ml++; 
      }
      else {break;}
    }
    //if ((long)Ml-(long)old_M != 0){
    //printf("%u: Filled half neighborhood with %u of %u edges.\n",
    //  idx, Ml-old_M, nd-st);
    //}
    IAl[idx+1] = Ml;
  }
  printf("Got lower triangular (%u of %u edges)\n", IAl[N], IAf[N]);
}


// Inputs: IA array, JA array, and N=number of vertices. (IA is size N+1)
void sort_neighborhoods(uint32_t * IA, uint32_t * JA, uint32_t N){
  for (uint32_t idx = 0; idx < N; idx++){
    uint32_t st = IA[idx];
    uint32_t nd = IA[idx+1]; 
    qsort(JA+st, nd-st, sizeof(uint32_t), cmp_u32);
  }
}

// See detailed description of this function in header.
void csr_to_center_csr(uint32_t * IA, uint32_t * JA, 
  int32_t ** IA_cent, int32_t ** JA_cent, uint32_t N){
  
  uint32_t M = IA[N]; // Retrieve num edges
  
  // IA_cent and JA_cent should not already be allocated.
  assert(!(*IA_cent) && !(*JA_cent));

  // Now allocate
  int32_t * IA_int = (int32_t *)malloc(sizeof(int32_t)*(N+1));
  int32_t * JA_int = (int32_t *)malloc(sizeof(int32_t)*M);
  assert(IA_int!=NULL && JA_int!=NULL);

  // Set pointers to centered offsets
  *IA_cent = IA_int + (N+1) / 2;
  *JA_cent = JA_int + M / 2;

  // Copy IA and JA into offset arrays, modifying values as needed.
  for (uint32_t idx = 0; idx < N+1; idx++){
    IA_int[idx] = (int32_t)((int64_t)IA[idx] - (int64_t)(M/2));
  }

  for (uint32_t idx = 0; idx < M; idx++){
    JA_int[idx] = (int32_t)((int64_t)JA[idx] - (int64_t)((N+1)/2));
  }
}
