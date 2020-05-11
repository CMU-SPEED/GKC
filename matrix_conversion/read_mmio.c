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
 
/*   15/01/2020 Mark Blanco, SPEED Group, CMU
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and write it out to binary-encoded data files for I, J coords (as 32-bit
 *   uints) and another file for values (as WTYPEs).
 *
 *   Each file begins with a 32-bit entry giving the number of remaining
 *   entries in that file. 
 *
 *   Usage:  a.out [filename] <optional output prefix>
 *
 *       
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
 */

#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include <algorithm>
#include <vector>
#include <set>
#include <climits>
#include <cassert>
#include <typeinfo>

typedef uint32_t VTYPE;
typedef uint32_t WTYPE;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
typedef std::pair<VTYPE, WTYPE> j_val;
bool cmp(const j_val & a, const j_val & b){
  return a.first < b.first;
}

// Sort elements in vector, then remove duplicates and resize the vector to
// match the new number of elements
// Return value is the number of new elements remaining.
uint32_t sort_and_dedup_vec(std::vector<j_val> & vec){
  if (vec.size() == 0) return 0;
  std::sort(vec.begin(), vec.end(), cmp);
  uint32_t new_size = 1; // First element is by definition not a duplicate.
  j_val old_val=vec[0];
  for (uint32_t idx = 1; idx < vec.size(); idx++){
    // Note: we don't care about uniqueness or lack thereof
    // for weights. Just that there is only one of each edge.
    // Hence just compare first entry in each pair.
    if ( vec[idx].first != old_val.first){
      old_val = vec[idx];
      vec[new_size] = vec[idx];
      new_size++;
    } 
  }
  // Update size of vector
  vec.resize(new_size);
  return new_size;
}


int main(int argc, char *argv[])
{
  int ret_code;
  bool symm = false;
  MM_typecode matcode;
  FILE *f;
  FILE *ia, *ja;
  uint64_t M, N, nz;   
  uint64_t i;
  WTYPE *val;
  char basename[256];
  bool symmcheck;


  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s [martix-market-filename] (optional:base-output-name) (optional:symmetry-flag[0|1])\n", argv[0]);
    fprintf(stderr, "NOTE: this program can symmetrize a non-symmetric matrix, but will write out the full matrix. It cannot triangularize any input.\n");
    exit(1);
  }
  if ((f = fopen(argv[1], "r")) == NULL) exit(1);
  
  if (mm_read_banner(f, &matcode) != 0)
  {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }
  symm =  mm_is_symmetric(matcode);
  if (!symm) symmcheck = false; // Default to false for non-symmetric
  else symmcheck = true;	// Default to true for symmetric graphs

  if (argc >= 3){
    // Read the optional base name:
    snprintf(basename, 256, "%s", argv[2]);	
    // Then read in the user-provided symmetry override, if present:
    if ( argc == 4 ) { // Check for symmetry flag
      symmcheck = (bool)(atoi(argv[3]) > 0);
    }
  }
  else
  {
    snprintf(basename, 256, "binary_graph");	
  }


  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
      mm_is_sparse(matcode) )
  {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f,&M, &N, &nz)) !=0)
    exit(1);

  uint64_t num_vertices = MAX(M, N);
  
  if (symm && symmcheck){
    std::cout << "Matrix is symmetric. Duplicating all edges each way..." << std::endl;
  }
  if (!symm && symmcheck) {
    std::cout << "WARNING: Matrix is NOT symmetric. You have elected to symmetrize it." << std::endl;
  }

  if (symm && !symmcheck) {
    std::cout << "WARNING: Matrix is symmetric (likely lower trianglular as is the standard for matrix market format). \n"
      "However, you have elected NOT to symmetrize it." << std::endl;
  }



  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  std::vector<std::vector<j_val>> CSR_MATRIX(num_vertices);
  
  VTYPE idx, jdx;
  WTYPE val_w;
  for (i=0; i<nz; i++)
  {
    // This line needs to be changed if the datatypes are modified!
    fscanf(f, "%u %u %u\n", &idx, &jdx, &val_w);
    idx--;
    jdx--;
    if (idx == jdx) continue; // Skip self edges.

    CSR_MATRIX[idx].push_back((j_val){jdx, val_w});
    if (symmcheck){
      CSR_MATRIX[jdx].push_back((j_val){idx, val_w});
    }
  }
  if (symmcheck && symm){
    std::cout << "WARNING: all edges duplicated due to symmetric flag in mtx input. Make sure this is what you wanted!\n";
  }
  if (symmcheck && !symm){
    std::cout << "WARNING: all edges duplicated due to user override. MTX file did not indicate symmetry." << std::endl;
  }

  if (f !=stdin) fclose(f);



  /********************************************************************/
  /* Convert multimap representation to arrays (IA, JA, VA) in memory */
  /********************************************************************/
  nz *= 2; // Assume everything doubled (this is the worst case when symmetrizing).
  uint64_t sz_ia = num_vertices + 1;
  uint64_t csum = 0;
  uint64_t jsum = 0;
  // First allocate memory:
  std::vector<VTYPE> IA(sz_ia);
  std::vector<j_val> JVA(nz);

  
  std::cout << "Moving 2D vectors into CSR vectors..." << std::endl;
  // Store edges into CSR arrays
  for ( uint64_t v_idx = 0; v_idx < num_vertices; v_idx++){
    // Check if v_idx is the source node for any edges:
    auto n_vidx = CSR_MATRIX[v_idx].size();
    // If there are some sort & copy them to the JVA array
    if ( n_vidx > 0){ 
      n_vidx = sort_and_dedup_vec( CSR_MATRIX[v_idx] );
      // Refresh begin and end:
      auto lb = CSR_MATRIX[v_idx].begin();
      auto ub = CSR_MATRIX[v_idx].end();
      for (auto st = lb; st != ub; st++) {
       JVA[jsum] = *st;
       jsum ++;
     }
   }
   IA[v_idx] = csum;
   csum += n_vidx;
 }
 IA[num_vertices] = csum;
 // Update number of edges actually in processed graph. 
 nz = csum; 

#ifdef DEBUG
 for ( uint64_t v_idx = 0; v_idx < num_vertices; v_idx++){
    auto n_vidx = CSR_MATRIX[v_idx].size();
    auto lb = CSR_MATRIX[v_idx].begin();
    auto ub = CSR_MATRIX[v_idx].end();
    for (auto st = lb; st != ub; st++) {
      uint32_t j_idx = st->first;
      printf("%u %u\n", v_idx, j_idx); 
   }
 }
 #endif



  /************************/
  /* now write out matrix */
  /************************/

  std::cout << "Opening IA, JA, and VA files..." << std::endl;
  std::ofstream foutIA;
  std::ofstream foutJA;
  std::ofstream foutVA;

  char IAfname[256];
  char JAfname[256];
  char VAfname[256];

  snprintf(IAfname, 256, "%s_ia.bin", basename);
  snprintf(JAfname, 256, "%s_ja.bin", basename);
  snprintf(VAfname, 256, "%s_va.bin", basename);

  foutIA.open(IAfname, std::ofstream::out | std::ofstream::binary );
  foutJA.open(JAfname, std::ofstream::out | std::ofstream::binary );
  foutVA.open(VAfname, std::ofstream::out | std::ofstream::binary );

  if ( !foutIA.is_open() || !foutJA.is_open() || !foutVA.is_open() ) {
    std::cerr << "ERROR: could not open IA, JA, VA files for writing." << std::endl;
    return 1;
  }

  // Sanity checks on final sizes, given that we store to disk with 32-bit
  // unsigned integers:
  assert(sz_ia < UINT_MAX);
  assert(nz < UINT_MAX);
  // If the above checks fail, then writing out as uint32_t won't work.

  std::cout << "Writing out to IA, JA, and VA files..." << std::endl;
  printf("Going to write %lu nodes and %lu edges...\n", num_vertices, nz);
  // First write out the size in entries of each file, and then file contents.
  // Write IA:
  VTYPE IA_WO_SZ = (VTYPE)sz_ia;
  foutIA.write(reinterpret_cast<const char*>(&IA_WO_SZ), sizeof(IA_WO_SZ));
  for ( uint64_t v_idx = 0; v_idx < sz_ia; v_idx++){
    VTYPE data = IA[v_idx];
    foutIA.write(reinterpret_cast<const char*>(&data), sizeof(data));
  }
  // And write JA:
  VTYPE JVA_WO_SZ = (VTYPE)nz;
  foutJA.write(reinterpret_cast<const char*>(&JVA_WO_SZ), sizeof(JVA_WO_SZ) );
  for ( uint64_t e_idx = 0; e_idx < nz; e_idx++){
    VTYPE data = JVA[e_idx].first;
    foutJA.write(reinterpret_cast<const char*>(&data), sizeof(data));
  }
  // and write VA:
  foutVA.write(reinterpret_cast<const char*>(&JVA_WO_SZ), sizeof(JVA_WO_SZ) );
  for ( uint64_t e_idx = 0; e_idx < nz; e_idx++){
    WTYPE data = JVA[e_idx].second;
    foutVA.write(reinterpret_cast<const char*>(&data), sizeof(data));
  }
  std::cout << "All done!\n" << std::endl;
  std::cout << "Wrote out " << sz_ia-1 << " vertices as "   << typeid(VTYPE).name() <<  std::endl;
  std::cout << "Wrote out " << nz <<      " edges as "      << typeid(VTYPE).name() << std::endl;
  std::cout << "Wrote out " << nz <<      " weights using " << typeid(WTYPE).name() << std::endl;


  foutIA.close();
  foutJA.close();
  foutVA.close();

  return 0;
}

