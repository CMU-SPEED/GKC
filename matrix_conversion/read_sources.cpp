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
 
/* Mark Blanco, 8 Jan 2020, Speed Group, CMU
 * Program to read source vertices out of mtx files and make them 
 * 0-indexed rather than 1-indexed.
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

int main(int argc, char *argv[])
{
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  uint64_t M, N;   
  char basename[256];

  if (argc < 2 || argc > 3)
  {
    fprintf(stderr, "Usage: %s [martix-market-sources-filename] (optional:base-output-name)\n", argv[0]);
    exit(1);
  }
  else    
  { 
    if ((f = fopen(argv[1], "r")) == NULL) exit(1);

    if (argc == 3){
      snprintf(basename, 256, "%s", argv[2]);	
    }
    else
    {
      snprintf(basename, 256, "graph_vertices");	
    }
  }

  if (mm_read_banner(f, &matcode) != 0)
  {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }

  if ((ret_code = mm_read_mtx_array_size(f,&M, &N)) !=0)
    exit(1);

  // Note: the mtx files store vertex ids as doubles in some cases.
  // So they need to be read as doubles and converted.
  std::vector<uint32_t> vertices;
  uint32_t vidx = 0;
  double raw_num = 0;
  for (uint64_t idx = 0; idx < M; idx++){  
    fscanf(f, "%lg\n", &raw_num);
    vidx = (uint32_t)raw_num;
    std::cout << vidx << std::endl;
    vidx --; // From ones to zero based indexing
    vertices.push_back(vidx);
  }
  std::cout << "Registered " << M << " vertices" << std::endl;

  if (f !=stdin) fclose(f);

  /************************/
  /* now write out array */
  /************************/

  std::ofstream foutA;

  char Afname[256];

  snprintf(Afname, 256, "%s.sources", basename);

  foutA.open(Afname, std::ofstream::out);

  if ( !foutA.is_open() ) {
    std::cerr << "ERROR: could not open file for writing." << std::endl;
    return 1;
  }

  std::cout << "Writing out " << vertices.size() << " vertices." << std::endl;
  for ( auto itr = vertices.begin(); itr != vertices.end(); itr++){
    //foutA.write(reinterpret_cast<const char*>(&(*itr)), sizeof(*itr));
    foutA << *itr << std::endl;
  }

  foutA.close();

  return 0;
}

