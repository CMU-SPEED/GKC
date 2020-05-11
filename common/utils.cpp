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
 
#include "utils.h"

unsigned long long rdtsc()
{
 unsigned a, d;

 __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

 return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

int cmp_u32(const void* p1, const void* p2){
  uint32_t v1 = *(uint32_t*)p1;
  uint32_t v2 = *(uint32_t*)p2;
  if (v1 < v2) return -1;
  if (v1 > v2) return 1;
  return 0;
}

#define BROOM_SIZE 1024*1024*30 // 30 MB
auto broom = std::vector<char>(BROOM_SIZE);
void clean_caches(){
#pragma omp parallel
  {
    for (uint32_t i = 0; i < BROOM_SIZE; i++)
    {
      broom[i] = broom[i] + 2;
    }
  }
}

char* truncate_fname(char* input_fname){
	char *output = (char*) malloc(256);
	uint32_t MAX_LEN = 255;
	uint32_t idx=0; 
	while (idx < 255 && input_fname[idx] != '\0') idx++;
	MAX_LEN = (idx-1 >= 0) ? idx-1 : 0;	
	uint32_t s_idx = 0, e_idx=MAX_LEN;
	// find start (find last / or take beginning):
	idx = 0;
	while (idx < MAX_LEN-1)
	{
		if (input_fname[idx] == '/') s_idx = idx+1;
		idx++;	
	}
	// Find end (find first '.' from end or end):
	idx = MAX_LEN;
	while (idx > 0) 
	{
		if (input_fname[idx] == '.') break;
		idx--;
	}
	e_idx = idx;
	strncpy(output, input_fname+s_idx, e_idx - s_idx);
 output[e_idx - s_idx] = '\0';
	return output;
}

