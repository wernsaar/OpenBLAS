/***************************************************************************
Copyright (c) 2013-2017, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* 2017/12/07 wernsaar@googlemail.com				     */

#include <stdio.h>
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){

  BLASLONG i, j;

  FLOAT *aoffset;
  FLOAT *aoffset1, *aoffset2, *aoffset3, *aoffset4;
  FLOAT *aoffset5, *aoffset6;

  FLOAT *boffset,  *boffset1, *boffset2, *boffset3, *boffset4;

  FLOAT ctemp01, ctemp02, ctemp03, ctemp04;
  FLOAT ctemp05, ctemp06, ctemp07, ctemp08;
  FLOAT ctemp09, ctemp10, ctemp11, ctemp12;
  FLOAT ctemp13, ctemp14, ctemp15, ctemp16;
  FLOAT ctemp17, ctemp18, ctemp19, ctemp20;
  FLOAT ctemp21, ctemp22, ctemp23, ctemp24;
  FLOAT ctemp25, ctemp26, ctemp27, ctemp28;
  FLOAT ctemp29, ctemp30, ctemp31, ctemp32;
  FLOAT ctemp33, ctemp34, ctemp35, ctemp36;

  aoffset   = a;
  boffset   = b;

  BLASLONG mdiv = m/6;
  BLASLONG mmod = m%6;
  BLASLONG ndiv = n/6;
  BLASLONG nmod = n%6;

#if 0
  fprintf(stderr, "M = %d N = %d\n", m, n);
#endif

  boffset2  = boffset  + m  * ndiv * 6;
  boffset3  = boffset2 + m  * (nmod & 4 );
  boffset4  = boffset3 + m  * (nmod & 2 );

  j = mdiv;
  if (j > 0){
    do{
      aoffset1  = aoffset;
      aoffset2  = aoffset1 + lda;
      aoffset3  = aoffset2 + lda;
      aoffset4  = aoffset3 + lda;
      aoffset5  = aoffset4 + lda;
      aoffset6  = aoffset5 + lda;
      aoffset += 6 * lda;

      boffset1  = boffset;
      boffset  += 36;

      i = ndiv;
      if (i > 0){
	do{

	  ctemp01 = *(aoffset1 + 0);
	  ctemp02 = *(aoffset1 + 1);
	  ctemp03 = *(aoffset1 + 2);
	  ctemp04 = *(aoffset1 + 3);
	  ctemp05 = *(aoffset1 + 4);
	  ctemp06 = *(aoffset1 + 5);
	  aoffset1 += 6;

	  ctemp07 = *(aoffset2 + 0);
	  ctemp08 = *(aoffset2 + 1);
	  ctemp09 = *(aoffset2 + 2);
	  ctemp10 = *(aoffset2 + 3);
	  ctemp11 = *(aoffset2 + 4);
	  ctemp12 = *(aoffset2 + 5);
	  aoffset2 += 6;

	  ctemp13 = *(aoffset3 + 0);
	  ctemp14 = *(aoffset3 + 1);
	  ctemp15 = *(aoffset3 + 2);
	  ctemp16 = *(aoffset3 + 3);
	  ctemp17 = *(aoffset3 + 4);
	  ctemp18 = *(aoffset3 + 5);
	  aoffset3 += 6;

	  ctemp19 = *(aoffset4 + 0);
	  ctemp20 = *(aoffset4 + 1);
	  ctemp21 = *(aoffset4 + 2);
	  ctemp22 = *(aoffset4 + 3);
	  ctemp23 = *(aoffset4 + 4);
	  ctemp24 = *(aoffset4 + 5);
	  aoffset4 += 6;

	  ctemp25 = *(aoffset5 + 0);
	  ctemp26 = *(aoffset5 + 1);
	  ctemp27 = *(aoffset5 + 2);
	  ctemp28 = *(aoffset5 + 3);
	  ctemp29 = *(aoffset5 + 4);
	  ctemp30 = *(aoffset5 + 5);
	  aoffset5 += 6;

	  ctemp31 = *(aoffset6 + 0);
	  ctemp32 = *(aoffset6 + 1);
	  ctemp33 = *(aoffset6 + 2);
	  ctemp34 = *(aoffset6 + 3);
	  ctemp35 = *(aoffset6 + 4);
	  ctemp36 = *(aoffset6 + 5);
	  aoffset6 += 6;

	  *(boffset1 +  0) = ctemp01;
	  *(boffset1 +  1) = ctemp02;
	  *(boffset1 +  2) = ctemp03;
	  *(boffset1 +  3) = ctemp04;
	  *(boffset1 +  4) = ctemp05;
	  *(boffset1 +  5) = ctemp06;

	  *(boffset1 +  6) = ctemp07;
	  *(boffset1 +  7) = ctemp08;
	  *(boffset1 +  8) = ctemp09;
	  *(boffset1 +  9) = ctemp10;
	  *(boffset1 + 10) = ctemp11;
	  *(boffset1 + 11) = ctemp12;

	  *(boffset1 + 12) = ctemp13;
	  *(boffset1 + 13) = ctemp14;
	  *(boffset1 + 14) = ctemp15;
	  *(boffset1 + 15) = ctemp16;
	  *(boffset1 + 16) = ctemp17;
	  *(boffset1 + 17) = ctemp18;

	  *(boffset1 + 18) = ctemp19;
	  *(boffset1 + 19) = ctemp20;
	  *(boffset1 + 20) = ctemp21;
	  *(boffset1 + 21) = ctemp22;
	  *(boffset1 + 22) = ctemp23;
	  *(boffset1 + 23) = ctemp24;

	  *(boffset1 + 24) = ctemp25;
	  *(boffset1 + 25) = ctemp26;
	  *(boffset1 + 26) = ctemp27;
	  *(boffset1 + 27) = ctemp28;
	  *(boffset1 + 28) = ctemp29;
	  *(boffset1 + 29) = ctemp30;

	  *(boffset1 + 30) = ctemp31;
	  *(boffset1 + 31) = ctemp32;
	  *(boffset1 + 32) = ctemp33;
	  *(boffset1 + 33) = ctemp34;
	  *(boffset1 + 34) = ctemp35;
	  *(boffset1 + 35) = ctemp36;


	  boffset1 += m * 6;
	  i --;
	}while(i > 0);
      }

      if (nmod & 4){
	ctemp01 = *(aoffset1 + 0);
	ctemp02 = *(aoffset1 + 1);
	ctemp03 = *(aoffset1 + 2);
	ctemp04 = *(aoffset1 + 3);
	aoffset1 += 4;

	ctemp05 = *(aoffset2 + 0);
	ctemp06 = *(aoffset2 + 1);
	ctemp07 = *(aoffset2 + 2);
	ctemp08 = *(aoffset2 + 3);
	aoffset2 += 4;

	ctemp09 = *(aoffset3 + 0);
	ctemp10 = *(aoffset3 + 1);
	ctemp11 = *(aoffset3 + 2);
	ctemp12 = *(aoffset3 + 3);
	aoffset3 += 4;

	ctemp13 = *(aoffset4 + 0);
	ctemp14 = *(aoffset4 + 1);
	ctemp15 = *(aoffset4 + 2);
	ctemp16 = *(aoffset4 + 3);
	aoffset4 += 4;

	ctemp17 = *(aoffset5 + 0);
	ctemp18 = *(aoffset5 + 1);
	ctemp19 = *(aoffset5 + 2);
	ctemp20 = *(aoffset5 + 3);
	aoffset5 += 4;

	ctemp21 = *(aoffset6 + 0);
	ctemp22 = *(aoffset6 + 1);
	ctemp23 = *(aoffset6 + 2);
	ctemp24 = *(aoffset6 + 3);
	aoffset6 += 4;

	*(boffset2 +  0) = ctemp01;
	*(boffset2 +  1) = ctemp02;
	*(boffset2 +  2) = ctemp03;
	*(boffset2 +  3) = ctemp04;

	*(boffset2 +  4) = ctemp05;
	*(boffset2 +  5) = ctemp06;
	*(boffset2 +  6) = ctemp07;
	*(boffset2 +  7) = ctemp08;

	*(boffset2 +  8) = ctemp09;
	*(boffset2 +  9) = ctemp10;
	*(boffset2 + 10) = ctemp11;
	*(boffset2 + 11) = ctemp12;

	*(boffset2 + 12) = ctemp13;
	*(boffset2 + 13) = ctemp14;
	*(boffset2 + 14) = ctemp15;
	*(boffset2 + 15) = ctemp16;

	*(boffset2 + 16) = ctemp17;
	*(boffset2 + 17) = ctemp18;
	*(boffset2 + 18) = ctemp19;
	*(boffset2 + 19) = ctemp20;

	*(boffset2 + 20) = ctemp21;
	*(boffset2 + 21) = ctemp22;
	*(boffset2 + 22) = ctemp23;
	*(boffset2 + 23) = ctemp24;

	boffset2 += 24;
      }

      if (nmod & 2){
	ctemp01 = *(aoffset1 + 0);
	ctemp02 = *(aoffset1 + 1);
	aoffset1 += 2;

	ctemp03 = *(aoffset2 + 0);
	ctemp04 = *(aoffset2 + 1);
	aoffset2 += 2;

	ctemp05 = *(aoffset3 + 0);
	ctemp06 = *(aoffset3 + 1);
	aoffset3 += 2;

	ctemp07 = *(aoffset4 + 0);
	ctemp08 = *(aoffset4 + 1);
	aoffset4 += 2;

	ctemp09 = *(aoffset5 + 0);
	ctemp10 = *(aoffset5 + 1);
	aoffset5 += 2;

	ctemp11 = *(aoffset6 + 0);
	ctemp12 = *(aoffset6 + 1);
	aoffset6 += 2;

	*(boffset3 +  0) = ctemp01;
	*(boffset3 +  1) = ctemp02;

	*(boffset3 +  2) = ctemp03;
	*(boffset3 +  3) = ctemp04;

	*(boffset3 +  4) = ctemp05;
	*(boffset3 +  5) = ctemp06;

	*(boffset3 +  6) = ctemp07;
	*(boffset3 +  7) = ctemp08;

	*(boffset3 +  8) = ctemp09;
	*(boffset3 +  9) = ctemp10;

	*(boffset3 + 10) = ctemp11;
	*(boffset3 + 11) = ctemp12;

	boffset3 += 12;
      }

      if (nmod & 1){
	ctemp01 = *(aoffset1 + 0);
	aoffset1 ++;
	ctemp02 = *(aoffset2 + 0);
	aoffset2 ++;
	ctemp03 = *(aoffset3 + 0);
	aoffset3 ++;
	ctemp04 = *(aoffset4 + 0);
	aoffset4 ++;
	ctemp05 = *(aoffset5 + 0);
	aoffset5 ++;
	ctemp06 = *(aoffset6 + 0);
	aoffset6 ++;

	*(boffset4 +  0) = ctemp01;
	*(boffset4 +  1) = ctemp02;
	*(boffset4 +  2) = ctemp03;
	*(boffset4 +  3) = ctemp04;
	*(boffset4 +  4) = ctemp05;
	*(boffset4 +  5) = ctemp06;
	boffset4 += 6;
      }

      j--;
    }while(j > 0);
  }

  if (mmod & 4){

    aoffset1  = aoffset;
    aoffset2  = aoffset1 + lda;
    aoffset3  = aoffset2 + lda;
    aoffset4  = aoffset3 + lda;
    aoffset += 4 * lda;

    boffset1  = boffset;
    boffset  += 24;

    i = ndiv;
    if (i > 0){

      do{

	  ctemp01 = *(aoffset1 + 0);
	  ctemp02 = *(aoffset1 + 1);
	  ctemp03 = *(aoffset1 + 2);
	  ctemp04 = *(aoffset1 + 3);
	  ctemp05 = *(aoffset1 + 4);
	  ctemp06 = *(aoffset1 + 5);
	  aoffset1 += 6;

	  ctemp07 = *(aoffset2 + 0);
	  ctemp08 = *(aoffset2 + 1);
	  ctemp09 = *(aoffset2 + 2);
	  ctemp10 = *(aoffset2 + 3);
	  ctemp11 = *(aoffset2 + 4);
	  ctemp12 = *(aoffset2 + 5);
	  aoffset2 += 6;

	  ctemp13 = *(aoffset3 + 0);
	  ctemp14 = *(aoffset3 + 1);
	  ctemp15 = *(aoffset3 + 2);
	  ctemp16 = *(aoffset3 + 3);
	  ctemp17 = *(aoffset3 + 4);
	  ctemp18 = *(aoffset3 + 5);
	  aoffset3 += 6;

	  ctemp19 = *(aoffset4 + 0);
	  ctemp20 = *(aoffset4 + 1);
	  ctemp21 = *(aoffset4 + 2);
	  ctemp22 = *(aoffset4 + 3);
	  ctemp23 = *(aoffset4 + 4);
	  ctemp24 = *(aoffset4 + 5);
	  aoffset4 += 6;

	*(boffset1 +  0) = ctemp01;
	*(boffset1 +  1) = ctemp02;
	*(boffset1 +  2) = ctemp03;
	*(boffset1 +  3) = ctemp04;
	*(boffset1 +  4) = ctemp05;
	*(boffset1 +  5) = ctemp06;

	*(boffset1 +  6) = ctemp07;
	*(boffset1 +  7) = ctemp08;
	*(boffset1 +  8) = ctemp09;
	*(boffset1 +  9) = ctemp10;
	*(boffset1 + 10) = ctemp11;
	*(boffset1 + 11) = ctemp12;

	*(boffset1 + 12) = ctemp13;
	*(boffset1 + 13) = ctemp14;
	*(boffset1 + 14) = ctemp15;
	*(boffset1 + 15) = ctemp16;
	*(boffset1 + 16) = ctemp17;
	*(boffset1 + 17) = ctemp18;

	*(boffset1 + 18) = ctemp19;
	*(boffset1 + 19) = ctemp20;
	*(boffset1 + 20) = ctemp21;
	*(boffset1 + 21) = ctemp22;
	*(boffset1 + 22) = ctemp23;
	*(boffset1 + 23) = ctemp24;

	boffset1 += 6 * m;
	i --;
      }while(i > 0);
    }

    if (nmod & 4) {
      ctemp01 = *(aoffset1 + 0);
      ctemp02 = *(aoffset1 + 1);
      ctemp03 = *(aoffset1 + 2);
      ctemp04 = *(aoffset1 + 3);
      aoffset1 += 4;

      ctemp05 = *(aoffset2 + 0);
      ctemp06 = *(aoffset2 + 1);
      ctemp07 = *(aoffset2 + 2);
      ctemp08 = *(aoffset2 + 3);
      aoffset2 += 4;

      ctemp09 = *(aoffset3 + 0);
      ctemp10 = *(aoffset3 + 1);
      ctemp11 = *(aoffset3 + 2);
      ctemp12 = *(aoffset3 + 3);
      aoffset3 += 4;

      ctemp13 = *(aoffset4 + 0);
      ctemp14 = *(aoffset4 + 1);
      ctemp15 = *(aoffset4 + 2);
      ctemp16 = *(aoffset4 + 3);
      aoffset4 += 4;

      *(boffset2 +  0) = ctemp01;
      *(boffset2 +  1) = ctemp02;
      *(boffset2 +  2) = ctemp03;
      *(boffset2 +  3) = ctemp04;
      *(boffset2 +  4) = ctemp05;
      *(boffset2 +  5) = ctemp06;
      *(boffset2 +  6) = ctemp07;
      *(boffset2 +  7) = ctemp08;

      *(boffset2 +  8) = ctemp09;
      *(boffset2 +  9) = ctemp10;
      *(boffset2 + 10) = ctemp11;
      *(boffset2 + 11) = ctemp12;
      *(boffset2 + 12) = ctemp13;
      *(boffset2 + 13) = ctemp14;
      *(boffset2 + 14) = ctemp15;
      *(boffset2 + 15) = ctemp16;
      boffset2 += 16;
    }

    if (nmod & 2){
      ctemp01 = *(aoffset1 + 0);
      ctemp02 = *(aoffset1 + 1);
      aoffset1 += 2;

      ctemp03 = *(aoffset2 + 0);
      ctemp04 = *(aoffset2 + 1);
      aoffset2 += 2;

      ctemp05 = *(aoffset3 + 0);
      ctemp06 = *(aoffset3 + 1);
      aoffset3 += 2;

      ctemp07 = *(aoffset4 + 0);
      ctemp08 = *(aoffset4 + 1);
      aoffset4 += 2;

      *(boffset3 +  0) = ctemp01;
      *(boffset3 +  1) = ctemp02;
      *(boffset3 +  2) = ctemp03;
      *(boffset3 +  3) = ctemp04;
      *(boffset3 +  4) = ctemp05;
      *(boffset3 +  5) = ctemp06;
      *(boffset3 +  6) = ctemp07;
      *(boffset3 +  7) = ctemp08;
      boffset3 += 8;
    }

    if (nmod & 1){
      ctemp01 = *(aoffset1 + 0);
      aoffset1 ++;
      ctemp02 = *(aoffset2 + 0);
      aoffset2 ++;
      ctemp03 = *(aoffset3 + 0);
      aoffset3 ++;
      ctemp04 = *(aoffset4 + 0);
      aoffset4 ++;

      *(boffset4 +  0) = ctemp01;
      *(boffset4 +  1) = ctemp02;
      *(boffset4 +  2) = ctemp03;
      *(boffset4 +  3) = ctemp04;
      boffset4 += 4;
    }
  }

  if (mmod & 2){
    aoffset1  = aoffset;
    aoffset2  = aoffset1 + lda;
    aoffset += 2 * lda;

    boffset1  = boffset;
    boffset  += 12;

    i = ndiv;
    if (i > 0){
      do{

	  ctemp01 = *(aoffset1 + 0);
	  ctemp02 = *(aoffset1 + 1);
	  ctemp03 = *(aoffset1 + 2);
	  ctemp04 = *(aoffset1 + 3);
	  ctemp05 = *(aoffset1 + 4);
	  ctemp06 = *(aoffset1 + 5);
	  aoffset1 += 6;

	  ctemp07 = *(aoffset2 + 0);
	  ctemp08 = *(aoffset2 + 1);
	  ctemp09 = *(aoffset2 + 2);
	  ctemp10 = *(aoffset2 + 3);
	  ctemp11 = *(aoffset2 + 4);
	  ctemp12 = *(aoffset2 + 5);
	  aoffset2 += 6;

	*(boffset1 +  0) = ctemp01;
	*(boffset1 +  1) = ctemp02;
	*(boffset1 +  2) = ctemp03;
	*(boffset1 +  3) = ctemp04;
	*(boffset1 +  4) = ctemp05;
	*(boffset1 +  5) = ctemp06;

	*(boffset1 +  6) = ctemp07;
	*(boffset1 +  7) = ctemp08;
	*(boffset1 +  8) = ctemp09;
	*(boffset1 +  9) = ctemp10;
	*(boffset1 + 10) = ctemp11;
	*(boffset1 + 11) = ctemp12;

	boffset1 += 6 * m;
	i --;
      }while(i > 0);
    }

    if (nmod & 4){
      ctemp01 = *(aoffset1 + 0);
      ctemp02 = *(aoffset1 + 1);
      ctemp03 = *(aoffset1 + 2);
      ctemp04 = *(aoffset1 + 3);
      aoffset1 += 4;

      ctemp05 = *(aoffset2 + 0);
      ctemp06 = *(aoffset2 + 1);
      ctemp07 = *(aoffset2 + 2);
      ctemp08 = *(aoffset2 + 3);
      aoffset2 += 4;

      *(boffset2 +  0) = ctemp01;
      *(boffset2 +  1) = ctemp02;
      *(boffset2 +  2) = ctemp03;
      *(boffset2 +  3) = ctemp04;
      *(boffset2 +  4) = ctemp05;
      *(boffset2 +  5) = ctemp06;
      *(boffset2 +  6) = ctemp07;
      *(boffset2 +  7) = ctemp08;
      boffset2 += 8;
    }

    if (nmod & 2){
      ctemp01 = *(aoffset1 + 0);
      ctemp02 = *(aoffset1 + 1);
      aoffset1 += 2;

      ctemp03 = *(aoffset2 + 0);
      ctemp04 = *(aoffset2 + 1);
      aoffset2 += 2;

      *(boffset3 +  0) = ctemp01;
      *(boffset3 +  1) = ctemp02;
      *(boffset3 +  2) = ctemp03;
      *(boffset3 +  3) = ctemp04;
      boffset3 += 4;
    }

    if (nmod & 1){
      ctemp01 = *(aoffset1 + 0);
      aoffset1 ++;
      ctemp02 = *(aoffset2 + 0);
      aoffset2 ++;

      *(boffset4 +  0) = ctemp01;
      *(boffset4 +  1) = ctemp02;
      boffset4 += 2;
    }
  }

  if (mmod & 1){
    aoffset1  = aoffset;
    aoffset += lda;

    boffset1  = boffset;
    boffset  += 6;

    i = ndiv;
    if (i > 0){
      do{
	ctemp01 = *(aoffset1 + 0);
	ctemp02 = *(aoffset1 + 1);
	ctemp03 = *(aoffset1 + 2);
	ctemp04 = *(aoffset1 + 3);
	ctemp05 = *(aoffset1 + 4);
	ctemp06 = *(aoffset1 + 5);
	aoffset1 += 6;

	*(boffset1 +  0) = ctemp01;
	*(boffset1 +  1) = ctemp02;
	*(boffset1 +  2) = ctemp03;
	*(boffset1 +  3) = ctemp04;
	*(boffset1 +  4) = ctemp05;
	*(boffset1 +  5) = ctemp06;

	boffset1 += 6 * m;
	 i --;
       }while(i > 0);
     }

     if (nmod & 4){
       ctemp01 = *(aoffset1 + 0);
       ctemp02 = *(aoffset1 + 1);
       ctemp03 = *(aoffset1 + 2);
       ctemp04 = *(aoffset1 + 3);
       aoffset1 += 4;

       *(boffset2 +  0) = ctemp01;
       *(boffset2 +  1) = ctemp02;
       *(boffset2 +  2) = ctemp03;
       *(boffset2 +  3) = ctemp04;
       boffset2 += 4;
     }

     if (nmod & 2){
       ctemp01 = *(aoffset1 + 0);
       ctemp02 = *(aoffset1 + 1);
       aoffset1 += 2;

       *(boffset3 +  0) = ctemp01;
       *(boffset3 +  1) = ctemp02;
       boffset3 += 2;
     }

     if (nmod & 1){
       ctemp01 = *(aoffset1 + 0);
       aoffset1 ++;
      *(boffset4 +  0) = ctemp01;
      boffset4 ++;
    }
  }

  return 0;
}
