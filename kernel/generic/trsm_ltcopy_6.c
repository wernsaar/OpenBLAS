/***************************************************************************
Copyright (c) 2013, The OpenBLAS Project
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

/* 2017/15/01 wernsaar@googlemail.com */

#include <stdio.h>
#include "common.h"

#ifndef UNIT
#define INV(a) (ONE / (a))
#else
#define INV(a) (ONE)
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG offset, FLOAT *b){

  BLASLONG i, ii, j, jj;

  FLOAT data01, data02, data03, data04, data05, data06, data07, data08;
  FLOAT data09, data10, data11, data12, data13, data14, data15, data16;
  FLOAT data17, data18, data19, data20, data21, data22, data23, data24;
  FLOAT data25, data26, data27, data28, data29, data30, data31, data32;
  FLOAT data33, data34, data35, data36;

  FLOAT *a1, *a2, *a3, *a4, *a5, *a6;

  BLASLONG ndiv = n/6;
  BLASLONG nmod = n%6;
  BLASLONG mdiv = m/6;
  BLASLONG mmod = m%6;


  jj = offset;

  j = ndiv;
  while (j > 0){

    a1 = a + 0 * lda;
    a2 = a + 1 * lda;
    a3 = a + 2 * lda;
    a4 = a + 3 * lda;
    a5 = a + 4 * lda;
    a6 = a + 5 * lda;

    ii = 0;
    i = mdiv;
    while (i > 0) {

      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

#ifndef UNIT
	data08 = *(a2 + 1);
#endif
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

#ifndef UNIT
	data15 = *(a3 + 2);
#endif
	data16 = *(a3 + 3);
	data17 = *(a3 + 4);
	data18 = *(a3 + 5);

#ifndef UNIT
	data22 = *(a4 + 3);
#endif
	data23 = *(a4 + 4);
	data24 = *(a4 + 5);

#ifndef UNIT
	data29 = *(a5 + 4);
#endif
	data30 = *(a5 + 5);

#ifndef UNIT
	data36 = *(a6 + 5);
#endif


	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  7) = INV(data08);
	*(b +  8) = data09;
	*(b +  9) = data10;
	*(b + 10) = data11;
	*(b + 11) = data12;

	*(b + 14) = INV(data15);
	*(b + 15) = data16;
	*(b + 16) = data17;
	*(b + 17) = data18;

	*(b + 21) = INV(data22);
	*(b + 22) = data23;
	*(b + 23) = data24;

	*(b + 28) = INV(data29);
	*(b + 29) = data30;

	*(b + 35) = INV(data36);

      }

      if (ii < jj) {

	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

	data07 = *(a2 + 0);
	data08 = *(a2 + 1);
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

	data13 = *(a3 + 0);
	data14 = *(a3 + 1);
	data15 = *(a3 + 2);
	data16 = *(a3 + 3);
	data17 = *(a3 + 4);
	data18 = *(a3 + 5);

	data19 = *(a4 + 0);
	data20 = *(a4 + 1);
	data21 = *(a4 + 2);
	data22 = *(a4 + 3);
	data23 = *(a4 + 4);
	data24 = *(a4 + 5);

	data25 = *(a5 + 0);
	data26 = *(a5 + 1);
	data27 = *(a5 + 2);
	data28 = *(a5 + 3);
	data29 = *(a5 + 4);
	data30 = *(a5 + 5);

	data31 = *(a6 + 0);
	data32 = *(a6 + 1);
	data33 = *(a6 + 2);
	data34 = *(a6 + 3);
	data35 = *(a6 + 4);
	data36 = *(a6 + 5);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  6) = data07;
	*(b +  7) = data08;
	*(b +  8) = data09;
	*(b +  9) = data10;
	*(b + 10) = data11;
	*(b + 11) = data12;

	*(b + 12) = data13;
	*(b + 13) = data14;
	*(b + 14) = data15;
	*(b + 15) = data16;
	*(b + 16) = data17;
	*(b + 17) = data18;

	*(b + 18) = data19;
	*(b + 19) = data20;
	*(b + 20) = data21;
	*(b + 21) = data22;
	*(b + 22) = data23;
	*(b + 23) = data24;

	*(b + 24) = data25;
	*(b + 25) = data26;
	*(b + 26) = data27;
	*(b + 27) = data28;
	*(b + 28) = data29;
	*(b + 29) = data30;

	*(b + 30) = data31;
	*(b + 31) = data32;
	*(b + 32) = data33;
	*(b + 33) = data34;
	*(b + 34) = data35;
	*(b + 35) = data36;

      }

      a1 += 6 * lda;
      a2 += 6 * lda;
      a3 += 6 * lda;
      a4 += 6 * lda;
      a5 += 6 * lda;
      a6 += 6 * lda;
      b += 36;

      i  --;
      ii += 6;
    }

    if (mmod & 4) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

#ifndef UNIT
	data08 = *(a2 + 1);
#endif
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

#ifndef UNIT
	data15 = *(a3 + 2);
#endif
	data16 = *(a3 + 3);
	data17 = *(a3 + 4);
	data18 = *(a3 + 5);

#ifndef UNIT
	data22 = *(a4 + 3);
#endif
	data23 = *(a4 + 4);
	data24 = *(a4 + 5);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  7) = INV(data08);
	*(b +  8) = data09;
	*(b +  9) = data10;
	*(b + 10) = data11;
	*(b + 11) = data12;

	*(b + 14) = INV(data15);
	*(b + 15) = data16;
	*(b + 16) = data17;
	*(b + 17) = data18;

	*(b + 21) = INV(data22);
	*(b + 22) = data23;
	*(b + 23) = data24;

      }

      if (ii < jj) {

	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

	data07 = *(a2 + 0);
	data08 = *(a2 + 1);
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

	data13 = *(a3 + 0);
	data14 = *(a3 + 1);
	data15 = *(a3 + 2);
	data16 = *(a3 + 3);
	data17 = *(a3 + 4);
	data18 = *(a3 + 5);

	data19 = *(a4 + 0);
	data20 = *(a4 + 1);
	data21 = *(a4 + 2);
	data22 = *(a4 + 3);
	data23 = *(a4 + 4);
	data24 = *(a4 + 5);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  6) = data07;
	*(b +  7) = data08;
	*(b +  8) = data09;
	*(b +  9) = data10;
	*(b + 10) = data11;
	*(b + 11) = data12;

	*(b + 12) = data13;
	*(b + 13) = data14;
	*(b + 14) = data15;
	*(b + 15) = data16;
	*(b + 16) = data17;
	*(b + 17) = data18;

	*(b + 18) = data19;
	*(b + 19) = data20;
	*(b + 20) = data21;
	*(b + 21) = data22;
	*(b + 22) = data23;
	*(b + 23) = data24;

      }

      a1 += 4 * lda;
      a2 += 4 * lda;
      a3 += 4 * lda;
      a4 += 4 * lda;
      b += 24;

      ii += 4;
    }

    if (mmod & 2) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

#ifndef UNIT
	data08 = *(a2 + 1);
#endif
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  7) = INV(data08);
	*(b +  8) = data09;
	*(b +  9) = data12;
	*(b + 10) = data11;
	*(b + 11) = data12;
      }

      if (ii < jj) {

	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

	data07 = *(a2 + 0);
	data08 = *(a2 + 1);
	data09 = *(a2 + 2);
	data10 = *(a2 + 3);
	data11 = *(a2 + 4);
	data12 = *(a2 + 5);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;

	*(b +  6) = data07;
	*(b +  7) = data08;
	*(b +  8) = data09;
	*(b +  9) = data10;
	*(b + 10) = data11;
	*(b + 11) = data12;

      }

      a1 += 2 * lda;
      a2 += 2 * lda;
      b += 12;

      ii += 2;
    }

    if (mmod & 1) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);
	data05 = *(a1 + 4);
	data06 = *(a1 + 5);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data05;
	*(b +  5) = data06;
      }
      b += 6;
    }
    a += 6;
    jj += 6;
    j  --;
  }

  if (nmod & 4) {

    a1 = a + 0 * lda;
    a2 = a + 1 * lda;
    a3 = a + 2 * lda;
    a4 = a + 3 * lda;

    ii = 0;
    i = (m >> 2);
    while (i > 0) {

      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

#ifndef UNIT
	data10 = *(a2 + 1);
#endif
	data11 = *(a2 + 2);
	data12 = *(a2 + 3);


#ifndef UNIT
	data19 = *(a3 + 2);
#endif
	data20 = *(a3 + 3);


#ifndef UNIT
	data28 = *(a4 + 3);
#endif

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;

	*(b +  5) = INV(data10);
	*(b +  6) = data11;
	*(b +  7) = data12;

	*(b + 10) = INV(data19);
	*(b + 11) = data20;
	*(b + 15) = INV(data28);
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

	data09 = *(a2 + 0);
	data10 = *(a2 + 1);
	data11 = *(a2 + 2);
	data12 = *(a2 + 3);

	data17 = *(a3 + 0);
	data18 = *(a3 + 1);
	data19 = *(a3 + 2);
	data20 = *(a3 + 3);

	data25 = *(a4 + 0);
	data26 = *(a4 + 1);
	data27 = *(a4 + 2);
	data28 = *(a4 + 3);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data09;
	*(b +  5) = data10;
	*(b +  6) = data11;
	*(b +  7) = data12;

	*(b +  8) = data17;
	*(b +  9) = data18;
	*(b + 10) = data19;
	*(b + 11) = data20;
	*(b + 12) = data25;
	*(b + 13) = data26;
	*(b + 14) = data27;
	*(b + 15) = data28;
      }

      a1 += 4 * lda;
      a2 += 4 * lda;
      a3 += 4 * lda;
      a4 += 4 * lda;
      b += 16;

      i  --;
      ii += 4;
    }

    if (m & 2) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

#ifndef UNIT
	data10 = *(a2 + 1);
#endif
	data11 = *(a2 + 2);
	data12 = *(a2 + 3);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;

	*(b +  6) = INV(data10);
	*(b +  7) = data11;
	*(b +  8) = data12;
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

	data09 = *(a2 + 0);
	data10 = *(a2 + 1);
	data11 = *(a2 + 2);
	data12 = *(a2 + 3);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	*(b +  4) = data09;
	*(b +  5) = data10;
	*(b +  6) = data11;
	*(b +  7) = data12;
      }

      a1 += 2 * lda;
      a2 += 2 * lda;
      b += 8;
      ii += 2;
    }

    if (m & 1) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
	// *(b +  4) = data05;
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data03 = *(a1 + 2);
	data04 = *(a1 + 3);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data03;
	*(b +  3) = data04;
      }
      b += 4;
    }
    a += 4;
    jj += 4;
  }

  if (nmod & 2) {

    a1 = a + 0 * lda;
    a2 = a + 1 * lda;

    ii = 0;
    i = (m >> 1);
    while (i > 0) {

      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);

#ifndef UNIT
	data10 = *(a2 + 1);
#endif

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
	*(b +  3) = INV(data10);
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);
	data09 = *(a2 + 0);
	data10 = *(a2 + 1);

	*(b +  0) = data01;
	*(b +  1) = data02;
	*(b +  2) = data09;
	*(b +  3) = data10;
      }

      a1 += 2 * lda;
      a2 += 2 * lda;
      b += 4;

      i  --;
      ii += 2;
    }

    if (m & 1) {
      if (ii == jj) {

#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	data02 = *(a1 + 1);

	*(b +  0) = INV(data01);
	*(b +  1) = data02;
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	data02 = *(a1 + 1);

	*(b +  0) = data01;
	*(b +  1) = data02;
      }
      b += 2;
    }
    a += 2;
    jj += 2;
  }

  if (nmod & 1) {

    a1 = a + 0 * lda;

    ii = 0;
    i = m;
    while (i > 0) {

      if (ii == jj) {
#ifndef UNIT
	data01 = *(a1 + 0);
#endif
	*(b +  0) = INV(data01);
      }

      if (ii < jj) {
	data01 = *(a1 + 0);
	*(b +  0) = data01;
      }

      a1 += lda;
      b += 1;

      i  --;
      ii += 1;
    }

  }

  return 0;
}


