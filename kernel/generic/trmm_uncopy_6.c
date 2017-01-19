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

/* 2017/14/01 wernsaar@googlemail.com */

#include <stdio.h>
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

  BLASLONG i, js;
  BLASLONG X;

  FLOAT data01, data02, data03, data04, data05, data06, data07, data08;
  FLOAT data09, data10, data11, data12, data13, data14, data15, data16;
  FLOAT data17, data18, data19, data20, data21, data22, data23, data24;
  FLOAT data25, data26, data27, data28, data29, data30, data31, data32;
  FLOAT data33, data34, data35, data36;

  FLOAT *ao1, *ao2, *ao3, *ao4, *ao5, *ao6;

  BLASLONG ndiv = n/6;
  BLASLONG nmod = n%6;
  BLASLONG mdiv = m/6;
  BLASLONG mmod = m%6;

  js = ndiv;

  if (js > 0){
    do {
      X = posX;

      if (posX <= posY) {
	ao1 = a + posX + (posY + 0) * lda;
	ao2 = a + posX + (posY + 1) * lda;
	ao3 = a + posX + (posY + 2) * lda;
	ao4 = a + posX + (posY + 3) * lda;
	ao5 = a + posX + (posY + 4) * lda;
	ao6 = a + posX + (posY + 5) * lda;
      } else {
	ao1 = a + posY + (posX + 0) * lda;
	ao2 = a + posY + (posX + 1) * lda;
	ao3 = a + posY + (posX + 2) * lda;
	ao4 = a + posY + (posX + 3) * lda;
	ao5 = a + posY + (posX + 4) * lda;
	ao6 = a + posY + (posX + 5) * lda;
      }

      i = mdiv;
      if (i > 0) {
	do {
	  if (X < posY) {

	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);
	    data03 = *(ao1 + 2);
	    data04 = *(ao1 + 3);
	    data05 = *(ao1 + 4);
	    data06 = *(ao1 + 5);

	    data07 = *(ao2 + 0);
	    data08 = *(ao2 + 1);
	    data09 = *(ao2 + 2);
	    data10 = *(ao2 + 3);
	    data11 = *(ao2 + 4);
	    data12 = *(ao2 + 5);

	    data13 = *(ao3 + 0);
	    data14 = *(ao3 + 1);
	    data15 = *(ao3 + 2);
	    data16 = *(ao3 + 3);
	    data17 = *(ao3 + 4);
	    data18 = *(ao3 + 5);

	    data19 = *(ao4 + 0);
	    data20 = *(ao4 + 1);
	    data21 = *(ao4 + 2);
	    data22 = *(ao4 + 3);
	    data23 = *(ao4 + 4);
	    data24 = *(ao4 + 5);

	    data25 = *(ao5 + 0);
	    data26 = *(ao5 + 1);
	    data27 = *(ao5 + 2);
	    data28 = *(ao5 + 3);
	    data29 = *(ao5 + 4);
	    data30 = *(ao5 + 5);

	    data31 = *(ao6 + 0);
	    data32 = *(ao6 + 1);
	    data33 = *(ao6 + 2);
	    data34 = *(ao6 + 3);
	    data35 = *(ao6 + 4);
	    data36 = *(ao6 + 5);

	    b[ 0] = data01;
	    b[ 1] = data07;
	    b[ 2] = data13;
	    b[ 3] = data19;
	    b[ 4] = data25;
	    b[ 5] = data31;

	    b[ 6] = data02;
	    b[ 7] = data08;
	    b[ 8] = data14;
	    b[ 9] = data20;
	    b[10] = data26;
	    b[11] = data32;

	    b[12] = data03;
	    b[13] = data09;
	    b[14] = data15;
	    b[15] = data21;
	    b[16] = data27;
	    b[17] = data33;

	    b[18] = data04;
	    b[19] = data10;
	    b[20] = data16;
	    b[21] = data22;
	    b[22] = data28;
	    b[23] = data34;

	    b[24] = data05;
	    b[25] = data11;
	    b[26] = data17;
	    b[27] = data23;
	    b[28] = data29;
	    b[29] = data35;

	    b[30] = data06;
	    b[31] = data12;
	    b[32] = data18;
	    b[33] = data24;
	    b[34] = data30;
	    b[35] = data36;

	    ao1 += 6;
	    ao2 += 6;
	    ao3 += 6;
	    ao4 += 6;
	    ao5 += 6;
	    ao6 += 6;

	    b += 36;

	  } else
	    if (X > posY) {
	      ao1 += 6 * lda;
	      ao2 += 6 * lda;
	      ao3 += 6 * lda;
	      ao4 += 6 * lda;
	      ao5 += 6 * lda;
	      ao6 += 6 * lda;

	      b += 36;

	    } else {

#ifndef UNIT
	      data01 = *(ao1 + 0);
#endif

	      data07 = *(ao2 + 0);
#ifndef UNIT
	      data08 = *(ao2 + 1);
#endif

	      data13 = *(ao3 + 0);
	      data14 = *(ao3 + 1);
#ifndef UNIT
	      data15 = *(ao3 + 2);
#endif

	      data19 = *(ao4 + 0);
	      data20 = *(ao4 + 1);
	      data21 = *(ao4 + 2);
#ifndef UNIT
	      data22 = *(ao4 + 3);
#endif

	      data25 = *(ao5 + 0);
	      data26 = *(ao5 + 1);
	      data27 = *(ao5 + 2);
	      data28 = *(ao5 + 3);
#ifndef UNIT
	      data29 = *(ao5 + 4);
#endif

	      data31 = *(ao6 + 0);
	      data32 = *(ao6 + 1);
	      data33 = *(ao6 + 2);
	      data34 = *(ao6 + 3);
	      data35 = *(ao6 + 4);
#ifndef UNIT
	      data36 = *(ao6 + 5);
#endif


#ifdef UNIT
	      b[ 0] = ONE;
#else
	      b[ 0] = data01;
#endif
	      b[ 1] = data07;
	      b[ 2] = data13;
	      b[ 3] = data19;
	      b[ 4] = data25;
	      b[ 5] = data31;

	      b[ 6] = ZERO;
#ifdef UNIT
	      b[ 7] = ONE;
#else
	      b[ 7] = data08;
#endif
	      b[ 8] = data14;
	      b[ 9] = data20;
	      b[10] = data26;
	      b[11] = data32;

	      b[12] = ZERO;
	      b[13] = ZERO;
#ifdef UNIT
	      b[14] = ONE;
#else
	      b[14] = data15;
#endif
	      b[15] = data21;
	      b[16] = data27;
	      b[17] = data33;

	      b[18] = ZERO;
	      b[19] = ZERO;
	      b[20] = ZERO;
#ifdef UNIT
	      b[21] = ONE;
#else
	      b[21] = data22;
#endif
	      b[22] = data28;
	      b[23] = data34;

	      b[24] = ZERO;
	      b[25] = ZERO;
	      b[26] = ZERO;
	      b[27] = ZERO;
#ifdef UNIT
	      b[28] = ONE;
#else
	      b[28] = data29;
#endif
	      b[29] = data35;


	      b[30] = ZERO;
	      b[31] = ZERO;
	      b[32] = ZERO;
	      b[33] = ZERO;
	      b[34] = ZERO;
#ifdef UNIT
	      b[35] = ONE;
#else
	      b[35] = data36;
#endif

	      ao1 += 6 * lda;
	      ao2 += 6 * lda;
	      ao3 += 6 * lda;
	      ao4 += 6 * lda;
	      ao5 += 6 * lda;
	      ao6 += 6 * lda;

	      b += 36;
	    }

	  X += 6;
	  i --;
	} while (i > 0);
      }

      i = mmod;
      if (i) {

	if (X < posY) {

	  if (m & 4) {

	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);
	    data03 = *(ao1 + 2);
	    data04 = *(ao1 + 3);

	    data07 = *(ao2 + 0);
	    data08 = *(ao2 + 1);
	    data09 = *(ao2 + 2);
	    data10 = *(ao2 + 3);

	    data13 = *(ao3 + 0);
	    data14 = *(ao3 + 1);
	    data15 = *(ao3 + 2);
	    data16 = *(ao3 + 3);

	    data19 = *(ao4 + 0);
	    data20 = *(ao4 + 1);
	    data21 = *(ao4 + 2);
	    data22 = *(ao4 + 3);

	    data25 = *(ao5 + 0);
	    data26 = *(ao5 + 1);
	    data27 = *(ao5 + 2);
	    data28 = *(ao5 + 3);

	    data31 = *(ao6 + 0);
	    data32 = *(ao6 + 1);
	    data33 = *(ao6 + 2);
	    data34 = *(ao6 + 3);

	    b[ 0] = data01;
	    b[ 1] = data07;
	    b[ 2] = data13;
	    b[ 3] = data19;
	    b[ 4] = data25;
	    b[ 5] = data31;

	    b[ 6] = data02;
	    b[ 7] = data08;
	    b[ 8] = data14;
	    b[ 9] = data20;
	    b[10] = data26;
	    b[11] = data32;

	    b[12] = data03;
	    b[13] = data09;
	    b[14] = data15;
	    b[15] = data21;
	    b[16] = data27;
	    b[17] = data33;

	    b[18] = data04;
	    b[19] = data10;
	    b[20] = data16;
	    b[21] = data22;
	    b[22] = data28;
	    b[23] = data34;

	    ao1 += 4;
	    ao2 += 4;
	    ao3 += 4;
	    ao4 += 4;
	    ao5 += 4;
	    ao6 += 4;

	  b += 24;
	  }

	  if (mmod & 2) {
	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);

	    data07 = *(ao2 + 0);
	    data08 = *(ao2 + 1);

	    data13 = *(ao3 + 0);
	    data14 = *(ao3 + 1);

	    data19 = *(ao4 + 0);
	    data20 = *(ao4 + 1);

	    data25 = *(ao5 + 0);
	    data26 = *(ao5 + 1);

	    data31 = *(ao6 + 0);
	    data32 = *(ao6 + 1);

	    b[ 0] = data01;
	    b[ 1] = data07;
	    b[ 2] = data13;
	    b[ 3] = data19;
	    b[ 4] = data25;
	    b[ 5] = data31;

	    b[ 6] = data02;
	    b[ 7] = data08;
	    b[ 8] = data14;
	    b[ 9] = data20;
	    b[10] = data26;
	    b[11] = data32;

	    ao1 += 2;
	    ao2 += 2;
	    ao3 += 2;
	    ao4 += 2;
	    ao5 += 2;
	    ao6 += 2;

	    b += 12;
	  }

	  if (mmod & 1) {
	    data01 = *(ao1 + 0);
	    data07 = *(ao2 + 0);
	    data13 = *(ao3 + 0);
	    data19 = *(ao4 + 0);
	    data25 = *(ao5 + 0);
	    data31 = *(ao6 + 0);

	    b[ 0] = data01;
	    b[ 1] = data07;
	    b[ 2] = data13;
	    b[ 3] = data19;
	    b[ 4] = data25;
	    b[ 5] = data31;

	    b += 6;
	  }
	} else
	  if (X > posY) {
	    if (mmod & 4) {
	      ao1 += 4 * lda;
	      ao2 += 4 * lda;
	      ao3 += 4 * lda;
	      ao4 += 4 * lda;

	      b += 24;
	    }

	    if (mmod & 2) {
	      ao1 += 2 * lda;
	      b += 12;
	    }

	    if (mmod & 1) {
	      b += 6;
	    }
	  } else {

#ifndef UNIT
	    data01 = *(ao1 + 0);
#endif
	    data07 = *(ao2 + 0);
	    data13 = *(ao3 + 0);
	    data19 = *(ao4 + 0);
	    data25 = *(ao5 + 0);
	    data31 = *(ao6 + 0);

	    if (i >= 2) {
#ifndef UNIT
	      data08 = *(ao2 + 1);
#endif
	      data14 = *(ao3 + 1);
	      data20 = *(ao4 + 1);
	      data26 = *(ao5 + 1);
	      data32 = *(ao6 + 1);
	    }

	    if (i >= 3) {
#ifndef UNIT
	      data15 = *(ao3 + 2);
#endif
	      data21 = *(ao4 + 2);
	      data27 = *(ao5 + 2);
	      data33 = *(ao6 + 2);
	    }

	    if (i >= 4) {
#ifndef UNIT
	      data22 = *(ao4 + 3);
#endif
	      data28 = *(ao5 + 3);
	      data34 = *(ao6 + 3);
	    }

	    if (i >= 5) {
#ifndef UNIT
	      data29 = *(ao5 + 4);
#endif
	      data35 = *(ao6 + 4);
	    }


#ifdef UNIT
	    b[ 0] = ONE;
#else
	    b[ 0] = data01;
#endif
	    b[ 1] = data07;
	    b[ 2] = data13;
	    b[ 3] = data19;
	    b[ 4] = data25;
	    b[ 5] = data31;
	    b += 6;

	    if(i >= 2) {
	      b[ 0] = ZERO;
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = data08;
#endif
	      b[ 2] = data14;
	      b[ 3] = data20;
	      b[ 4] = data26;
	      b[ 5] = data32;
	      b += 6;
	    }

	    if (i >= 3) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = data15;
#endif
	      b[ 3] = data21;
	      b[ 4] = data27;
	      b[ 5] = data33;
	      b += 6;
	    }

	    if (i >= 4) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
#ifdef UNIT
	      b[ 3] = ONE;
#else
	      b[ 3] = data22;
#endif
	      b[ 4] = data28;
	      b[ 5] = data34;
	      b += 6;
	    }

	    if (i >= 5) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
	      b[ 2] = ZERO;
	      b[ 3] = ZERO;
#ifdef UNIT
	      b[ 4] = ONE;
#else
	      b[ 4] = data29;
#endif
	      b[ 5] = data35;
	      b += 6;
	    }
	  }
      }

      posY += 6;
      js --;
    } while (js > 0);
  } /* End of main loop */


  if (nmod & 4){
      X = posX;

      if (posX <= posY) {
	ao1 = a + posX + (posY + 0) * lda;
	ao2 = a + posX + (posY + 1) * lda;
	ao3 = a + posX + (posY + 2) * lda;
	ao4 = a + posX + (posY + 3) * lda;
      } else {
	ao1 = a + posY + (posX + 0) * lda;
	ao2 = a + posY + (posX + 1) * lda;
	ao3 = a + posY + (posX + 2) * lda;
	ao4 = a + posY + (posX + 3) * lda;
      }

      i = (m >> 2);
      if (i > 0) {
	do {
	  if (X < posY) {
	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);
	    data03 = *(ao1 + 2);
	    data04 = *(ao1 + 3);

	    data09 = *(ao2 + 0);
	    data10 = *(ao2 + 1);
	    data11 = *(ao2 + 2);
	    data12 = *(ao2 + 3);

	    data17 = *(ao3 + 0);
	    data18 = *(ao3 + 1);
	    data19 = *(ao3 + 2);
	    data20 = *(ao3 + 3);

	    data25 = *(ao4 + 0);
	    data26 = *(ao4 + 1);
	    data27 = *(ao4 + 2);
	    data28 = *(ao4 + 3);

	    b[ 0] = data01;
	    b[ 1] = data09;
	    b[ 2] = data17;
	    b[ 3] = data25;

	    b[ 4] = data02;
	    b[ 5] = data10;
	    b[ 6] = data18;
	    b[ 7] = data26;

	    b[ 8] = data03;
	    b[ 9] = data11;
	    b[10] = data19;
	    b[11] = data27;

	    b[12] = data04;
	    b[13] = data12;
	    b[14] = data20;
	    b[15] = data28;

	    ao1 += 4;
	    ao2 += 4;
	    ao3 += 4;
	    ao4 += 4;

	    b += 16;

	  } else
	    if (X > posY) {
	      ao1 += 4 * lda;
	      ao2 += 4 * lda;
	      ao3 += 4 * lda;
	      ao4 += 4 * lda;
	      b += 16;

	    } else {

#ifdef UNIT
	      data09 = *(ao2 + 0);

	      data17 = *(ao3 + 0);
	      data18 = *(ao3 + 1);

	      data25 = *(ao4 + 0);
	      data26 = *(ao4 + 1);
	      data27 = *(ao4 + 2);

	      b[ 0] = ONE;
	      b[ 1] = data09;
	      b[ 2] = data17;
	      b[ 3] = data25;

	      b[ 4] = ZERO;
	      b[ 5] = ONE;
	      b[ 6] = data18;
	      b[ 7] = data26;

	      b[ 8] = ZERO;
	      b[ 9] = ZERO;
	      b[10] = ONE;
	      b[11] = data27;

	      b[12] = ZERO;
	      b[13] = ZERO;
	      b[14] = ZERO;
	      b[15] = ONE;
#else
	      data01 = *(ao1 + 0);

	      data09 = *(ao2 + 0);
	      data10 = *(ao2 + 1);

	      data17 = *(ao3 + 0);
	      data18 = *(ao3 + 1);
	      data19 = *(ao3 + 2);

	      data25 = *(ao4 + 0);
	      data26 = *(ao4 + 1);
	      data27 = *(ao4 + 2);
	      data28 = *(ao4 + 3);

	      b[ 0] = data01;
	      b[ 1] = data09;
	      b[ 2] = data17;
	      b[ 3] = data25;

	      b[ 4] = ZERO;
	      b[ 5] = data10;
	      b[ 6] = data18;
	      b[ 7] = data26;

	      b[ 8] = ZERO;
	      b[ 9] = ZERO;
	      b[10] = data19;
	      b[11] = data27;

	      b[12] = ZERO;
	      b[13] = ZERO;
	      b[14] = ZERO;
	      b[15] = data28;
#endif
	      ao1 += 4 * lda;
	      ao2 += 4 * lda;
	      ao3 += 4 * lda;
	      ao4 += 4 * lda;

	      b += 16;
	    }

	  X += 4;
	  i --;
	} while (i > 0);
      }

      i = (m & 3);
      if (i) {

	if (X < posY) {

	  if (m & 2) {
	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);

	    data09 = *(ao2 + 0);
	    data10 = *(ao2 + 1);

	    data17 = *(ao3 + 0);
	    data18 = *(ao3 + 1);

	    data25 = *(ao4 + 0);
	    data26 = *(ao4 + 1);

	    b[ 0] = data01;
	    b[ 1] = data09;
	    b[ 2] = data17;
	    b[ 3] = data25;

	    b[ 4] = data02;
	    b[ 5] = data10;
	    b[ 6] = data18;
	    b[ 7] = data26;

	    ao1 += 2;
	    ao2 += 2;
	    ao3 += 2;
	    ao4 += 2;

	    b += 8;
	  }

	  if (m & 1) {
	    data01 = *(ao1 + 0);
	    data09 = *(ao2 + 0);
	    data17 = *(ao3 + 0);
	    data25 = *(ao4 + 0);

	    b[ 0] = data01;
	    b[ 1] = data09;
	    b[ 2] = data17;
	    b[ 3] = data25;

	    b += 4;
	  }
	} else
	  if (X > posY) {
	    if (m & 2) {
	      ao1 += 2 * lda;
	      b += 8;
	    }

	    if (m & 1) {
	      b += 4;
	    }
	  } else {

#ifndef UNIT
	    data01 = *(ao1 + 0);
#endif
	    data09 = *(ao2 + 0);
	    data17 = *(ao3 + 0);
	    data25 = *(ao4 + 0);

	    if (i >= 2) {
#ifndef UNIT
	      data10 = *(ao2 + 1);
#endif
	      data18 = *(ao3 + 1);
	      data26 = *(ao4 + 1);
	    }

	    if (i >= 3) {
#ifndef UNIT
	      data19 = *(ao3 + 2);
#endif
	      data27 = *(ao4 + 2);
	    }

#ifdef UNIT
	    b[ 0] = ONE;
#else
	    b[ 0] = data01;
#endif
	    b[ 1] = data09;
	    b[ 2] = data17;
	    b[ 3] = data25;
	    b += 4;

	    if(i >= 2) {
	      b[ 0] = ZERO;
#ifdef UNIT
	      b[ 1] = ONE;
#else
	      b[ 1] = data10;
#endif
	      b[ 2] = data18;
	      b[ 3] = data26;
	      b += 4;
	    }

	    if (i >= 3) {
	      b[ 0] = ZERO;
	      b[ 1] = ZERO;
#ifdef UNIT
	      b[ 2] = ONE;
#else
	      b[ 2] = data19;
#endif
	      b[ 3] = data27;
	      b += 4;
	    }
	  }
      }

      posY += 4;
  }

  if (nmod & 2){
      X = posX;

      if (posX <= posY) {
	ao1 = a + posX + (posY + 0) * lda;
	ao2 = a + posX + (posY + 1) * lda;
      } else {
	ao1 = a + posY + (posX + 0) * lda;
	ao2 = a + posY + (posX + 1) * lda;
      }

      i = (m >> 1);
      if (i > 0) {
	do {
	  if (X < posY) {
	    data01 = *(ao1 + 0);
	    data02 = *(ao1 + 1);

	    data09 = *(ao2 + 0);
	    data10 = *(ao2 + 1);

	    b[ 0] = data01;
	    b[ 1] = data09;
	    b[ 2] = data02;
	    b[ 3] = data10;

	    ao1 += 2;
	    ao2 += 2;
	    b += 4;

	  } else
	    if (X > posY) {
	      ao1 += 2 * lda;
	      ao2 += 2 * lda;
	      b += 4;

	    } else {

#ifdef UNIT
	      data09 = *(ao2 + 0);

	      b[ 0] = ONE;
	      b[ 1] = data09;
	      b[ 2] = ZERO;
	      b[ 3] = ONE;
#else
	      data01 = *(ao1 + 0);

	      data09 = *(ao2 + 0);
	      data10 = *(ao2 + 1);

	      b[ 0] = data01;
	      b[ 1] = data09;
	      b[ 2] = ZERO;
	      b[ 3] = data10;
#endif
	      ao1 += 2 * lda;
	      ao2 += 2 * lda;
	      b += 4;
	    }

	  X += 2;
	  i --;
	} while (i > 0);
      }

      if (m & 1) {

	if (X < posY) {
	  data01 = *(ao1 + 0);
	  data09 = *(ao2 + 0);

	  b[ 0] = data01;
	  b[ 1] = data09;
	  b += 2;
	} else
	  if (X > posY) {
	    b += 2;
	  } else {
#ifdef UNIT
	    data09 = *(ao2 + 0);
	    b[ 0] = ONE;
	    b[ 1] = data09;
#else
	    data01 = *(ao1 + 0);
	    data09 = *(ao2 + 0);
	    b[ 0] = data01;
	    b[ 1] = data09;
#endif
	    b += 2;
	  }
      }
      posY += 2;
  }

  if (nmod & 1){
      X = posX;

      if (posX <= posY) {
	ao1 = a + posX + (posY + 0) * lda;
      } else {
	ao1 = a + posY + (posX + 0) * lda;
      }

      i = m;
      if (m > 0) {
	do {
	  if (X < posY) {
	    data01 = *(ao1 + 0);
	    b[ 0] = data01;
	    ao1 += 1;
	    b += 1;
	  } else
	    if (X > posY) {
	      ao1 += lda;
	      b += 1;
	    } else {
#ifdef UNIT
	      b[ 0] = ONE;
#else
	      data01 = *(ao1 + 0);
	      b[ 0] = data01;
#endif
	      ao1 += lda;
	      b += 1;
	    }

	  X += 1;
	  i --;
	} while (i > 0);
      }
  }

  return 0;
}

