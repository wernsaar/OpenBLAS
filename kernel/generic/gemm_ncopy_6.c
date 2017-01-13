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
/* 2017/13/01 wernsaar@googlemail.com */

#include <stdio.h>
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){
  BLASLONG i, j;

  FLOAT *aoffset;
  FLOAT *aoffset1, *aoffset2, *aoffset3, *aoffset4;
  FLOAT *aoffset5, *aoffset6;

  FLOAT *boffset;
  FLOAT ctemp01, ctemp02, ctemp03, ctemp04;
  FLOAT ctemp05, ctemp06, ctemp07, ctemp08;
  FLOAT ctemp09, ctemp10, ctemp11, ctemp12;
  FLOAT ctemp13, ctemp14, ctemp15, ctemp16;
  FLOAT ctemp17, ctemp18, ctemp19, ctemp20;
  FLOAT ctemp21, ctemp22, ctemp23, ctemp24;
  FLOAT ctemp25, ctemp26, ctemp27, ctemp28;
  FLOAT ctemp29, ctemp30, ctemp31, ctemp32;
  FLOAT ctemp33, ctemp34, ctemp35, ctemp36;

  BLASLONG ndiv = n/6;
  BLASLONG nmod = n%6;
  BLASLONG mdiv = m/6;
  BLASLONG mmod = m%6;

  aoffset = a;
  boffset = b;

  j = ndiv;
  if (j > 0){
    do{
      aoffset1  = aoffset;
      aoffset2  = aoffset1 + lda;
      aoffset3  = aoffset2 + lda;
      aoffset4  = aoffset3 + lda;
      aoffset5  = aoffset4 + lda;
      aoffset6  = aoffset5 + lda;
      aoffset += 6 * lda;

      i = mdiv;
      if (i > 0){
	do{
	  ctemp01 = *(aoffset1 +  0);
	  ctemp02 = *(aoffset1 +  1);
	  ctemp03 = *(aoffset1 +  2);
	  ctemp04 = *(aoffset1 +  3);
	  ctemp05 = *(aoffset1 +  4);
	  ctemp06 = *(aoffset1 +  5);

	  ctemp07 = *(aoffset2 +  0);
	  ctemp08 = *(aoffset2 +  1);
	  ctemp09 = *(aoffset2 +  2);
	  ctemp10 = *(aoffset2 +  3);
	  ctemp11 = *(aoffset2 +  4);
	  ctemp12 = *(aoffset2 +  5);

	  ctemp13 = *(aoffset3 +  0);
	  ctemp14 = *(aoffset3 +  1);
	  ctemp15 = *(aoffset3 +  2);
	  ctemp16 = *(aoffset3 +  3);
	  ctemp17 = *(aoffset3 +  4);
	  ctemp18 = *(aoffset3 +  5);

	  ctemp19 = *(aoffset4 +  0);
	  ctemp20 = *(aoffset4 +  1);
	  ctemp21 = *(aoffset4 +  2);
	  ctemp22 = *(aoffset4 +  3);
	  ctemp23 = *(aoffset4 +  4);
	  ctemp24 = *(aoffset4 +  5);

	  ctemp25 = *(aoffset5 +  0);
	  ctemp26 = *(aoffset5 +  1);
	  ctemp27 = *(aoffset5 +  2);
	  ctemp28 = *(aoffset5 +  3);
	  ctemp29 = *(aoffset5 +  4);
	  ctemp30 = *(aoffset5 +  5);

	  ctemp31 = *(aoffset6 +  0);
	  ctemp32 = *(aoffset6 +  1);
	  ctemp33 = *(aoffset6 +  2);
	  ctemp34 = *(aoffset6 +  3);
	  ctemp35 = *(aoffset6 +  4);
	  ctemp36 = *(aoffset6 +  5);


	  *(boffset +  0) = ctemp01;
	  *(boffset +  1) = ctemp07;
	  *(boffset +  2) = ctemp13;
	  *(boffset +  3) = ctemp19;
	  *(boffset +  4) = ctemp25;
	  *(boffset +  5) = ctemp31;

	  *(boffset +  6) = ctemp02;
	  *(boffset +  7) = ctemp08;
	  *(boffset +  8) = ctemp14;
	  *(boffset +  9) = ctemp20;
	  *(boffset + 10) = ctemp26;
	  *(boffset + 11) = ctemp32;

	  *(boffset + 12) = ctemp03;
	  *(boffset + 13) = ctemp09;
	  *(boffset + 14) = ctemp15;
	  *(boffset + 15) = ctemp21;
	  *(boffset + 16) = ctemp27;
	  *(boffset + 17) = ctemp33;

	  *(boffset + 18) = ctemp04;
	  *(boffset + 19) = ctemp10;
	  *(boffset + 20) = ctemp16;
	  *(boffset + 21) = ctemp22;
	  *(boffset + 22) = ctemp28;
	  *(boffset + 23) = ctemp34;

	  *(boffset + 24) = ctemp05;
	  *(boffset + 25) = ctemp11;
	  *(boffset + 26) = ctemp17;
	  *(boffset + 27) = ctemp23;
	  *(boffset + 28) = ctemp29;
	  *(boffset + 29) = ctemp35;


	  *(boffset + 30) = ctemp06;
	  *(boffset + 31) = ctemp12;
	  *(boffset + 32) = ctemp18;
	  *(boffset + 33) = ctemp24;
	  *(boffset + 34) = ctemp30;
	  *(boffset + 35) = ctemp36;


	  aoffset1 +=  6;
	  aoffset2 +=  6;
	  aoffset3 +=  6;
	  aoffset4 +=  6;
	  aoffset5 +=  6;
	  aoffset6 +=  6;
	  boffset  += 36;
	  i --;
	}while(i > 0);
      }

      i = mmod;
      if (i > 0){
	do{
	  ctemp01 = *(aoffset1 +  0);
	  ctemp02 = *(aoffset2 +  0);
	  ctemp03 = *(aoffset3 +  0);
	  ctemp04 = *(aoffset4 +  0);
	  ctemp05 = *(aoffset5 +  0);
	  ctemp06 = *(aoffset6 +  0);

	  *(boffset +  0) = ctemp01;
	  *(boffset +  1) = ctemp02;
	  *(boffset +  2) = ctemp03;
	  *(boffset +  3) = ctemp04;
	  *(boffset +  4) = ctemp05;
	  *(boffset +  5) = ctemp06;

	  aoffset1 ++;
	  aoffset2 ++;
	  aoffset3 ++;
	  aoffset4 ++;
	  aoffset5 ++;
	  aoffset6 ++;

	  boffset += 6;
	  i --;
	}while(i > 0);
      }
      j--;
    }while(j > 0);
  } /* end of if(j > 0) */

  if (nmod & 4){
    aoffset1  = aoffset;
    aoffset2  = aoffset1 + lda;
    aoffset3  = aoffset2 + lda;
    aoffset4  = aoffset3 + lda;
    aoffset += 4 * lda;

    i = (m >> 2);
    if (i > 0){
      do{
	ctemp01 = *(aoffset1 +  0);
	ctemp02 = *(aoffset1 +  1);
	ctemp03 = *(aoffset1 +  2);
	ctemp04 = *(aoffset1 +  3);

	ctemp05 = *(aoffset2 +  0);
	ctemp06 = *(aoffset2 +  1);
	ctemp07 = *(aoffset2 +  2);
	ctemp08 = *(aoffset2 +  3);

	ctemp09 = *(aoffset3 +  0);
	ctemp10 = *(aoffset3 +  1);
	ctemp11 = *(aoffset3 +  2);
	ctemp12 = *(aoffset3 +  3);

	ctemp13 = *(aoffset4 +  0);
	ctemp14 = *(aoffset4 +  1);
	ctemp15 = *(aoffset4 +  2);
	ctemp16 = *(aoffset4 +  3);

	*(boffset +  0) = ctemp01;
	*(boffset +  1) = ctemp05;
	*(boffset +  2) = ctemp09;
	*(boffset +  3) = ctemp13;

	*(boffset +  4) = ctemp02;
	*(boffset +  5) = ctemp06;
	*(boffset +  6) = ctemp10;
	*(boffset +  7) = ctemp14;

	*(boffset +  8) = ctemp03;
	*(boffset +  9) = ctemp07;
	*(boffset + 10) = ctemp11;
	*(boffset + 11) = ctemp15;

	*(boffset + 12) = ctemp04;
	*(boffset + 13) = ctemp08;
	*(boffset + 14) = ctemp12;
	*(boffset + 15) = ctemp16;

	aoffset1 +=  4;
	aoffset2 +=  4;
	aoffset3 +=  4;
	aoffset4 +=  4;
	boffset  +=  16;
	i --;
      }while(i > 0);
    }

    i = (m & 3);
    if (i > 0){
      do{
	ctemp01 = *(aoffset1 +  0);
	ctemp02 = *(aoffset2 +  0);
	ctemp03 = *(aoffset3 +  0);
	ctemp04 = *(aoffset4 +  0);

	*(boffset +  0) = ctemp01;
	*(boffset +  1) = ctemp02;
	*(boffset +  2) = ctemp03;
	*(boffset +  3) = ctemp04;

	aoffset1 ++;
	aoffset2 ++;
	aoffset3 ++;
	aoffset4 ++;

	boffset += 4;
	i --;
      }while(i > 0);
    }
  } /* end of if(j > 0) */

  if (nmod & 2){
    aoffset1  = aoffset;
    aoffset2  = aoffset1 + lda;
    aoffset += 2 * lda;

    i = (m >> 1);
    if (i > 0){
      do{
	ctemp01 = *(aoffset1 +  0);
	ctemp02 = *(aoffset1 +  1);
	ctemp03 = *(aoffset2 +  0);
	ctemp04 = *(aoffset2 +  1);

	*(boffset +  0) = ctemp01;
	*(boffset +  1) = ctemp03;
	*(boffset +  2) = ctemp02;
	*(boffset +  3) = ctemp04;

	aoffset1 +=  2;
	aoffset2 +=  2;
	boffset  +=  4;
	i --;
      }while(i > 0);
    }

    if (m & 1){
      ctemp01 = *(aoffset1 +  0);
      ctemp02 = *(aoffset2 +  0);

      *(boffset +  0) = ctemp01;
      *(boffset +  1) = ctemp02;

      aoffset1 ++;
      aoffset2 ++;
      boffset += 2;
    }
  } /* end of if(j > 0) */

  if (nmod & 1){
    aoffset1  = aoffset;

    i = m;
    if (i > 0){
      do{
	ctemp01 = *(aoffset1 +  0);

	*(boffset +  0) = ctemp01;

	aoffset1 ++;
	boffset  ++;
	i --;
      }while(i > 0);
    }

  } /* end of if(j > 0) */

  return 0;
}
