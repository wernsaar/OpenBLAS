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

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

  BLASLONG i, js, offset;

  FLOAT data01, data02, data03, data04, data05, data06, data07, data08;
  FLOAT *ao1, *ao2, *ao3, *ao4, *ao5, *ao6, *ao7, *ao8;

  BLASLONG ndiv = n/6;
  BLASLONG nmod = n%6;

  js = ndiv;
  while (js > 0){

    offset = posX - posY;

    if (offset >  0) ao1 = a + posX + 0 + posY * lda; else ao1 = a + posY + (posX + 0) * lda;
    if (offset > -1) ao2 = a + posX + 1 + posY * lda; else ao2 = a + posY + (posX + 1) * lda;
    if (offset > -2) ao3 = a + posX + 2 + posY * lda; else ao3 = a + posY + (posX + 2) * lda;
    if (offset > -3) ao4 = a + posX + 3 + posY * lda; else ao4 = a + posY + (posX + 3) * lda;
    if (offset > -4) ao5 = a + posX + 4 + posY * lda; else ao5 = a + posY + (posX + 4) * lda;
    if (offset > -5) ao6 = a + posX + 5 + posY * lda; else ao6 = a + posY + (posX + 5) * lda;

    i     = m;

    while (i > 0) {
      data01 = *(ao1 + 0);
      data02 = *(ao2 + 0);
      data03 = *(ao3 + 0);
      data04 = *(ao4 + 0);
      data05 = *(ao5 + 0);
      data06 = *(ao6 + 0);

      if (offset >   0) ao1 += lda; else ao1 ++;
      if (offset >  -1) ao2 += lda; else ao2 ++;
      if (offset >  -2) ao3 += lda; else ao3 ++;
      if (offset >  -3) ao4 += lda; else ao4 ++;
      if (offset >  -4) ao5 += lda; else ao5 ++;
      if (offset >  -5) ao6 += lda; else ao6 ++;

      b[ 0] = data01;
      b[ 1] = data02;
      b[ 2] = data03;
      b[ 3] = data04;
      b[ 4] = data05;
      b[ 5] = data06;

      b += 6;

      offset --;
      i --;
    }

    posX += 6;
    js --;
  }

  if (nmod & 4) {
    offset = posX - posY;

    if (offset >  0) ao1 = a + posX + 0 + posY * lda; else ao1 = a + posY + (posX + 0) * lda;
    if (offset > -1) ao2 = a + posX + 1 + posY * lda; else ao2 = a + posY + (posX + 1) * lda;
    if (offset > -2) ao3 = a + posX + 2 + posY * lda; else ao3 = a + posY + (posX + 2) * lda;
    if (offset > -3) ao4 = a + posX + 3 + posY * lda; else ao4 = a + posY + (posX + 3) * lda;

    i     = m;

    while (i > 0) {
      data01 = *(ao1 + 0);
      data02 = *(ao2 + 0);
      data03 = *(ao3 + 0);
      data04 = *(ao4 + 0);

      if (offset >   0) ao1 += lda; else ao1 ++;
      if (offset >  -1) ao2 += lda; else ao2 ++;
      if (offset >  -2) ao3 += lda; else ao3 ++;
      if (offset >  -3) ao4 += lda; else ao4 ++;

      b[ 0] = data01;
      b[ 1] = data02;
      b[ 2] = data03;
      b[ 3] = data04;

      b += 4;

      offset --;
      i --;
    }

    posX += 4;
  }

  if (nmod & 2) {
    offset = posX - posY;

    if (offset >  0) ao1 = a + posX + 0 + posY * lda; else ao1 = a + posY + (posX + 0) * lda;
    if (offset > -1) ao2 = a + posX + 1 + posY * lda; else ao2 = a + posY + (posX + 1) * lda;

    i     = m;

    while (i > 0) {
      data01 = *(ao1 + 0);
      data02 = *(ao2 + 0);

      if (offset >   0) ao1 += lda; else ao1 ++;
      if (offset >  -1) ao2 += lda; else ao2 ++;

      b[ 0] = data01;
      b[ 1] = data02;

      b += 2;

      offset --;
      i --;
    }

    posX += 2;
  }

  if (nmod & 1) {
    offset = posX - posY;

    if (offset >  0) ao1 = a + posX + 0 + posY * lda; else ao1 = a + posY + (posX + 0) * lda;

    i     = m;

    while (i > 0) {
      data01 = *(ao1 + 0);

      if (offset >   0) ao1 += lda; else ao1 ++;

      b[ 0] = data01;

      b ++;

      offset --;
      i --;
    }
  }

  return 0;
}

