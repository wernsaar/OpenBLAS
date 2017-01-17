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

/* 2017/17/01 wernsaar@googlemail.com */

#include "common.h"

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc ,BLASLONG offset)
{

   BLASLONG i,j,k;
   FLOAT *C0,*C1,*C2,*C3,*ptrba,*ptrbb;

   FLOAT res0_0;
   FLOAT res0_1;
   FLOAT res0_2;
   FLOAT res0_3;
   FLOAT res0_4;
   FLOAT res0_5;

   FLOAT res1_0;
   FLOAT res1_1;
   FLOAT res1_2;
   FLOAT res1_3;
   FLOAT res1_4;
   FLOAT res1_5;

   FLOAT res2_0;
   FLOAT res2_1;
   FLOAT res2_2;
   FLOAT res2_3;
   FLOAT res2_4;
   FLOAT res2_5;

   FLOAT res3_0;
   FLOAT res3_1;
   FLOAT res3_2;
   FLOAT res3_3;
   FLOAT res3_4;
   FLOAT res3_5;

   FLOAT a0;
   FLOAT a1;

   FLOAT b0;
   FLOAT b1;
   FLOAT b2;
   FLOAT b3;


   BLASLONG mdiv = bm/6;
   BLASLONG mmod = bm%6;


   for (j=0; j<bn/4; j+=1)
   {
        C0 = C;
        C1 = C0+ldc;
        C2 = C1+ldc;
        C3 = C2+ldc;

        ptrba = ba;

        for (i=0; i<mdiv; i+=1)
        {

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;
		res0_4 = 0;
		res0_5 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;
		res1_4 = 0;
		res1_5 = 0;

		res2_0 = 0;
		res2_1 = 0;
		res2_2 = 0;
		res2_3 = 0;
		res2_4 = 0;
		res2_5 = 0;

		res3_0 = 0;
		res3_1 = 0;
		res3_2 = 0;
		res3_3 = 0;
		res3_4 = 0;
		res3_5 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;
			res2_1 += a1*b2;
			res3_1 += a1*b3;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;
			res2_2 += a0*b2;
			res3_2 += a0*b3;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;
			res2_3 += a1*b2;
			res3_3 += a1*b3;

			a0 = ptrba[4];
			res0_4 += a0*b0;
			res1_4 += a0*b1;
			res2_4 += a0*b2;
			res3_4 += a0*b3;

			a1 = ptrba[5];
			res0_5 += a1*b0;
			res1_5 += a1*b1;
			res2_5 += a1*b2;
			res3_5 += a1*b3;

			ptrba = ptrba+6;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;
		res0_4 *= alpha;
		res0_5 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;
		res1_4 *= alpha;
		res1_5 *= alpha;

		res2_0 *= alpha;
		res2_1 *= alpha;
		res2_2 *= alpha;
		res2_3 *= alpha;
		res2_4 *= alpha;
		res2_5 *= alpha;

		res3_0 *= alpha;
		res3_1 *= alpha;
		res3_2 *= alpha;
		res3_3 *= alpha;
		res3_4 *= alpha;
		res3_5 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;
		C0[4] += res0_4;
		C0[5] += res0_5;

		C1[0] += res1_0;
		C1[1] += res1_1;
		C1[2] += res1_2;
		C1[3] += res1_3;
		C1[4] += res1_4;
		C1[5] += res1_5;

		C2[0] += res2_0;
		C2[1] += res2_1;
		C2[2] += res2_2;
		C2[3] += res2_3;
		C2[4] += res2_4;
		C2[5] += res2_5;

		C3[0] += res3_0;
		C3[1] += res3_1;
		C3[2] += res3_2;
		C3[3] += res3_3;
		C3[4] += res3_4;
		C3[5] += res3_5;


		C0 = C0+6;
		C1 = C1+6;
		C2 = C2+6;
		C3 = C3+6;
	}

	if ( mmod & 4 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;

		res2_0 = 0;
		res2_1 = 0;
		res2_2 = 0;
		res2_3 = 0;

		res3_0 = 0;
		res3_1 = 0;
		res3_2 = 0;
		res3_3 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;
			res2_1 += a1*b2;
			res3_1 += a1*b3;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;
			res2_2 += a0*b2;
			res3_2 += a0*b3;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;
			res2_3 += a1*b2;
			res3_3 += a1*b3;

			ptrba = ptrba+4;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;

		res2_0 *= alpha;
		res2_1 *= alpha;
		res2_2 *= alpha;
		res2_3 *= alpha;

		res3_0 *= alpha;
		res3_1 *= alpha;
		res3_2 *= alpha;
		res3_3 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;

		C1[0] += res1_0;
		C1[1] += res1_1;
		C1[2] += res1_2;
		C1[3] += res1_3;

		C2[0] += res2_0;
		C2[1] += res2_1;
		C2[2] += res2_2;
		C2[3] += res2_3;

		C3[0] += res3_0;
		C3[1] += res3_1;
		C3[2] += res3_2;
		C3[3] += res3_3;


		C0 = C0+4;
		C1 = C1+4;
		C2 = C2+4;
		C3 = C3+4;

	}

	if ( mmod & 2 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;

		res1_0 = 0;
		res1_1 = 0;

		res2_0 = 0;
		res2_0 = 0;

		res2_1 = 0;
		res2_1 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;
			res2_1 += a1*b2;
			res3_1 += a1*b3;

			ptrba = ptrba+2;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;

		res2_0 *= alpha;
		res2_0 *= alpha;

		res3_1 *= alpha;
		res3_1 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;

		C1[0] += res1_0;
		C1[1] += res1_1;

		C2[0] += res2_0;
		C2[0] += res2_0;

		C3[1] += res3_1;
		C3[1] += res3_1;


		C0 = C0+2;
		C1 = C1+2;
		C2 = C2+2;
		C3 = C3+2;

	}

	if ( mmod & 1 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res1_0 = 0;
		res2_0 = 0;
		res3_0 = 0;


		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			ptrba = ptrba+1;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res1_0 *= alpha;
		res2_0 *= alpha;
		res3_0 *= alpha;

		C0[0] += res0_0;
		C1[0] += res1_0;
		C2[0] += res2_0;
		C3[0] += res3_0;


		C0 = C0+1;
		C1 = C1+1;
		C2 = C2+1;
		C3 = C3+1;

	}

	k = bk << 2;
	bb = bb + k;
        i = (ldc<<2);
        C = C+i;
    }



   for (j=0; j<bn&2; j+=1)
   {
        C0 = C;
        C1 = C0+ldc;

        ptrba = ba;

        for (i=0; i<mdiv; i+=1)
        {

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;
		res0_4 = 0;
		res0_5 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;
		res1_4 = 0;
		res1_5 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;

			a0 = ptrba[4];
			res0_4 += a0*b0;
			res1_4 += a0*b1;

			a1 = ptrba[5];
			res0_5 += a1*b0;
			res1_5 += a1*b1;

			ptrba = ptrba+6;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;
		res0_4 *= alpha;
		res0_5 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;
		res1_4 *= alpha;
		res1_5 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;
		C0[4] += res0_4;
		C0[5] += res0_5;

		C1[0] += res1_0;
		C1[1] += res1_1;
		C1[2] += res1_2;
		C1[3] += res1_3;
		C1[4] += res1_4;
		C1[5] += res1_5;


		C0 = C0+6;
		C1 = C1+6;
	}

	if ( mmod & 4 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;


		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;

			ptrba = ptrba+4;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;

		C1[0] += res1_0;
		C1[1] += res1_1;
		C1[2] += res1_2;
		C1[3] += res1_3;


		C0 = C0+4;
		C1 = C1+4;

	}

	if ( mmod & 2 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;

		res1_0 = 0;
		res1_1 = 0;


		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			ptrba = ptrba+2;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;

		C1[0] += res1_0;
		C1[1] += res1_1;


		C0 = C0+2;
		C1 = C1+2;

	}

	if ( mmod & 1 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res1_0 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			ptrba = ptrba+1;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;

		res1_0 *= alpha;

		C0[0] += res0_0;

		C1[0] += res1_0;


		C0 = C0+1;
		C1 = C1+1;

	}


        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
    }



   for (j=0; j<(bn&1); j+=1)
   {
        C0 = C;

        ptrba = ba;

        for (i=0; i<mdiv; i+=1)
        {

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;
		res0_4 = 0;
		res0_5 = 0;


		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			a0 = ptrba[2];
			res0_2 += a0*b0;

			a1 = ptrba[3];
			res0_3 += a1*b0;

			a0 = ptrba[4];
			res0_4 += a0*b0;

			a1 = ptrba[5];
			res0_5 += a1*b0;


			ptrba = ptrba+6;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;
		res0_4 *= alpha;
		res0_5 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;
		C0[4] += res0_4;
		C0[5] += res0_5;

		C0 = C0+6;
	}

	if ( mmod & 4 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;


		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			a0 = ptrba[2];
			res0_2 += a0*b0;

			a1 = ptrba[3];
			res0_3 += a1*b0;

			ptrba = ptrba+4;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;

		C0 = C0+4;

	}

	if ( mmod & 2 )
	{

		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			ptrba = ptrba+2;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;

		C0 = C0+2;

	}

	if ( mmod & 1 )
	{

		ptrbb = bb;

		res0_0 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			ptrba = ptrba+1;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;

		C0[0] += res0_0;

		C0 = C0+1;

	}

        k = (bk<<0);
        bb = bb+k;
        C = C+ldc;
   }
   return 0;
}

