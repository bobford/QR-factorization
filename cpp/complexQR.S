
/*

        Copyright © 2021 Robert A. Ford. All rights reserved.

        No further distribution is authorized without the expressed written consent of the copyright holder.

        IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
        INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE COPYRIGHT HOLDER HAS BEEN ADVISED
        OF THE POSSIBILITY OF SUCH DAMAGE.
        THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
        AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".
        THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


        QR decomposition of a complex matrix using Householder technique in Armv8.2-A assembly language.

        see:    http://www.seas.ucla.edu/~vandenbe/133A/lectures/qr.pdf
        and     https://codereview.stackexchange.com/questions/120978/householder-transformation
        and     http://arith.cs.ucla.edu/publications/House-Asil06.pdf implementing the method for the Householder vector

        Execution time ~ 3600 microseconds for this routine (R only, not including Q) for 192x120 matrix of random numbers on [-1.0, +1.0]
		approximately 60% efficience on 192x120 matrix
        Execution time ~ 18 milliseconds including .cpp calling program for adding computation of Q.
        Execution time ~ 40.3 milliseconds for this routine (R only, not including Q) for 384x240 matrix (222 millisec for both)

.
    Calling sequence:
    long result = complexQR((double) A, (double) v, (double) Q, (int) number of rows in A, (int) number of columns in A, (int) flag);
	result is always zero
	
    Rows of matrices must be contiguous: A[][] is stored as a one dimensional array: A[], likewise Q

    r0  = (double) A[]                                  size = 2 * m * n; "times 2" is for the complex data
    r1  = (double) v[]                                  size = 2 * m
    r2  = (double) Q[]                                  size = 2 * m * m; Q must be set to the identity matrix prior to calling this routine.
    r3  = (int) number of rows in A
    r4  = (int) number of columns in A
    r5  = (int) flag: 0 => do not compute Q, 1 => compute Q
	
	complex numbers are stored as real and imaginary in sequential locations


    scratch registers

        r6  = pointer to current value of A[][]
        r7  = working register for v[], = r1
        r8  = scratch
        r9  = scratch
        r10  = scratch
        r11 = offset into A[][]
        r12 = counter
        r13 = counter
        r14 = counter, outer loop, used for addressing
        r15 = counter, outer loop, counts down

        d0-d7

	This does not use the SIMD registers.

*/


.text

#include "C:/Users/bob/StudioProjects/ComplexQR/app/src/main/cpp/registers.h"       // this file (complexQR.S, not registers.h!) must have "S" suffix, not "s"


//      name the registers

loop    .req x15
incr    .req x11
addr    .req x1
v       .req x7
A       .req x6
rows    .req x3
cols    .req x4
count   .req x10

        .global complexQR
        .p2align 4
        .type    complexQR, %function

complexQR:
//        save_registers                            // not needed

        mov x14, cols                               // loop over all columns
        mov loop, cols
        cmp rows, cols                              // rows == cols? Is matrix square?
        bne lp
        sub loop, loop, #1                          // square matrix does not to do the last element
lp:
        mov x7, x1                                  // use v pointer
        fsub d2, d2, d2                             // clear register for sum
        fsub d3, d3, d3                             // clear register for sum
        sub x10, x4, x14                            // column index, k = 0, 1, 2, ... n-1

        mul x11, x10, x4                            // k*n
        add x11, x11, x10                           // k*n + k = offset to A[k][k] or Q[k][k]
        add x11, x11, x11                           // *2 for complex numbers

        mov x6, x11                                 // and a copy of where it came from

        add x8, x0, x11, lsl #3                     // increment to real part of complex matrix, both A an Q
        ldp d0, d1, [x8]                            // d0 = real A[k][k] and d1 = imag A[k][k]
        fmov d6, d0                                 // keep a copy of A[0][0] for later
        fmov d7, d1
        sub x12, x3, x10                            // m-k, this is the loop counter

vHv:                                                // this loop computes x Hermitian times x, where x = subcolumn of A
        fmadd d2, d0, d0, d2                        // d2 = d0 * d0 + d2;   sum += A[i][k]*A[i][k] real part
        fmadd d2, d1, d1, d2                        // d2 = d1 * d1 + d2;   sum += A[i][k]*A[i][k] imag part
        add x11, x11, x4, lsl #1                    // next row for complex elements
        add x8, x0, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]

        subs x12, x12, #1                           // decrement loop counter
        cbnz x12, vHv                               // and continue loop if not zero. Result left in d2


        fsqrt d3, d2                                // sqrt of xH * x = norm (x)

        fsub d4, d4, d4                             // zero
        fsub d4, d4, d3                             // -sqrt
        fsub d5, d5, d5                             // zero
        fcmp d6, d5                                 // compare real (A0) with zero
        fcsel d5, d4, d3, lt                        // double sqrSum = sign*Math.sqrt(sum) = sign*norm(x)

        fsub d2, d2, d2                             // clear accumulator
        fmadd d2, d6, d6, d2                        // add re*re
        fmadd d2, d7, d7, d2                        // add im*im
        fsqrt d2, d2                                // sqrt( re*re + im*im)
        fdiv d3, d3, d2                             // norm(x) / sqrt()
        fmadd d0, d6, d3, d6
        fmadd d1, d7, d3, d7

        stp d0, d1, [x7], #16                       // debug v0

        fsub d2, d2, d2                             // clear register for sum
        fmadd d2, d0, d0, d2                        // start norm of v
        fmadd d2, d1, d1, d2

        sub x12, x3, x10                            // m-k
        sub x12, x12, #1                            // m-k-1 = m-(k+1)

        mov x11, x6                                 // restore copy of first element of A[][]
        add x11, x11, x4, lsl #1                    // offset to A[k+1][k]

vnorm:
        add x8, x0, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]
        fmadd d2, d0, d0, d2
        fmadd d2, d1, d1, d2
        stp d0, d1, [x7], #16                       // v[k]
        add x11, x11, x4, lsl #1

        subs x12, x12, #1                           // decrement loop counter
        cbnz x12, vnorm                             // and continue loop if not zero

        fsqrt d3, d2
		fmov d4, #0x70								// d4 = 1.0
		fdiv d3, d4, d3								// 1. / d3 = reciprocal of d3 so a multiply can be used in the loop

        mov x11, x6                                 // restore copy of first element of A[][]
        sub x12, x3, x10                            // m-k, restore the loop counter for loop on rows
        mov x7, x1                                  // use v pointer
norm:
        ldp d0, d1, [x7]
        fmul d0, d0, d3								// replaces original: fdiv d0, d0, d3
        fmul d1, d1, d3								// replaces original: fdiv d1, d1, d3
        stp d0, d1, [x7], #16                       // v[k]

        subs x12, x12, #1                           // decrement loop counter
        cbnz x12, norm                              // and continue loop if not zero


//      have Householder vector, now apply it

        sub x12, x4, x10                            // loop over the rest of the columns

house:
        fsub d4, d4, d4                             // clear register for sum real
        fsub d5, d5, d5                             // clear register for sum imag

        sub x13, x3, x10                            // m-k, the rest of the rows
        mov x7, x1                                  // load v[]
        mov x11, x6                                 // restore copy of first element of A[][]

sumR:                                               // this computes the complex scalar v Hermitian times x, where x = subcolumn of A
        add x8, x0, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]                            // d0 = real A[k][k] and d1 = imag A[k][k]
        ldp d2, d3, [x7], #16                       // load v

        fmadd d4, d0, d2, d4                        // d4 = d0 * d2 + d4;   sum real +=  Real(v[i]) * Real(A[i][j])
        fmadd d4, d1, d3, d4                        // d4 = d1 * d3 + d4;   sum real +=  Imag(v[i]) * Imag(A[i][j])
        fmsub d5, d0, d3, d5                        // d5 = -d0 * d3 + d5;  sum imag -=  Imag(v[i]) * Real(A[i][j])
        fmadd d5, d1, d2, d5                        // d5 = d1 * d2 + d5;   sum imag +=  Real(v[i]) * Imag(A[i][j])
        add x11, x11, x4, lsl #1                    // next row

        subs x13, x13, #1                           // decrement loop counter
        cbnz x13, sumR                              // and continue loop if not zero

        fadd d4, d4, d4                             // 2*sum
        fadd d5, d5, d5                             // 2*sum

        sub x13, x3, x10                            // m-k
        mov x7, x1
        mov x11, x6                                 // first element of A[][]

Aij:
        add x8, x0, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]                            // d0 = real A[k][k] and d1 = imag A[k][k]
        ldp d2, d3, [x7], #16                       // load v

        fmsub d0, d4, d2, d0                        // d0 = d0 - d4 * d2;   real A[i][j] -= 2 * Real(sum) * Real(v[i])
        fmadd d0, d5, d3, d0                        // d0 = d0 + d5 * d3;   real A[i][j] += 2 * Imag(sum) * Imag(v[i])
        fmsub d1, d4, d3, d1                        // d0 = d0 - d4 * d2;   imag A[i][j] -= 2 * Real(sum) * Imag(v[i])
        fmsub d1, d5, d2, d1                        // d0 = d0 - d5 * d3;   imag A[i][j] -= 2 * Imag(sum) * Real(v[i])

        stp d0, d1, [x8]                            // back to same place
        add x11, x11, x4, lsl #1                    // next row

        subs x13, x13, #1                           // decrement loop counter
        cbnz x13, Aij                               // and continue loop if not zero

        add x6, x6, #2                              // next column, complex data

        subs x12, x12, #1                           // decrement loop counter
        cbnz x12, house                             // and continue loop if not zero


//      now get the Q matrix, if requested.  Q Conjugate Transpose times Q is the identity matrix

        cbz x5, noQ

        mov x12, x3                                 // loop over all columns of Q: Q is mxm

        eor x6, x6, x6                              // starting offset for columns of Q[][]: x2 + x6
        mul x9, x10, x3                             // starting offset for rows of Q[][]: k*m

houseQ:
        fsub d4, d4, d4                             // clear register for sum real
        fsub d5, d5, d5                             // clear register for sum imag

        sub x13, x3, x10                            // m-k
        mov x7, x1                                  // load v[]
        add x11, x9, x6                             // j*m + j = offset to Q[j]
        add x11, x11, x11                           // complex

sumQ:
        add x8, x2, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]                            // d0 = real Q[k][k] and d1 = imag Q[k][k]
        ldp d2, d3, [x7], #16                       // load v

        fmadd d4, d0, d2, d4                        // d4 = d0 * d2 + d4;   sum real += Real(v[i]) * Real(Q[i][j])
        fmadd d4, d1, d3, d4                        // d4 = d1 * d3 + d4;   sum real += Imag(v[i]) * Imag(Q[i][j])
        fmsub d5, d0, d3, d5                        // d5 = -d0 * d3 + d5;   sum imag -= Imag(v[i]) * Real(Q[i][j])
        fmadd d5, d1, d2, d5                        // d5 = d1 * d2 + d5;   sum imag += Real(v[i]) * Imag(Q[i][j])

        add x11, x11, x3, lsl #1                    // next row x3, not x4 because Q is mxm

        subs x13, x13, #1                           // decrement loop counter
        cbnz x13, sumQ                              // and continue loop if not zero

        fadd d4, d4, d4                             // 2*sum real
        fadd d5, d5, d5                             // 2*sum imag

        sub x13, x3, x10                            // m-k
        mov x7, x1                                  // load v[]
        add x11, x9, x6                             // j*m + j = offset to Q[j]
        add x11, x11, x11                           // complex


Qij:
        add x8, x2, x11, lsl #3                     // increment to real part of complex
        ldp d0, d1, [x8]                            // d0 = real Q[k][k] and d1 = imag Q[k][k]
        ldp d2, d3, [x7], #16                       // load v

        fmsub d0, d4, d2, d0                        // d0 = d0 - d4 * d2;   real A[i][j] -= 2 * Real(sum) * Real(v[i])
        fmadd d0, d5, d3, d0                        // d0 = d0 + d5 * d3;   real A[i][j] += 2 * Imag(sum) * Imag(v[i])
        fmsub d1, d4, d3, d1                        // d0 = d0 - d4 * d2;   imag A[i][j] -= 2 * Real(sum) * Imag(v[i])
        fmsub d1, d5, d2, d1                        // d0 = d0 - d5 * d3;   imag A[i][j] -= 2 * Imag(sum) * Real(v[i])

        stp d0, d1, [x8]                            // back to same place
        add x11, x11, x3, lsl #1                    // next row

        subs x13, x13, #1                           // decrement loop counter
        cbnz x13, Qij                               // and continue loop if not zero

        add x6, x6, #1                              // next column, complex entries

        subs x12, x12, #1                           // decrement loop counter
        cbnz x12, houseQ                            // and continue loop if not zero

noQ:
        subs x14, x14, #1                           // decrement index counter
        subs loop, loop, #1                         // this is the loop counter
        cbnz loop, lp                               // and continue loop if not zero

        eor x0, x0, x0

//        restore_registers

        ret

