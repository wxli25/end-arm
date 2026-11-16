/*
 * END KEM -- Cortex-M4 assembly macros shared by the NTT kernels (ntt7_m4.S, ntt12289_m4.S).
 *
 * Register roles (fixed in every macro of this file):
 *   r0  pointer into the polynomial          r10 twiddle / scratch
 *   r1  twiddle pointer (scratch in gs4_col  r11 qa = q * 2^alpha
 *       and ict3_light)                      r12 q in the top half
 *   r2-r9 eight packed coefficient pairs     r14 scratch
 *   s0-s7   eight parked words               s23  final constant of the Gentleman-Sande inverse
 *   s8-s22  fifteen twiddles of four layers  s24  plant(1), used for in-register reductions
 *   (ict3_light and ict3_twist use other FPU registers, see there)
 *
 * Four layers are merged on 16 words w_0..w_15 that lie STEP bytes apart: three layers on the eight odd words,
 * three layers on the eight even words, and the layer that pairs (w_2i, w_2i+1) as one add/sub per pair, with
 * the other operand parked in s0-s7 (the method of the pqm4 ML-KEM code).
 */

.macro ct a0, a1
    ct_bfly \a0, \a1, r10, r14, r12, r11
.endm
.macro gs a0, a1
    gs_bfly \a0, \a1, r10, r14, r12, r11
.endm
.macro pl a
    doubleplant \a, r10, r14, r12, r11
.endm

.macro load8 off, stride
    ldr.w   r2, [r0, #(\off)]
    ldr.w   r3, [r0, #(\off+\stride)]
    ldr.w   r4, [r0, #(\off+2*\stride)]
    ldr.w   r5, [r0, #(\off+3*\stride)]
    ldr.w   r6, [r0, #(\off+4*\stride)]
    ldr.w   r7, [r0, #(\off+5*\stride)]
    ldr.w   r8, [r0, #(\off+6*\stride)]
    ldr.w   r9, [r0, #(\off+7*\stride)]
.endm
.macro store8 off, stride
    str.w   r2, [r0, #(\off)]
    str.w   r3, [r0, #(\off+\stride)]
    str.w   r4, [r0, #(\off+2*\stride)]
    str.w   r5, [r0, #(\off+3*\stride)]
    str.w   r6, [r0, #(\off+4*\stride)]
    str.w   r7, [r0, #(\off+5*\stride)]
    str.w   r8, [r0, #(\off+6*\stride)]
    str.w   r9, [r0, #(\off+7*\stride)]
.endm
.macro pl8
    pl      r2
    pl      r3
    pl      r4
    pl      r5
    pl      r6
    pl      r7
    pl      r8
    pl      r9
.endm

/* ---------------------------------------------------------------------------------------------------------
 * forward (Cooley-Tukey)
 * --------------------------------------------------------------------------------------------------------- */

/* three layers on r2-r9, twiddles s8 | s9 s10 | s11-s14 */
.macro ct3_fpu
    vmov    r10, s8
    ct      r2, r6
    ct      r3, r7
    ct      r4, r8
    ct      r5, r9
    vmov    r10, s9
    ct      r2, r4
    ct      r3, r5
    vmov    r10, s10
    ct      r6, r8
    ct      r7, r9
    vmov    r10, s11
    ct      r2, r3
    vmov    r10, s12
    ct      r4, r5
    vmov    r10, s13
    ct      r6, r7
    vmov    r10, s14
    ct      r8, r9
.endm

/* multiply an odd word by its last-layer twiddle and park it */
.macro ct_park c, stw, sd
    vmov    r10, \stw
    pl      \c
    vmov    \sd, \c
.endm

/* last of the four layers: even word c (register) +- parked odd word */
.macro ct_last c, sodd, off, STEP
    vmov    r10, \sodd
    ssub16  r14, \c, r10
    sadd16  \c, \c, r10
    str.w   \c, [r0, #(\off)]
    str.w   r14, [r0, #(\off+\STEP)]
.endm

/* four CT layers on the 16 words at r0, r0 + STEP, ..., r0 + 15*STEP; twiddles in s8-s22.
 * RED = 1: every word is first reduced to [-h, h] (multiplied by s24 = plant(1)). */
.macro ct4_col STEP, RED
    load8   \STEP, (2*\STEP)            /* odd words */
    .if \RED
    vmov    r10, s24
    pl8
    .endif
    ct3_fpu
    ct_park r2, s15, s0
    ct_park r3, s16, s1
    ct_park r4, s17, s2
    ct_park r5, s18, s3
    ct_park r6, s19, s4
    ct_park r7, s20, s5
    ct_park r8, s21, s6
    ct_park r9, s22, s7
    load8   0, (2*\STEP)                /* even words */
    .if \RED
    vmov    r10, s24
    pl8
    .endif
    ct3_fpu
    ct_last r2, s0, 0, \STEP
    ct_last r3, s1, (2*\STEP), \STEP
    ct_last r4, s2, (4*\STEP), \STEP
    ct_last r5, s3, (6*\STEP), \STEP
    ct_last r6, s4, (8*\STEP), \STEP
    ct_last r7, s5, (10*\STEP), \STEP
    ct_last r8, s6, (12*\STEP), \STEP
    ct_last r9, s7, (14*\STEP), \STEP
.endm

/* ---------------------------------------------------------------------------------------------------------
 * inverse of the 7-layer transforms: Cooley-Tukey butterflies on the bit-reversed input, then a twist
 * (the method of the pqm4 ML-KEM code)
 *
 * With the 128 leaf positions P, psi the 256th root of unity of the forward transform and omega = psi^2:
 *     a_P = psi^-P / 128 * sum_i A_i omega^(-brv7(i) P).
 * The sum is a decimation-in-time transform. Its layer with half-size m = 1, 2, ..., 64 pairs the positions
 * P and P + m (P mod 2m < m) with the twiddle omega^(-(64/m) (P mod m)), which depends on the position only.
 * So the layers m = 1, 2, 4, 8 use the same eight constants e^i (e = omega^-8) in every block of 16 positions,
 * and seven of the twelve butterflies of the layers m = 1, 2, 4 on eight positions have the twiddle 1.
 * --------------------------------------------------------------------------------------------------------- */

/* layers m = 1, 2, 4 on r2-r9 = eight consecutive positions; s8 = plant(1), s10 = e^2, s12 = e^4, s14 = e^6.
 * The butterflies with the twiddle 1 are a plain add/sub: nothing is reduced and the sums double.
 * CUT = 0: for |in| <= B the results are bounded by (8B, 2B+2h, 4B+h, 2B+2h, 8B, 2B+2h, 4B+h, 2B+2h).
 * CUT = 1: the butterfly (c0, c4) of the layer m = 4 multiplies by plant(1), so c0 and c4 are bounded by 4B+h.
 * r1, r10 and r14 are scratch. */
.macro ict3_light CUT
    sadd16  r14, r2, r3                 /* m = 1: c0 -> r14 */
    ssub16  r3, r2, r3
    sadd16  r1, r4, r5                  /*        c2 -> r1  */
    ssub16  r5, r4, r5
    sadd16  r2, r6, r7                  /*        c4 -> r2  */
    ssub16  r7, r6, r7
    sadd16  r4, r8, r9                  /*        c6 -> r4  */
    ssub16  r9, r8, r9
    sadd16  r8, r14, r1                 /* m = 2: c0 -> r8  */
    ssub16  r1, r14, r1                 /*        c2 -> r1  */
    sadd16  r6, r2, r4                  /*        c4 -> r6  */
    ssub16  r4, r2, r4                  /*        c6 -> r4  */
    vmov    r10, s12
    ct      r3, r5
    ct      r7, r9
    .if \CUT                            /* m = 4: c0 -> r2, c4 -> r6 */
    vmov    r10, s8
    ct_bfly_to r2, r6, r8, r6, r10, r14, r12, r11
    .else
    sadd16  r2, r8, r6
    ssub16  r6, r8, r6
    .endif
    vmov    r10, s10
    ct      r3, r7
    vmov    r10, s12
    ct_bfly_to r4, r8, r1, r4, r10, r14, r12, r11       /* c2 -> r4, c6 -> r8 */
    vmov    r10, s14
    ct      r5, r9
.endm

/* layers m = 16, 32, 64 and the twist on r2-r9 = the positions j, j + 16, ..., j + 112 of a column j.
 * Constants of the column in s0 | s1 s2 | s3-s6 (twiddles) and s7-s14 (twist psi^-P / 128).
 * Every layer adds at most h; the twist returns [-h, h]. CANON = 1: mapped from [-128, 128] to [0, 257). */
.macro ict3_twist CANON
    vmov    r10, s0
    ct      r2, r3
    ct      r4, r5
    ct      r6, r7
    ct      r8, r9
    vmov    r10, s1
    ct      r2, r4
    ct      r6, r8
    vmov    r10, s2
    ct      r3, r5
    ct      r7, r9
    vmov    r10, s3
    ct      r2, r6
    vmov    r10, s4
    ct      r3, r7
    vmov    r10, s5
    ct      r4, r8
    vmov    r10, s6
    ct      r5, r9
    vmov    r10, s7
    pl      r2
    vmov    r10, s8
    pl      r3
    vmov    r10, s9
    pl      r4
    vmov    r10, s10
    pl      r5
    vmov    r10, s11
    pl      r6
    vmov    r10, s12
    pl      r7
    vmov    r10, s13
    pl      r8
    vmov    r10, s14
    pl      r9
    .if \CANON
    mov.w   r10, #0x01010101
    canon_pair r2, r14, r10
    canon_pair r3, r14, r10
    canon_pair r4, r14, r10
    canon_pair r5, r14, r10
    canon_pair r6, r14, r10
    canon_pair r7, r14, r10
    canon_pair r8, r14, r10
    canon_pair r9, r14, r10
    .endif
.endm

/* ---------------------------------------------------------------------------------------------------------
 * inverse of the complete transform mod 12289 (Gentleman-Sande)
 * --------------------------------------------------------------------------------------------------------- */

/* first of the four layers on the words at off, off + STEP: the sum stays in c, the product is parked in sd */
.macro gs_first c, off, stw, sd, STEP
    ldr.w   \c, [r0, #(\off)]
    ldr.w   r1, [r0, #(\off+\STEP)]
    ssub16  r10, \c, r1
    sadd16  \c, \c, r1
    vmov    r1, \stw
    doubleplant r10, r1, r14, r12, r11
    vmov    \sd, r10
.endm

/* three layers on r2-r9, twiddles s16-s19 | s20 s21 | s22; the sums of the first layer (r2, r4, r6, r8) are
 * reduced with s24 = plant(1), and at the end the sums r2-r5 are multiplied by s23. */
.macro gs3_fpu
    vmov    r10, s16
    gs      r2, r3
    vmov    r10, s17
    gs      r4, r5
    vmov    r10, s18
    gs      r6, r7
    vmov    r10, s19
    gs      r8, r9
    vmov    r10, s24
    pl      r2
    pl      r4
    pl      r6
    pl      r8
    vmov    r10, s20
    gs      r2, r4
    gs      r3, r5
    vmov    r10, s21
    gs      r6, r8
    gs      r7, r9
    vmov    r10, s22
    gs      r2, r6
    gs      r3, r7
    gs      r4, r8
    gs      r5, r9
    vmov    r10, s23
    pl      r2
    pl      r3
    pl      r4
    pl      r5
.endm

/* four GS layers on the 16 words at r0, r0 + STEP, ..., r0 + 15*STEP; constants in s8-s23 (r1 is clobbered).
 * REDSUM = 1: the eight sums of the first layer are reduced with s24 before the other three layers. */
.macro gs4_col STEP, REDSUM
    gs_first r2, 0, s8, s0, \STEP
    gs_first r3, (2*\STEP), s9, s1, \STEP
    gs_first r4, (4*\STEP), s10, s2, \STEP
    gs_first r5, (6*\STEP), s11, s3, \STEP
    gs_first r6, (8*\STEP), s12, s4, \STEP
    gs_first r7, (10*\STEP), s13, s5, \STEP
    gs_first r8, (12*\STEP), s14, s6, \STEP
    gs_first r9, (14*\STEP), s15, s7, \STEP
    .if \REDSUM
    vmov    r10, s24
    pl8
    .endif
    gs3_fpu
    store8  0, (2*\STEP)
    vmov    r2, s0
    vmov    r3, s1
    vmov    r4, s2
    vmov    r5, s3
    vmov    r6, s4
    vmov    r7, s5
    vmov    r8, s6
    vmov    r9, s7
    gs3_fpu
    store8  \STEP, (2*\STEP)
.endm
