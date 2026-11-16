/*
 * END KEM -- Cortex-M4 kernel: leaf products of the 7-layer NTT domains (mod 3329 and mod 257), Thumb-2 + DSP.
 *
 *   r = a * b  leaf by leaf, each leaf in Z_q[x]/(x^L - zeta_leaf), L = END_LEAF (4 for END-512, 8 for END-1024).
 *
 * Method (as in the pqm4 ML-KEM code, extended to leaves of 4 and 8 coefficients):
 *   1. b is scaled once per leaf with Plantard multiplications: bc = b*c and zb = zeta*b*c, c = -2^32 mod q.
 *      Both are reduced to [-h, h], h = (q-1)/2, whatever int16 values b holds.
 *   2. every output coefficient is one dot product of a with bc / zb, accumulated in 32 bits with the dual
 *      16-bit multiply-accumulate instructions (SMUAD, SMUADX, SMLAD, SMLADX);
 *   3. one Plantard reduction per coefficient: (sum) * (-2^-32) cancels the factor c.
 * Inputs may be ANY int16 values (in particular the unreduced output of the forward NTT): a dot product is
 * bounded by L * 32768 * h <= 4.4e8, inside the valid interval of the reduction (plant.i).
 * Output: [-h, h] for q = 3329, canonical [0, 257) for q = 257.  r may alias a or b (a leaf is read before it is written).
 *
 * void name(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N], const int32_t tab[129]);
 *   tab = end_plant_base<q>: plant(c), then plant(+-zeta_leaf * c) for the 128 leaves.
 */
#include "asm_defs.h"
#ifdef END_M4_ASM

    .syntax unified
    .thumb
#include "plant.i"

/* register roles: r0 r, r1 a, r2 b, r3 leaf constants, r11 qa, r12 q (top half), r14 q^-1 mod 2^32; r4-r10 work */

.macro red a
    plant_red \a, r14, r12, r11
.endm

#if END_S == 2
/* ---------------------------------------------------------------------------------------------------------
 * leaves of 4 coefficients:  A01 = [a0|a1], A23 = [a2|a3] (low half first), likewise bc and zb
 *   r0 = a0 bc0 + a1 zb3 + a2 zb2 + a3 zb1        r1 = a0 bc1 + a1 bc0 + a2 zb3 + a3 zb2
 *   r2 = a0 bc2 + a1 bc1 + a2 bc0 + a3 zb3        r3 = a0 bc3 + a1 bc2 + a2 bc1 + a3 bc0
 * --------------------------------------------------------------------------------------------------------- */
.macro basemul name, Q, QA, QINV, CANON
    .global \name
    .type   \name, %function
    .align  2
\name:
    push    {r4-r11, lr}
    ldr.w   r14, [r3], #4               /* plant(c) */
    vmov    s0, r14
    add.w   r14, r3, #512               /* end of the 128 leaf constants */
    vmov    s1, r14
    movw    r11, #\QA
    movw    r12, #0
    movt    r12, #\Q
    movw    r14, #(\QINV & 0xFFFF)
    movt    r14, #(\QINV >> 16)
1:
    ldr.w   r6, [r2], #4                /* [b0|b1] */
    ldr.w   r7, [r2], #4                /* [b2|b3] */
    ldr.w   r4, [r3], #4                /* plant(zeta c) */
    smulwb  r5, r4, r7                  /* r8 = [zb2|zb3] */
    smulwt  r8, r4, r7
    smlabt  r5, r5, r12, r11
    smlabt  r8, r8, r12, r11
    pkhtb   r8, r8, r5, asr #16
    smulwt  r9, r4, r6                  /* r9.hi = zb1 */
    smlabt  r9, r9, r12, r11
    vmov    r4, s0
    doubleplant r6, r4, r5, r12, r11    /* r6 = [bc0|bc1] */
    doubleplant r7, r4, r5, r12, r11    /* r7 = [bc2|bc3] */
    ldr.w   r4, [r1], #4                /* [a0|a1] */
    ldr.w   r5, [r1], #4                /* [a2|a3] */
    pkhbt   r10, r6, r8                 /* [bc0|zb3] */
    pkhbt   r9, r8, r9                  /* [zb2|zb1] */
    smuad   r10, r4, r10
    smlad   r10, r5, r9, r10
    red     r10                         /* r0 */
    smuadx  r9, r4, r6
    smladx  r9, r5, r8, r9
    red     r9                          /* r1 */
    pkhtb   r10, r9, r10, asr #16       /* [r0|r1] */
    pkhbt   r8, r6, r8                  /* [bc0|zb3] */
    pkhbt   r9, r7, r6                  /* [bc2|bc1] */
    smuad   r9, r4, r9
    smlad   r9, r5, r8, r9
    red     r9                          /* r2 */
    smuadx  r8, r4, r7
    smladx  r8, r5, r6, r8
    red     r8                          /* r3 */
    pkhtb   r9, r8, r9, asr #16         /* [r2|r3] */
    .if \CANON
    mov.w   r4, #0x01010101
    canon_pair r10, r5, r4
    canon_pair r9, r5, r4
    .endif
    str.w   r10, [r0], #4
    str.w   r9, [r0], #4
    vmov    r4, s1
    cmp     r3, r4
    bne.w   1b
    pop     {r4-r11, pc}
    .size   \name, .-\name
.endm

    .text
    basemul END_SYM(m4_basemul3329), 3329, 26632, 0x6ba8f301, 0
    basemul END_SYM(m4_basemul257), 257, 16448, 0xff00ff01, 1
#endif /* END_S == 2 */

#if END_S == 4
/* ---------------------------------------------------------------------------------------------------------
 * leaves of 8 coefficients.  With e[j] = bc_j (j >= 0), e[j] = zb_{j+8} (j < 0):
 *     r_k = sum_{m=0..3}  a_{2m} e[k-2m] + a_{2m+1} e[k-2m-1]
 *   odd  k = 2t+1:  sum_m SMUADX(A_m, E_{t-m}),   E_j = [e_{2j} | e_{2j+1}]   (the natural words of zb and bc)
 *   even k = 2t  :  sum_m SMUAD (A_m, P_{t-m}),   P_j = [e_{2j} | e_{2j-1}] = PKHBT(E_j, E_{j-1})
 * The eight E words (E_-4..E_3 = zb, bc) and the seven P words (P_-3..P_3) of a leaf live in a 64-byte stack frame;
 * the a, b and constant pointers are parked in s2-s4 so that the dot products have A_0..A_3, four operands and
 * two accumulators in registers.
 *   frame:  sp+0..28  E_-4 .. E_3        sp+32..56  P_-3 .. P_3
 * --------------------------------------------------------------------------------------------------------- */
.macro dots8 t, CANON
    ldr.w   r5, [sp, #(32+4*\t)]        /* P_{t-3} .. P_t */
    ldr.w   r6, [sp, #(36+4*\t)]
    ldr.w   r7, [sp, #(40+4*\t)]
    ldr.w   r8, [sp, #(44+4*\t)]
    smuad   r9, r1, r8
    smlad   r9, r2, r7, r9
    smlad   r9, r3, r6, r9
    smlad   r9, r4, r5, r9
    red     r9                          /* r_{2t} */
    ldr.w   r5, [sp, #(4+4*\t)]         /* E_{t-3} .. E_t */
    ldr.w   r6, [sp, #(8+4*\t)]
    ldr.w   r7, [sp, #(12+4*\t)]
    ldr.w   r8, [sp, #(16+4*\t)]
    smuadx  r10, r1, r8
    smladx  r10, r2, r7, r10
    smladx  r10, r3, r6, r10
    smladx  r10, r4, r5, r10
    red     r10                         /* r_{2t+1} */
    pkhtb   r9, r10, r9, asr #16
    .if \CANON
    mov.w   r5, #0x01010101
    canon_pair r9, r6, r5
    .endif
    str.w   r9, [r0], #4
.endm

.macro basemul name, Q, QA, QINV, CANON
    .global \name
    .type   \name, %function
    .align  2
\name:
    push    {r4-r11, lr}
    sub     sp, sp, #64
    ldr.w   r14, [r3], #4               /* plant(c) */
    vmov    s0, r14
    add.w   r14, r3, #512               /* end of the 128 leaf constants */
    vmov    s1, r14
    vmov    s2, r1
    vmov    s3, r2
    vmov    s4, r3
    movw    r11, #\QA
    movw    r12, #0
    movt    r12, #\Q
    movw    r14, #(\QINV & 0xFFFF)
    movt    r14, #(\QINV >> 16)
1:
    vmov    r9, s3
    ldm     r9!, {r5-r8}                /* b0..b7 */
    vmov    s3, r9
    vmov    r10, s4
    ldr.w   r9, [r10], #4               /* plant(zeta c) */
    vmov    s4, r10
    smulwb  r10, r9, r5                 /* r1-r4 = zb */
    smulwt  r1, r9, r5
    smlabt  r10, r10, r12, r11
    smlabt  r1, r1, r12, r11
    pkhtb   r1, r1, r10, asr #16
    smulwb  r10, r9, r6
    smulwt  r2, r9, r6
    smlabt  r10, r10, r12, r11
    smlabt  r2, r2, r12, r11
    pkhtb   r2, r2, r10, asr #16
    smulwb  r10, r9, r7
    smulwt  r3, r9, r7
    smlabt  r10, r10, r12, r11
    smlabt  r3, r3, r12, r11
    pkhtb   r3, r3, r10, asr #16
    smulwb  r10, r9, r8
    smulwt  r4, r9, r8
    smlabt  r10, r10, r12, r11
    smlabt  r4, r4, r12, r11
    pkhtb   r4, r4, r10, asr #16
    str.w   r1, [sp, #0]
    str.w   r2, [sp, #4]
    str.w   r3, [sp, #8]
    str.w   r4, [sp, #12]
    pkhbt   r10, r2, r1                 /* P_-3 = [zb2|zb1] */
    str.w   r10, [sp, #32]
    pkhbt   r10, r3, r2                 /* P_-2 = [zb4|zb3] */
    str.w   r10, [sp, #36]
    pkhbt   r10, r4, r3                 /* P_-1 = [zb6|zb5] */
    str.w   r10, [sp, #40]
    vmov    r9, s0
    doubleplant r5, r9, r10, r12, r11   /* r5-r8 = bc */
    doubleplant r6, r9, r10, r12, r11
    doubleplant r7, r9, r10, r12, r11
    doubleplant r8, r9, r10, r12, r11
    pkhbt   r10, r5, r4                 /* P_0 = [bc0|zb7] */
    str.w   r10, [sp, #44]
    pkhbt   r10, r6, r5                 /* P_1 = [bc2|bc1] */
    str.w   r10, [sp, #48]
    pkhbt   r10, r7, r6                 /* P_2 = [bc4|bc3] */
    str.w   r10, [sp, #52]
    pkhbt   r10, r8, r7                 /* P_3 = [bc6|bc5] */
    str.w   r10, [sp, #56]
    str.w   r5, [sp, #16]
    str.w   r6, [sp, #20]
    str.w   r7, [sp, #24]
    str.w   r8, [sp, #28]
    vmov    r9, s2
    ldm     r9!, {r1-r4}                /* a0..a7 */
    vmov    s2, r9
    dots8   0, \CANON
    dots8   1, \CANON
    dots8   2, \CANON
    dots8   3, \CANON
    vmov    r9, s4
    vmov    r10, s1
    cmp     r9, r10
    bne.w   1b
    add     sp, sp, #64
    pop     {r4-r11, pc}
    .size   \name, .-\name
.endm

    .text
    basemul END_SYM(m4_basemul3329), 3329, 26632, 0x6ba8f301, 0
    basemul END_SYM(m4_basemul257), 257, 16448, 0xff00ff01, 1
#endif /* END_S == 4 */

#endif /* END_M4_ASM */
