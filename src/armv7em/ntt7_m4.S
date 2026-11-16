/*
 * END KEM -- Cortex-M4 kernel: the 7-layer NTTs mod 3329 and mod 257, forward and inverse (Thumb-2 + DSP + FPU registers).
 *
 * Data layout: int16 a[END_N], two coefficients per 32-bit word. The transform leaves are END_LEAF coefficients
 * (END_S = END_N/256 words), so every butterfly of every layer pairs whole words and the two packed halves always
 * share a twiddle.
 *
 * Forward (Cooley-Tukey, twiddles zeta[1..127] in the usual order):
 *   layers 1-4 : the 16 words of a column (8*END_S words apart) are transformed together (ct4_col in ntt_m4.i),
 *                the 15 twiddles stay in s8-s22;
 *   layers 5-7 : 8 words (END_S words apart) of a block, 7 twiddles per block read from memory.
 * Inverse (Cooley-Tukey butterflies on the bit-reversed input, then a twist; the method of the pqm4 ML-KEM code,
 *   see ntt_m4.i for the formula):
 *   layers m = 1, 2, 4, 8 : the 16 positions of a block. Three layers on the upper eight words, which are then
 *                multiplied by the twiddles e^i of the layer m = 8 and parked; three layers on the lower eight
 *                words; the layer m = 8 is one add/sub per pair. Seven of every twelve butterflies of the layers
 *                m = 1, 2, 4 have the twiddle 1 and need no multiplication.
 *   layers m = 16, 32, 64 : the 8 words of a column (16 positions apart), 7 twiddles per column, followed by
 *                the twist: every word is multiplied by psi^-P / 128. The 15 constants of a column stay in
 *                s0-s14 for its END_S words.
 *
 * Arithmetic: Plantard multiplication (plant.i). Products are in [-(q-1)/2, (q-1)/2] =: [-h, h] for any int16 input.
 * Range discipline (replayed by scripts/ntt_bounds.py, exercised by the extreme-input tests of test_kernels):
 *   forward : every layer adds at most h, nothing is reduced:  |out| <= |in| + 7h.
 *             mod 3329: |in| <= 16384 gives |out| <= 28032.  mod 257: the output is made canonical at the end.
 *   inverse : |in| <= B = 4095. A butterfly with a multiplication adds at most h, one without doubles.
 *             Upper half of a block: at most 8B = 32760 after three layers, then every word is multiplied (the
 *             first one by 1, only to reduce it). Lower half: the butterfly (c0, c4) of the layer m = 4 multiplies
 *             by 1, which keeps the words at 4B + h; after the layer m = 8 at most 4B + 2h.
 *             The last three layers add 3h: at most 4B + 5h = 24700 (q = 3329). The twist reduces to [-h, h].
 *             (pqm4 does without the two multiplications by 1 per block; that needs |in| <= 1735 for q = 3329,
 *             and the decoder passes values up to 3328.)
 *
 * Constant time: straight-line code, loop bounds and addresses are public.
 */
#include "asm_defs.h"
#ifdef END_M4_ASM

    .syntax unified
    .thumb
#include "plant.i"
#include "ntt_m4.i"

#define STEP1  (32*END_S)       /* bytes between the 16 words of a column (layers 1-4)   */
#define DS     (4*END_S)        /* bytes between the 8 words of a group   (layers 5-7)   */
#define COLB   (32*END_S)       /* bytes of all columns: 8*END_S words                    */
#define NBYTES (2*END_N)

/* any pair with |x| <= 17476 -> [0, 257): add 68*257, fold the two bytes (2^8 = -1), add 257 where negative.
 * r10 = 0x44444444 (68*257 twice), r12 = 0x01010101 (257 twice) */
.macro canon257_lazy c
    uadd16  \c, \c, r10
    uxtb16  r14, \c
    uxtb16  \c, \c, ror #8
    ssub16  \c, r14, \c
    uadd16  r14, \c, r12
    sel     \c, r14, \c
.endm

/* ------------------------------------------------------------------------------------------------------------
 * void name(int16_t a[END_N], const int32_t tw[127]);     tw = end_plant_ntt<q>
 * ------------------------------------------------------------------------------------------------------------ */
.macro ntt7_fwd name, Q, QA, CANON
    .global \name
    .type   \name, %function
    .align  2
\name:
    push    {r4-r11, lr}
    vpush   {s16-s25}
    movw    r11, #\QA
    movw    r12, #0
    movt    r12, #\Q
    /* ---- layers 1-4 ---- */
    vldm    r1!, {s8-s22}
    add.w   r14, r0, #COLB
    vmov    s25, r14
1:
    ct4_col STEP1, 0
    add.w   r0, r0, #4
    vmov    r14, s25
    cmp     r0, r14
    bne.w   1b
    /* ---- layers 5-7 ---- */
    sub.w   r0, r0, #COLB
    add.w   r14, r0, #NBYTES
    vmov    s25, r14
2:
    .set    j, 0
    .rept   END_S
    load8   (4*j), DS
    ldr.w   r10, [r1, #0]
    ct      r2, r6
    ct      r3, r7
    ct      r4, r8
    ct      r5, r9
    ldr.w   r10, [r1, #4]
    ct      r2, r4
    ct      r3, r5
    ldr.w   r10, [r1, #8]
    ct      r6, r8
    ct      r7, r9
    ldr.w   r10, [r1, #12]
    ct      r2, r3
    ldr.w   r10, [r1, #16]
    ct      r4, r5
    ldr.w   r10, [r1, #20]
    ct      r6, r7
    ldr.w   r10, [r1, #24]
    ct      r8, r9
    .if \CANON
    mov.w   r10, #0x44444444
    mov.w   r12, #0x01010101
    canon257_lazy r2
    canon257_lazy r3
    canon257_lazy r4
    canon257_lazy r5
    canon257_lazy r6
    canon257_lazy r7
    canon257_lazy r8
    canon257_lazy r9
    movw    r12, #0
    movt    r12, #\Q
    .endif
    store8  (4*j), DS
    .set    j, j+1
    .endr
    add.w   r1, r1, #28
    add.w   r0, r0, #(8*DS)
    vmov    r14, s25
    cmp     r0, r14
    bne.w   2b
    vpop    {s16-s25}
    pop     {r4-r11, pc}
    .size   \name, .-\name
.endm

/* ------------------------------------------------------------------------------------------------------------
 * void name(int16_t a[END_N], const int32_t tw[248]);     tw = end_plant_intt<q>
 * CANON : 1 = output in [0, q) (257), 0 = output centred in [-h, h]
 * The loop bodies handle two of the END_S words of a position; END-1024 runs them twice (unrolling all four
 * words costs 3 KB of code per function and saves 1.4% of its time).
 * ------------------------------------------------------------------------------------------------------------ */
#if END_S > 2
#define INV_SAVE {s16-s18}
#else
#define INV_SAVE {s16-s17}
#endif
.macro ntt7_inv name, Q, QA, CANON
    .global \name
    .type   \name, %function
    .align  2
\name:
    push    {r4-r11, lr}
    vpush   INV_SAVE
    movw    r11, #\QA
    movw    r12, #0
    movt    r12, #\Q
    /* ---- layers m = 1, 2, 4, 8: blocks of 16 positions ---- */
    vldm    r1!, {s8-s15}               /* e^0 .. e^7 */
    vmov    s16, r1                     /* r1 is scratch in ict3_light */
    add.w   r14, r0, #NBYTES
    vmov    s17, r14
1:
    .if END_S > 2
    add.w   r14, r0, #(4*END_S)
    vmov    s18, r14
3:
    .endif
    .set    j, 0
    .rept   2
    load8   (4*j+8*DS), DS              /* upper half: three layers, times e^i, parked */
    ict3_light 0
    ct_park r2, s8, s0
    ct_park r3, s9, s1
    ct_park r4, s10, s2
    ct_park r5, s11, s3
    ct_park r6, s12, s4
    ct_park r7, s13, s5
    ct_park r8, s14, s6
    ct_park r9, s15, s7
    load8   (4*j), DS                   /* lower half: three layers, then the layer m = 8 */
    ict3_light 1
    ct_last r2, s0, (4*j), (8*DS)
    ct_last r3, s1, (4*j+DS), (8*DS)
    ct_last r4, s2, (4*j+2*DS), (8*DS)
    ct_last r5, s3, (4*j+3*DS), (8*DS)
    ct_last r6, s4, (4*j+4*DS), (8*DS)
    ct_last r7, s5, (4*j+5*DS), (8*DS)
    ct_last r8, s6, (4*j+6*DS), (8*DS)
    ct_last r9, s7, (4*j+7*DS), (8*DS)
    .set    j, j+1
    .endr
    .if END_S > 2
    add.w   r0, r0, #8
    vmov    r14, s18
    cmp     r0, r14
    bne.w   3b
    add.w   r0, r0, #(16*DS-4*END_S)
    .else
    add.w   r0, r0, #(16*DS)
    .endif
    vmov    r14, s17
    cmp     r0, r14
    bne.w   1b
    /* ---- layers m = 16, 32, 64 and the twist: columns of 8 words ---- */
    sub.w   r0, r0, #NBYTES
    vmov    r1, s16
    add.w   r14, r0, #(16*DS)
    vmov    s15, r14
2:
    vldm    r1!, {s0-s14}
    .if END_S > 2
    vmov    s16, r1                     /* r1 = end of the column while its words are processed */
    add.w   r1, r0, #DS
4:
    .endif
    .set    j, 0
    .rept   2
    load8   (4*j), (16*DS)
    ict3_twist \CANON
    store8  (4*j), (16*DS)
    .set    j, j+1
    .endr
    add.w   r0, r0, #8
    .if END_S > 2
    cmp     r0, r1
    bne.w   4b
    vmov    r1, s16
    .endif
    vmov    r14, s15
    cmp     r0, r14
    bne.w   2b
    vpop    INV_SAVE
    pop     {r4-r11, pc}
    .size   \name, .-\name
.endm

    .text
    ntt7_fwd END_SYM(m4_ntt3329), 3329, 26632, 0
    ntt7_inv END_SYM(m4_intt3329), 3329, 26632, 0
    ntt7_fwd END_SYM(m4_ntt257), 257, 16448, 1
    ntt7_inv END_SYM(m4_intt257), 257, 16448, 1

#endif /* END_M4_ASM */
