/*
 * END KEM -- Cortex-M4 kernel: the complete NTT mod 12289 (forward, inverse) and the pointwise product (Thumb-2 + DSP + FPU registers).
 *
 * The transform has log2(END_N) layers (9 or 10). With two coefficients per word, all layers but the last pair
 * whole words; the last layer pairs the two halves of a word.
 *
 * Forward (Cooley-Tukey, twiddles zeta[1..END_N-1] in the usual order), |in| <= 8191:
 *   layers 1-4 : columns of 16 words, 15 twiddles in s8-s22 (ct4_col)                     -> <= 8191 + 4h = 32767
 *   layers 5-8 : 16 blocks, 16 words 1 (END-512) or 2 (END-1024) words apart, 15 twiddles
 *                per block; every word is reduced on entry                               -> <= 5h = 30720
 *   last       : per pair of adjacent words [x0|x1], [y0|y1].
 *                END-1024 only: layer 9 on the two words (first word reduced)             -> <= 2h
 *                in-word layer: U = [x0|y0], V = [x1|y1], one twiddle per half
 *                (END-512: U reduced first)                                               -> <= 2h (512), 3h (1024)
 *   The output is NOT reduced: any representative with |x| <= 3h. Only poly_pointwise12289 consumes it.
 * Inverse (Gentleman-Sande), |in| <= 8191, output in [-h, h]:
 *   first      : in-word layer per pair of words (sums <= 2|in|); END-1024: then the layer on the two words, sum reduced
 *   next four  : per block (gs4_col); END-512 reduces the eight sums of the first of them (<= 4|in| = 32764)
 *   layers 4-1 : per column (gs4_col), scaling by 1/n merged into the last layer
 *   Inside gs4_col the sums of the second layer are reduced, so that no value exceeds 4h = 24576.
 * h = 6144; Plantard arithmetic with alpha = 1 (plant.i). The schedule is replayed by scripts/ntt_bounds.py.
 */
#include "asm_defs.h"
#ifdef END_M4_ASM

    .syntax unified
    .thumb
#include "plant.i"
#include "ntt_m4.i"

#define STEP1  (32*END_S)       /* bytes between the 16 words of a column (layers 1-4)          */
#define COLB   (32*END_S)       /* bytes of all columns                                          */
#define STEPB  (2*END_S)        /* bytes between the 16 words of the next four layers            */
#define NCOLB  (END_S/2)        /* columns per block of the next four layers                     */
#define BLKB   (32*END_S)       /* bytes of a block                                              */
#define NBYTES (2*END_N)
#define Q      12289
#define QA     24578            /* q * 2^1 */
#if END_S == 2
#define B_REDSUM 1
#else
#define B_REDSUM 0
#endif

/* ------------------------------------------------------------------------------------------------------------
 * void m4_ntt12289(int16_t a[END_N], const int32_t *tw);     tw = end_plant_ntt12289_<n>
 * ------------------------------------------------------------------------------------------------------------ */
    .text
    .global END_SYM(m4_ntt12289)
    .type   END_SYM(m4_ntt12289), %function
    .align  2
END_SYM(m4_ntt12289):
    push    {r4-r11, lr}
    vpush   {s16-s25}
    movw    r11, #QA
    movw    r12, #0
    movt    r12, #Q
    vldr    s24, [r1]                   /* plant(1) */
    add.w   r1, r1, #4
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
    sub.w   r0, r0, #COLB
    /* ---- layers 5-8 ---- */
    add.w   r14, r0, #NBYTES
    vmov    s25, r14
2:
    vldm    r1!, {s8-s22}
    .rept   NCOLB
    ct4_col STEPB, 1
    add.w   r0, r0, #4
    .endr
    add.w   r0, r0, #(BLKB-4*NCOLB)
    vmov    r14, s25
    cmp     r0, r14
    bne.w   2b
    sub.w   r0, r0, #NBYTES
    /* ---- last layer(s), per pair of adjacent words ---- */
    vmov    r7, s24
    add.w   r8, r0, #NBYTES
3:
    ldr.w   r2, [r0]
    ldr.w   r3, [r0, #4]
#if END_S == 4
    ldr.w   r10, [r1], #4               /* layer 9 */
    doubleplant r2, r7, r14, r12, r11
    ct      r2, r3
#endif
    ldr.w   r5, [r1], #4
    ldr.w   r6, [r1], #4
    pkhbt   r4, r2, r3, lsl #16         /* U = [x0|y0] */
    pkhtb   r3, r3, r2, asr #16         /* V = [x1|y1] */
#if END_S == 2
    doubleplant r4, r7, r14, r12, r11
#endif
    smulwb  r14, r5, r3
    smulwt  r3, r6, r3
    smlabt  r14, r14, r12, r11
    smlabt  r3, r3, r12, r11
    pkhtb   r14, r3, r14, asr #16
    ssub16  r3, r4, r14
    sadd16  r4, r4, r14
    pkhbt   r2, r4, r3, lsl #16
    pkhtb   r3, r3, r4, asr #16
    str.w   r2, [r0], #4
    str.w   r3, [r0], #4
    cmp     r0, r8
    bne.w   3b
    vpop    {s16-s25}
    pop     {r4-r11, pc}
    .size   END_SYM(m4_ntt12289), .-END_SYM(m4_ntt12289)

/* ------------------------------------------------------------------------------------------------------------
 * void m4_intt12289(int16_t a[END_N], const int32_t *tw);     tw = end_plant_intt12289_<n>
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(m4_intt12289)
    .type   END_SYM(m4_intt12289), %function
    .align  2
END_SYM(m4_intt12289):
    push    {r4-r11, lr}
    vpush   {s16-s26}
    movw    r11, #QA
    movw    r12, #0
    movt    r12, #Q
    vldr    s24, [r1]                   /* plant(1) */
    add.w   r1, r1, #4
    /* ---- first layer(s), per pair of adjacent words ---- */
    vmov    r7, s24
    add.w   r8, r0, #NBYTES
1:
    ldr.w   r2, [r0]
    ldr.w   r3, [r0, #4]
    ldr.w   r5, [r1], #4
    ldr.w   r6, [r1], #4
    pkhbt   r4, r2, r3, lsl #16         /* U = [x0|y0] */
    pkhtb   r3, r3, r2, asr #16         /* V = [x1|y1] */
    ssub16  r14, r4, r3
    sadd16  r4, r4, r3
    smulwb  r3, r5, r14
    smulwt  r14, r6, r14
    smlabt  r3, r3, r12, r11
    smlabt  r14, r14, r12, r11
    pkhtb   r3, r14, r3, asr #16
    pkhbt   r2, r4, r3, lsl #16
    pkhtb   r3, r3, r4, asr #16
#if END_S == 4
    ldr.w   r10, [r1], #4               /* layer 9 */
    gs      r2, r3
    doubleplant r2, r7, r14, r12, r11
#endif
    str.w   r2, [r0], #4
    str.w   r3, [r0], #4
    cmp     r0, r8
    bne.w   1b
    sub.w   r0, r0, #NBYTES
    /* ---- next four layers, per block ---- */
    add.w   r14, r0, #NBYTES
    vmov    s25, r14
2:
    vldm    r1!, {s8-s23}
    vmov    s26, r1
    .rept   NCOLB
    gs4_col STEPB, B_REDSUM
    add.w   r0, r0, #4
    .endr
    vmov    r1, s26
    add.w   r0, r0, #(BLKB-4*NCOLB)
    vmov    r14, s25
    cmp     r0, r14
    bne.w   2b
    sub.w   r0, r0, #NBYTES
    /* ---- layers 4-1 ---- */
    vldm    r1, {s8-s23}
    add.w   r14, r0, #COLB
    vmov    s25, r14
3:
    gs4_col STEP1, 0
    add.w   r0, r0, #4
    vmov    r14, s25
    cmp     r0, r14
    bne.w   3b
    vpop    {s16-s26}
    pop     {r4-r11, pc}
    .size   END_SYM(m4_intt12289), .-END_SYM(m4_intt12289)

/* ------------------------------------------------------------------------------------------------------------
 * void m4_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N], int32_t plantc);
 *   r = a o b (true product). b is first multiplied by c = -2^32 mod q (plantc = plant(c)), which also reduces it
 *   to [-h, h]; the Plantard reduction of a_i * (b_i c) then returns a_i b_i in [-h, h]. a, b: any int16.
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(m4_pointwise12289)
    .type   END_SYM(m4_pointwise12289), %function
    .align  2
END_SYM(m4_pointwise12289):
    push    {r4-r11, lr}
    movw    r11, #QA
    movw    r12, #0
    movt    r12, #Q
    movw    r14, #0xd001                /* q^-1 mod 2^32 = 0x08ffd001 */
    movt    r14, #0x08ff
    add.w   r10, r0, #NBYTES
1:
    ldr.w   r4, [r1], #4
    ldr.w   r5, [r1], #4
    ldr.w   r6, [r2], #4
    ldr.w   r7, [r2], #4
    doubleplant r6, r3, r8, r12, r11
    doubleplant r7, r3, r8, r12, r11
    smulbb  r8, r4, r6
    smultt  r9, r4, r6
    plant_red r8, r14, r12, r11
    plant_red r9, r14, r12, r11
    pkhtb   r8, r9, r8, asr #16
    smulbb  r9, r5, r7
    smultt  r4, r5, r7
    plant_red r9, r14, r12, r11
    plant_red r4, r14, r12, r11
    pkhtb   r9, r4, r9, asr #16
    str.w   r8, [r0], #4
    str.w   r9, [r0], #4
    cmp     r0, r10
    bne.w   1b
    pop     {r4-r11, pc}
    .size   END_SYM(m4_pointwise12289), .-END_SYM(m4_pointwise12289)

#endif /* END_M4_ASM */
