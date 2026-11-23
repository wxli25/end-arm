/* Standalone Cortex-M DSP prototype for multiplication by the all-ones
 * polynomial mod 257. Input lanes are canonical residues. See bench_ones257.c. */
#include "../src/armv7em/asm_defs.h"

    .syntax unified
    .thumb
    .text
    .align 2
    .global END_SYM(bench_ones257_dsp)
    .type END_SYM(bench_ones257_dsp), %function
END_SYM(bench_ones257_dsp):
    push    {r4-r8, lr}
    mov     r2, r1
    add     r3, r1, #(2*END_N)
    mov     r4, #0x00010001
    movs    r5, #0
.Lones_total:
    .rept 4
    ldr     r6, [r2], #4
    smlad   r5, r6, r4, r5               /* add both coefficients */
    .endr
    cmp     r2, r3
    bne     .Lones_total

    /* Each prefix x = 257*N - total + 2*sum(a[0..i]) is positive and
     * bounded by 513*N. Multiplication by ceil(2^32/257) therefore
     * gives the exact quotient in its high word (valid for x < 2^24). */
    movw    r6, #((257*END_N) & 0xffff)
    movt    r6, #((257*END_N) >> 16)
    sub     r5, r6, r5
    lsl     r4, r4, #1                   /* two copies of 2 */
    movw    r7, #257
    movw    r8, #0x0100
    movt    r8, #0x00ff
.Lones_prefix:
    .rept 2
    ldr     r6, [r1], #4
    smlabb  r12, r6, r4, r5              /* prefix of first coefficient */
    smlad   r5, r6, r4, r5               /* carry raw prefix between pairs */
    smmul   r2, r12, r8
    smmul   lr, r5, r8
    mls     r12, r2, r7, r12
    mls     lr, lr, r7, r5
    pkhbt   r12, r12, lr, lsl #16
    str     r12, [r0], #4
    .endr
    cmp     r1, r3
    bne     .Lones_prefix
    pop     {r4-r8, pc}
    .size END_SYM(bench_ones257_dsp), .-END_SYM(bench_ones257_dsp)
