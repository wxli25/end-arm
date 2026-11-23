/*
 * END KEM -- backend configuration of the ARMv7E-M backend.
 *
 * deckern_m4.S provides CRT on integer DSP cores and the other overrides
 * when END_M4_ASM (BACKEND=m4f) enables the kernels using FPU registers.
 */
#ifndef END_BACKEND_H
#define END_BACKEND_H

#if defined(__ARM_FEATURE_DSP) && !defined(END_M_PORTABLE_CRT)
#define END_ASM_POLY_CRT_C2       1
#endif

#ifdef END_M4_ASM
#define END_ASM_POLY_LIFT257      1
#define END_ASM_POLY_CAND_ROUND   1
#define END_ASM_POLY_CAND_COMBINE 1
#define END_ASM_POLY_CAND_CHECK   1
#define END_ASM_POLY_CMOV32       1
#endif

#endif /* END_BACKEND_H */
