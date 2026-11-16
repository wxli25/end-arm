/*
 * END KEM -- backend configuration of the ARMv7E-M backend.
 *
 * With END_M4_ASM (BACKEND=m4f) the functions below are provided by the hand-written assembly in
 * deckern_m4.S, so the portable versions in src/common/end_deckern.c are left out.
 */
#ifndef END_BACKEND_H
#define END_BACKEND_H

#ifdef END_M4_ASM
#define END_ASM_POLY_LIFT257      1
#define END_ASM_POLY_CAND_ROUND   1
#define END_ASM_POLY_CAND_COMBINE 1
#define END_ASM_POLY_CAND_CHECK   1
#define END_ASM_POLY_CMOV32       1
#endif

#endif /* END_BACKEND_H */
