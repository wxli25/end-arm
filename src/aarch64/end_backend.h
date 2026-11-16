/*
 * END KEM -- backend configuration of the AArch64 NEON backend.
 *
 * The functions below are implemented with NEON in polyvec.c, so the portable versions in
 * src/common/end_deckern.c are left out.
 */
#ifndef END_BACKEND_H
#define END_BACKEND_H

#define END_ASM_FCL_SELECT        1
#define END_ASM_POLY_LIFT257      1
#define END_ASM_POLY_CAND_ROUND   1
#define END_ASM_POLY_CAND_COMBINE 1
#define END_ASM_POLY_CAND_CHECK   1
#define END_ASM_POLY_CMOV32       1

#endif /* END_BACKEND_H */
