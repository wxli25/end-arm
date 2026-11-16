/*
 * END KEM -- definitions shared by the Cortex-M4 assembly kernels.
 * Macros only: this file is included from .S files, where end_params.h (C declarations) cannot be used.
 */
#ifndef END_ASM_DEFS_H
#define END_ASM_DEFS_H

#ifndef END_N
#define END_N 512
#endif

/* END_S = END_N / 256: number of 32-bit words (coefficient pairs) of one leaf of the 7-layer transforms */
#if END_N == 512
#define END_SYM(s) end512_##s
#define END_S 2
#elif END_N == 1024
#define END_SYM(s) end1024_##s
#define END_S 4
#else
#error "END_N must be 512 or 1024"
#endif

#endif /* END_ASM_DEFS_H */
