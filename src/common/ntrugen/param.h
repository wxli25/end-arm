/* shim for the vendored ntrugen-derived files: map END_N onto the reference names */
#ifndef NTRUGEN_PARAM_SHIM_H
#define NTRUGEN_PARAM_SHIM_H
#ifndef END_N
#define END_N 512
#endif
#define DIM_N END_N
#if END_N == 512
#define LOGN 9
#else
#define LOGN 10
#endif
#endif
