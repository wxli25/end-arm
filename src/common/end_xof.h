/*
 * END KEM -- SHAKE256 XOF adapter used by key generation.
 *
 * The library needs three operations on one SHAKE256 state:
 *   end_xof_absorb_once  : absorb the seed and finalize
 *   end_xof_squeeze      : streaming output (continues at the current position)
 *   end_xof_squeezeblocks: permute and output whole blocks without touching the
 *                          streaming position (pq-crystals squeezeblocks semantics,
 *                          used by the END-1024 ternary sampler)
 *
 * Default build: pq-crystals fips202.c shipped in src/common.
 * pqm4 build (END_PQM4): pqm4's mupq/common fips202 on top of its Keccak assembly,
 * so that pqm4's "hashing" benchmark accounts for the hash time. Both produce the
 * same bytes (checked by the END KAT in pqm4's test/testvectors runs).
 */
#ifndef END_XOF_H
#define END_XOF_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "fips202.h"

#ifdef END_PQM4
typedef shake256incctx end_xof_state;

static inline void end_xof_absorb_once(end_xof_state *st, const uint8_t *in, size_t inlen)
{
    shake256_inc_init(st);
    shake256_inc_absorb(st, in, inlen);
    shake256_inc_finalize(st);
}
static inline void end_xof_squeeze(uint8_t *out, size_t outlen, end_xof_state *st)
{
    shake256_inc_squeeze(out, outlen, st);
}
static inline void end_xof_squeezeblocks(uint8_t *out, size_t nblocks, end_xof_state *st)
{
    shake256ctx blocks;   /* same Keccak state words, streaming position (ctx[25]) left as is */
    memcpy(blocks.ctx, st->ctx, sizeof blocks.ctx);
    shake256_squeezeblocks(out, nblocks, &blocks);
    memcpy(st->ctx, blocks.ctx, sizeof blocks.ctx);
}
static inline void end_xof_release(end_xof_state *st)
{
    shake256_inc_ctx_release(st);
}
#else
typedef keccak_state end_xof_state;

static inline void end_xof_absorb_once(end_xof_state *st, const uint8_t *in, size_t inlen)
{
    shake256_absorb_once(st, in, inlen);
}
static inline void end_xof_squeeze(uint8_t *out, size_t outlen, end_xof_state *st)
{
    shake256_squeeze(out, outlen, st);
}
static inline void end_xof_squeezeblocks(uint8_t *out, size_t nblocks, end_xof_state *st)
{
    shake256_squeezeblocks(out, nblocks, st);
}
static inline void end_xof_release(end_xof_state *st)
{
    (void)st;
}
#endif

#endif /* END_XOF_H */
