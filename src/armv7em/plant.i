/*
 * END KEM -- Cortex-M4 assembly macros: 16-bit Plantard arithmetic on packed coefficient pairs.
 *
 * One 32-bit register holds two consecutive int16 coefficients (low half = even index).
 * Plantard arithmetic follows Huang et al., "Improved Plantard Arithmetic for Lattice-based
 * Cryptography" (ePrint 2022/956), the instruction sequences are those of the pqm4 ML-KEM code.
 *
 * Register conventions of every macro:
 *   q   : the modulus in the TOP half of the register (bottom half ignored)
 *   qa  : q * 2^alpha as a 32-bit value (alpha = 6 / 3 / 1 for q = 257 / 3329 / 12289)
 *   tw  : plant(z) = ((z * (-2^32) mod q) * q^-1) mod 2^32 for a twiddle z (scripts/gen_tables.py)
 *
 * Facts verified exhaustively by `scripts/gen_tables.py --check`:
 *   - a twiddle multiplication accepts ANY int16 input and returns a*z mod q in [-(q-1)/2, (q-1)/2];
 *   - plant_red maps a 32-bit c in [q*2^16 - q*2^(16+alpha), 2^31) to c * (-2^-32) mod q, same range.
 * All instructions used here have data-independent timing on Cortex-M4.
 */

/* a <- a * z mod q for both halves */
.macro doubleplant a, tw, tmp, q, qa
    smulwb  \tmp, \tw, \a
    smulwt  \a, \tw, \a
    smlabt  \tmp, \tmp, \q, \qa
    smlabt  \a, \a, \q, \qa
    pkhtb   \a, \a, \tmp, asr #16
.endm

/* Cooley-Tukey butterfly with the results in other registers: (d0, d1) <- (a0 + z*a1, a0 - z*a1).
 * Each output is at most (q-1)/2 larger than a0. a1 is destroyed; d1 must not be a0. */
.macro ct_bfly_to d0, d1, a0, a1, tw, tmp, q, qa
    smulwb  \tmp, \tw, \a1
    smulwt  \a1, \tw, \a1
    smlabt  \tmp, \tmp, \q, \qa
    smlabt  \a1, \a1, \q, \qa
    pkhtb   \tmp, \a1, \tmp, asr #16
    ssub16  \d1, \a0, \tmp
    sadd16  \d0, \a0, \tmp
.endm

/* Cooley-Tukey butterfly in place: (a0, a1) <- (a0 + z*a1, a0 - z*a1) */
.macro ct_bfly a0, a1, tw, tmp, q, qa
    ct_bfly_to \a0, \a1, \a0, \a1, \tw, \tmp, \q, \qa
.endm

/* Gentleman-Sande butterfly: (a0, a1) <- (a0 + a1, (a0 - a1)*z); the sum doubles, the product is reduced */
.macro gs_bfly a0, a1, tw, tmp, q, qa
    ssub16  \tmp, \a0, \a1
    sadd16  \a0, \a0, \a1
    smulwb  \a1, \tw, \tmp
    smulwt  \tmp, \tw, \tmp
    smlabt  \a1, \a1, \q, \qa
    smlabt  \tmp, \tmp, \q, \qa
    pkhtb   \a1, \tmp, \a1, asr #16
.endm

/* a <- a * (-2^-32) mod q, result in the TOP half of a (bottom half is garbage) */
.macro plant_red a, qinv, q, qa
    mul     \a, \a, \qinv
    smlatt  \a, \a, \q, \qa
.endm

/* both halves in (-q, q) -> [0, q): qq holds q in both halves.
 * UADD16 sets the GE flags of a half exactly when the 16-bit addition carries, i.e. when the half is negative. */
.macro canon_pair a, tmp, qq
    uadd16  \tmp, \a, \qq
    sel     \a, \tmp, \a
.endm
