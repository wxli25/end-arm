#!/usr/bin/env python3
"""Worst-case value bounds of the lazily reduced Cortex-M4 assembly kernels (src/armv7em/*.S).

The assembly keeps coefficients in int16 halves and reduces as rarely as possible. This script replays
the butterfly schedule of every kernel on *bounds* (largest possible absolute value per word) and
fails if a value could leave int16, or if a 32-bit accumulator could leave the valid interval of the
Plantard reduction. It must be kept in step with the reduction schedule of the assembly.

    python3 scripts/ntt_bounds.py          # prints one line per kernel, exit status 1 on a violation

Arithmetic facts used (verified exhaustively by scripts/gen_tables.py --check):
  * a Plantard twiddle product is in [-h, h], h = (q-1)/2, for any int16 input;
  * the 32-bit Plantard reduction is valid for c in [q*2^16 - q*2^(16+alpha), 2^31).
"""
import sys

INT16 = 32767
ALPHA = {257: 6, 3329: 3, 12289: 1}
ok = True


def report(name, worst, limit, note=""):
    global ok
    good = worst <= limit
    ok &= good
    print("%-34s max %11d / %11d  (%5.1f%%)  %s %s" % (name, worst, limit, 100.0 * worst / limit, "ok" if good else "OVERFLOW", note))


class Track:
    """Bounds of the words of one polynomial; remembers the largest value ever produced."""

    def __init__(self, nwords, bound):
        self.b = [bound] * nwords
        self.worst = bound

    def see(self, v):
        self.worst = max(self.worst, v)
        return v

    def ct(self, i0, i1, h):          # (a0, a1) <- (a0 + z a1, a0 - z a1)
        self.b[i0] = self.b[i1] = self.see(self.b[i0] + h)

    def light(self, i0, i1):          # (a0, a1) <- (a0 + a1, a0 - a1), no multiplication
        self.b[i0] = self.b[i1] = self.see(self.b[i0] + self.b[i1])

    def gs(self, i0, i1, h):          # (a0, a1) <- (a0 + a1, (a0 - a1) z): sum and difference have the same bound
        self.b[i0] = self.see(self.b[i0] + self.b[i1])
        self.b[i1] = h

    def red(self, i, h):              # Plantard product by a constant
        self.b[i] = h


def fwd7(q, n, bound_in):
    """ntt7_fwd: layers 1-4 per column of 16 words, layers 5-7 per group of 8 words; no reduction."""
    h, s = (q - 1) // 2, n // 256
    t = Track(n // 2, bound_in)
    step = 8 * s
    for col in range(8 * s):
        w = [col + m * step for m in range(16)]
        for parity in (1, 0):                      # odd words first (they are parked), then the even ones
            c = [w[2 * k + parity] for k in range(8)]
            for (i, j) in ((0, 4), (1, 5), (2, 6), (3, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 1), (2, 3), (4, 5), (6, 7)):
                t.ct(c[i], c[j], h)
            if parity:
                for k in range(8):
                    t.red(c[k], h)                 # multiplied by the layer-4 twiddle
        for k in range(8):                         # layer 4: even +- parked odd
            t.b[w[2 * k]] = t.b[w[2 * k + 1]] = t.see(t.b[w[2 * k]] + h)
    for blk in range(16):
        for j in range(s):
            c = [blk * 8 * s + j + k * s for k in range(8)]
            for (i, k) in ((0, 4), (1, 5), (2, 6), (3, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 1), (2, 3), (4, 5), (6, 7)):
                t.ct(c[i], c[k], h)
    return t


def inv7(q, n, bound_in):
    """ntt7_inv: Cooley-Tukey layers m = 1, 2, 4, 8 per block of 16 positions (ict3_light on the upper half, which
    is then multiplied and parked, and on the lower half with the (c0, c4) butterfly multiplying by 1), then the
    layers m = 16, 32, 64 and the twist per column (ict3_twist)."""
    h, s = (q - 1) // 2, n // 256
    t = Track(n // 2, bound_in)
    for blk in range(8):
        for j in range(s):
            for upper in (1, 0):
                c = [(blk * 16 + upper * 8 + k) * s + j for k in range(8)]
                for (i, k) in ((0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (4, 6)):
                    t.light(c[i], c[k])
                t.ct(c[1], c[3], h)
                t.ct(c[5], c[7], h)
                if upper:
                    t.light(c[0], c[4])
                else:
                    t.ct(c[0], c[4], h)
                t.ct(c[1], c[5], h)
                t.ct(c[2], c[6], h)
                t.ct(c[3], c[7], h)
                if upper:
                    for k in range(8):
                        t.red(c[k], h)             # multiplied by e^k
            for k in range(8):                     # layer m = 8: lower +- parked upper
                t.light((blk * 16 + k) * s + j, (blk * 16 + 8 + k) * s + j)
    for col in range(16):
        for j in range(s):
            c = [(col + 16 * k) * s + j for k in range(8)]
            for (i, k) in ((0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (4, 6), (1, 3), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7)):
                t.ct(c[i], c[k], h)
            for k in range(8):
                t.red(c[k], h)                     # twist
    return t


CT3 = ((0, 4), (1, 5), (2, 6), (3, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 1), (2, 3), (4, 5), (6, 7))
GS3A = ((0, 1), (2, 3), (4, 5), (6, 7))
GS3B = ((0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7))


def ct4_col(t, w, h, red):
    """ct4_col of ntt_m4.i on the 16 word indices w."""
    for parity in (1, 0):
        c = [w[2 * k + parity] for k in range(8)]
        if red:
            for k in range(8):
                t.red(c[k], h)
        for (i, j) in CT3:
            t.ct(c[i], c[j], h)
        if parity:
            for k in range(8):
                t.red(c[k], h)
    for k in range(8):
        t.b[w[2 * k]] = t.b[w[2 * k + 1]] = t.see(t.b[w[2 * k]] + h)


def gs4_col(t, w, h, redsum):
    """gs4_col of ntt_m4.i on the 16 word indices w."""
    for k in range(8):
        t.gs(w[2 * k], w[2 * k + 1], h)
        if redsum:
            t.red(w[2 * k], h)
    for parity in (0, 1):
        c = [w[2 * k + parity] for k in range(8)]
        for (i, j) in GS3A:
            t.gs(c[i], c[j], h)
        for k in (0, 2, 4, 6):
            t.red(c[k], h)
        for (i, j) in GS3B:
            t.gs(c[i], c[j], h)
        for k in range(4):
            t.red(c[k], h)


def fwd12289(n, bound_in):
    """m4_ntt12289: layers 1-4 per column, four more layers per block with entry reduction, last layer(s) per word pair."""
    q = 12289
    h, s = (q - 1) // 2, n // 256
    nw = n // 2
    t = Track(nw, bound_in)
    for col in range(8 * s):
        ct4_col(t, [col + m * 8 * s for m in range(16)], h, False)
    for blk in range(16):
        for col in range(s // 2):
            ct4_col(t, [blk * 8 * s + col + m * (s // 2) for m in range(16)], h, True)
    for i in range(0, nw, 2):
        if n == 1024:                              # layer 9 on two words, first one reduced
            t.red(i, h)
            t.ct(i, i + 1, h)
            t.b[i] = t.b[i + 1] = t.see(max(t.b[i], t.b[i + 1]) + h)     # in-word layer: U +- product
        else:                                      # in-word layer: U reduced, then U +- product
            t.b[i] = t.b[i + 1] = t.see(h + h)
    return t


def inv12289(n, bound_in):
    """m4_intt12289: first layer(s) per word pair, four layers per block, layers 4-1 per column."""
    q = 12289
    h, s = (q - 1) // 2, n // 256
    nw = n // 2
    t = Track(nw, bound_in)
    for i in range(0, nw, 2):
        t.b[i] = t.b[i + 1] = t.see(2 * max(t.b[i], t.b[i + 1]))        # in-word layer: sums of two halves
        if n == 1024:
            t.gs(i, i + 1, h)
            t.red(i, h)
    for blk in range(16):
        for col in range(s // 2):
            gs4_col(t, [blk * 8 * s + col + m * (s // 2) for m in range(16)], h, n == 512)
    for col in range(8 * s):
        gs4_col(t, [col + m * 8 * s for m in range(16)], h, False)
    return t


def plant_interval(q):
    return q * (1 << 16) * ((1 << ALPHA[q]) - 1)   # magnitude of the (smaller) negative end


for n in (512, 1024):
    print("END-%d" % n)
    # forward transforms: contract |in| <= 16384
    t = fwd7(3329, n, 16384)
    report("  ntt3329  (|in| <= 16384)", t.worst, INT16)
    t = fwd7(257, n, 16384)
    report("  ntt257   (|in| <= 16384)", t.worst, 17476, "limit of the final fold: 68*257")
    # inverse transforms: contract |in| <= 4095
    t = inv7(3329, n, 4095)
    report("  intt3329 (|in| <= 4095)", t.worst, INT16)
    t = inv7(257, n, 4095)
    report("  intt257  (|in| <= 4095)", t.worst, INT16)
    # leaf products: L dual MACs of an arbitrary int16 with a reduced operand
    leaf = n // 128
    for q in (3329, 257):
        report("  basemul%d accumulator" % q, leaf * 32768 * ((q - 1) // 2), plant_interval(q), "any int16 inputs")
    # complete transform mod 12289: contract |in| <= 8191 in both directions
    t = fwd12289(n, 8191)
    report("  ntt12289  (|in| <= 8191)", t.worst, INT16, "output bound %d" % max(t.b))
    t = inv12289(n, 8191)
    report("  intt12289 (|in| <= 8191)", t.worst, INT16)
    report("  pointwise12289 product", 32768 * 6144, plant_interval(12289), "any int16 inputs")

sys.exit(0 if ok else 1)
