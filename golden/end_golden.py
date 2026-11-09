#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
END KEM golden model (exact integer arithmetic).

This file is the executable specification used to validate the C backends:
  * polynomial arithmetic in Z[x]/(x^n+1) (exact, no modular shortcuts)
  * the reference NTT-257 layout used by the public-key wire format
  * pk / ct codecs (bit-exact with the reference implementation, commit dde0a77)
  * ACWC message layout, compression, SHA3 hashes
  * the trapdoor decoder in the "doubled / centred" domain (V1-V5 of the design brief):
        (mu'', -s'')  = (mu', -s') - 2*e~*(G,F) - 2*round((C2 - 2*q*w*e~)/(2qQ))*(g,f)
    with Free Candidate Localization (FCL) and distributional verification.

Everything here is deliberately slow and simple; it is the oracle, not the product.
"""
import sys, os, hashlib

try:
    import numpy as np
    HAVE_NP = True
except Exception:  # pragma: no cover
    HAVE_NP = False

Q = 257            # small decryption modulus
QP = 12289         # precision modulus Q'
QQP = Q * QP       # 3158273
TWO_QQP = 2 * QQP  # 6316546

# ---------------------------------------------------------------------------
# parameter sets
# ---------------------------------------------------------------------------
class Params(object):
    def __init__(self, n):
        self.n = n
        self.logn = n.bit_length() - 1
        if n == 512:
            self.de = 64; self.delta = 1; self.messbytes = 32
            self.pklen = 514; self.ctlen = 384; self.leaf = 4
            self.hash1 = 'sha3_256'; self.gamma = 6
            self.B_gf = (2, 1); self.B_GF = (2, 1)
        elif n == 1024:
            self.de = 90; self.delta = 0; self.messbytes = 64
            self.pklen = 1027; self.ctlen = 832; self.leaf = 8
            self.hash1 = 'sha3_512'; self.gamma = 3
            self.B_gf = (2, 1); self.B_GF = (2, 1)   # reference uses 2-fold family in both stages
        else:
            raise ValueError(n)
        self.m0_gf = 2
        self.m0_GF = 2

P512 = Params(512)
P1024 = Params(1024)

# ---------------------------------------------------------------------------
# exact polynomial arithmetic in Z[x]/(x^n+1)
# ---------------------------------------------------------------------------
def pmul(a, b, n):
    """exact negacyclic product (int64 via numpy when the ranges allow, else python ints)"""
    if HAVE_NP:
        ma = max(1, max(abs(int(x)) for x in a)); mb = max(1, max(abs(int(x)) for x in b))
        if ma * mb * n < (1 << 62):
            A = np.array(a, dtype=np.int64); B = np.array(b, dtype=np.int64)
            C = np.convolve(A, B)
            r = C[:n].copy()
            r[:n - 1] -= C[n:]
            return [int(x) for x in r]
    r = [0] * n
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            k = i + j
            if k < n:
                r[k] += ai * bj
            else:
                r[k - n] -= ai * bj
    return r

def pmul_mod(a, b, n, m):
    return [x % m for x in pmul([x % m for x in a], [x % m for x in b], n)]

def xpow_mul(a, k, n):
    """x^k * a  (negacyclic), 0 <= k < n"""
    r = [0] * n
    for i in range(n):
        j = i + k
        if j < n:
            r[j] += a[i]
        else:
            r[j - n] -= a[i]
    return r

def center(x, m):
    """centred representative in [-m/2, m/2)"""
    x %= m
    if x >= (m + 1) // 2:   # for odd m: x >= (m+1)/2 -> negative
        x -= m
    return x

def center_even(x, m):
    """centred representative for even m, in [-m/2, m/2)"""
    x %= m
    if x >= m // 2:
        x -= m
    return x

def round_div(x, d):
    """round(x/d) with ties rounded towards +inf  ( floor((x + d/2)/d) )"""
    return (x + d // 2) // d

# ---------------------------------------------------------------------------
# NTT mod 257 in the reference layout (7 layers, leaves of 4 (n=512) or 8 (n=1024))
# twiddle k is 3^brv7(k)  (3 is a primitive 256-th root of unity mod 257)
# ---------------------------------------------------------------------------
def brv(x, bits):
    r = 0
    for _ in range(bits):
        r = (r << 1) | (x & 1); x >>= 1
    return r

ZETAS257 = [pow(3, brv(k, 7), Q) for k in range(128)]
assert ZETAS257[:8] == [1, 241, 64, 4, 249, 128, 2, 225]

def ntt257(a, P):
    """forward NTT, output coefficients in [0,257), reference layout"""
    n, leaf = P.n, P.leaf
    a = [x % Q for x in a]
    k = 1
    length = n >> 1
    while length >= leaf:
        start = 0
        while start < n:
            z = ZETAS257[k]; k += 1
            for j in range(start, start + length):
                t = (z * a[j + length]) % Q
                a[j + length] = (a[j] - t) % Q
                a[j] = (a[j] + t) % Q
            start += 2 * length
        length >>= 1
    return a

def intt257(a, P):
    n, leaf = P.n, P.leaf
    a = [x % Q for x in a]
    length = leaf
    while length <= n >> 1:
        nblocks = n // (2 * length)
        start = 0
        b = 0
        while start < n:
            z = ZETAS257[nblocks + b]; zi = pow(z, Q - 2, Q)
            for j in range(start, start + length):
                t = a[j]
                a[j] = (t + a[j + length]) % Q
                a[j + length] = ((t - a[j + length]) * zi) % Q
            start += 2 * length; b += 1
        length <<= 1
    inv = pow(n // leaf, Q - 2, Q)
    return [(x * inv) % Q for x in a]

def leaf_zeta(P, blk):
    """modulus of leaf block blk is x^leaf - leaf_zeta"""
    z = ZETAS257[64 + (blk >> 1)]
    return z if (blk & 1) == 0 else (Q - z) % Q

def ntt257_mul(a, b, P):
    """pointwise leaf multiplication in the NTT domain"""
    n, L = P.n, P.leaf
    r = [0] * n
    for blk in range(n // L):
        z = leaf_zeta(P, blk); o = blk * L
        for i in range(L):
            for j in range(L):
                k = i + j
                if k < L:
                    r[o + k] += a[o + i] * b[o + j]
                else:
                    r[o + k - L] += z * a[o + i] * b[o + j]
        for k in range(L):
            r[o + k] %= Q
    return r

def _mat_inv_mod(M, p):
    L = len(M)
    A = [row[:] + [1 if i == j else 0 for j in range(L)] for i, row in enumerate(M)]
    for c in range(L):
        piv = None
        for r in range(c, L):
            if A[r][c] % p:
                piv = r; break
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p)
        A[c] = [(x * inv) % p for x in A[c]]
        for r in range(L):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [row[L:] for row in A]

def ntt257_inv(a, P):
    """inverse in the NTT domain; returns None if not invertible"""
    n, L = P.n, P.leaf
    r = [0] * n
    for blk in range(n // L):
        z = leaf_zeta(P, blk); o = blk * L
        # matrix of multiplication by a_blk in Z_257[x]/(x^L - z): column j = a * x^j
        M = [[0] * L for _ in range(L)]
        for j in range(L):
            for i in range(L):
                k = i + j
                if k < L:
                    M[k][j] = (M[k][j] + a[o + i]) % Q
                else:
                    M[k - L][j] = (M[k - L][j] + z * a[o + i]) % Q
        Mi = _mat_inv_mod(M, Q)
        if Mi is None:
            return None
        for i in range(L):
            r[o + i] = Mi[i][0] % Q    # inverse = M^-1 * e_0
    return r

# ---------------------------------------------------------------------------
# compression / decompression and message layout
# ---------------------------------------------------------------------------
def compress(x, P):
    if P.n == 512:
        return (x * 255) >> 10
    return (x * 717) >> 11

def decompress(c, P):
    if P.n == 512:
        return 4 * c + 3
    return ((c * 731) >> 8) + 2

def s_from_mr(mprime, r, P):
    """s[i*MB + j] = bit i of mprime[j]; s[n/2 + i*MB + j] = bit i of r[j]"""
    n, MB = P.n, P.messbytes
    s = [0] * n
    for i in range(8):
        for j in range(MB):
            s[i * MB + j] = (mprime[j] >> i) & 1
            s[n // 2 + i * MB + j] = (r[j] >> i) & 1
    return s

def mr_from_s(s, P):
    n, MB = P.n, P.messbytes
    mp = bytearray(MB); r = bytearray(MB)
    for i in range(8):
        for j in range(MB):
            mp[j] |= (s[i * MB + j] & 1) << i
            r[j] |= (s[n // 2 + i * MB + j] & 1) << i
    return bytes(mp), bytes(r)

def H1(data, P):
    return hashlib.new(P.hash1, data).digest()

def xor(a, b):
    return bytes(x ^ y for x, y in zip(a, b))

# ---------------------------------------------------------------------------
# ciphertext codec
# ---------------------------------------------------------------------------
C_TABLE90 = [90 ** i for i in range(8)]

def encode_ct(c, P):
    n = P.n
    if n == 512:
        out = bytearray(384)
        j = 0
        for i in range(0, n, 4):
            out[j] = ((c[i + 3] & 3) << 6) | c[i]
            out[j + 1] = (((c[i + 3] >> 2) & 3) << 6) | c[i + 1]
            out[j + 2] = (((c[i + 3] >> 4) & 3) << 6) | c[i + 2]
            j += 3
        return bytes(out)
    tmp = []
    for i in range(0, n, 8):
        tmp.append(sum(c[i + k] * C_TABLE90[k] for k in range(8)))
    res = [0] * 104
    index = 0
    for i in range(0, 80, 5):
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFF)
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFF)
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFF)
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFF)
        index += 4
        res[64] = ((res[64] << 4) | (tmp[i + 4] >> 48)) & (2 ** 64 - 1)
    index = 65
    for i in range(80, 125, 5):
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFF)
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFF)
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFF)
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFF)
        index += 4
        res[101] = ((res[101] << 4) | (tmp[i + 4] >> 48)) & (2 ** 64 - 1)
    res[101] = ((res[101] << 28) | (tmp[125] & 0x0FFFFFFF)) & (2 ** 64 - 1)
    res[102] = (tmp[126] << 12) | ((tmp[125] >> 28) & 0xFFF)
    res[103] = (tmp[127] << 12) | ((tmp[125] >> 40) & 0xFFF)
    out = bytearray()
    for x in res:
        out += (x & (2 ** 64 - 1)).to_bytes(8, 'little')
    return bytes(out)

def decode_ct(b, P):
    n = P.n
    if n == 512:
        c = [0] * n; j = 0
        for i in range(0, 384, 3):
            c[j] = b[i] & 0x3f; c[j + 1] = b[i + 1] & 0x3f; c[j + 2] = b[i + 2] & 0x3f
            c[j + 3] = (((b[i + 2] >> 6) & 3) << 4) | (((b[i + 1] >> 6) & 3) << 2) | ((b[i] >> 6) & 3)
            j += 4
        return c
    res = [int.from_bytes(b[8 * i:8 * i + 8], 'little') for i in range(104)]
    tmp = [0] * 128
    tmp[125] = ((res[103] & 0xFFF) << 40) | ((res[102] & 0xFFF) << 28) | (res[101] & 0x0FFFFFFF)
    tmp[126] = (res[102] >> 12) & 0xFFFFFFFFFFFFF
    tmp[127] = (res[103] >> 12) & 0xFFFFFFFFFFFFF
    index = 65; t = res[101] >> 28; j = 0
    for i in range(80, 125, 5):
        tmp[i + 0] = res[index + 0] >> 12
        tmp[i + 1] = res[index + 1] >> 12
        tmp[i + 2] = res[index + 2] >> 12
        tmp[i + 3] = res[index + 3] >> 12
        tmp[i + 4] = (((t >> (32 - j * 4)) & 15) << 48) | ((res[index + 3] & 0xFFF) << 36) | \
                     ((res[index + 2] & 0xFFF) << 24) | ((res[index + 1] & 0xFFF) << 12) | (res[index + 0] & 0xFFF)
        index += 4; j += 1
    index = 0; t = res[64]; j = 0
    for i in range(0, 80, 5):
        tmp[i + 0] = res[index + 0] >> 12
        tmp[i + 1] = res[index + 1] >> 12
        tmp[i + 2] = res[index + 2] >> 12
        tmp[i + 3] = res[index + 3] >> 12
        tmp[i + 4] = (((t >> (60 - j * 4)) & 15) << 48) | ((res[index + 3] & 0xFFF) << 36) | \
                     ((res[index + 2] & 0xFFF) << 24) | ((res[index + 1] & 0xFFF) << 12) | (res[index + 0] & 0xFFF)
        index += 4; j += 1
    c = [0] * n
    for blk in range(128):
        x = tmp[blk]
        for k in range(8):
            c[8 * blk + k] = x % 90; x //= 90
    return c

# ---------------------------------------------------------------------------
# public-key codec (mixed radix, values taken modulo 258 = 86*3)
# ---------------------------------------------------------------------------
PK_T86 = [86 ** i for i in range(7)]
PK_T3 = [3 ** i for i in range(17)]
M64 = 2 ** 64 - 1

def encode_pk(h, P):
    """h: NTT-257 domain coefficients (any int16 representation; reduced mod 258 like the reference)"""
    n = P.n
    if n == 512:
        nw, nrem, hb_off, w72_off, last_off, m3_off, nw3, m3rem = 73, 1, 384, 405, 411, 412, 30, 2
    else:
        nw, nrem, hb_off, w72_off, last_off, m3_off, nw3, m3rem = 147, 2, 768, 810, None, 823, 60, 4
    vals = [((x + 33024) % 258) for x in h]
    r86 = [v % 86 for v in vals]; r3 = [v % 3 for v in vals]
    out = bytearray(P.pklen)
    nfull = (n - nrem) // 7
    tmp = [0] * nw
    for w in range(nfull):
        tmp[w] = sum(r86[7 * w + j] * PK_T86[j] for j in range(7))
    if n == 1024:
        tmp[146] = r86[n - 2] * PK_T86[0] + r86[n - 1] * PK_T86[1]
    j = 0
    nlow = nfull if n == 512 else 144   # words stored in the 64-bit low section
    if n == 512:
        nlow = 72
    for i in range(0, nlow, 3):
        r0 = ((tmp[i + 1] << 19) | (tmp[i] & 0x7FFFF)) & M64
        r1 = ((tmp[i + 2] << 19) | ((tmp[i] >> 19) & 0x7FFFF)) & M64
        out[8 * j:8 * j + 8] = r0.to_bytes(8, 'little')
        out[8 * (j + 1):8 * (j + 1) + 8] = r1.to_bytes(8, 'little')
        j += 2
    idx = hb_off
    for i in range(0, nlow, 24):
        for k in range(7):
            out[idx + k] = (((tmp[i + 3 * k] >> 38) << 1) | ((tmp[i + 21] >> (38 + k)) & 1)) & 0xFF
        idx += 7
    if n == 512:
        for i in range(6):
            out[idx] = (tmp[72] >> (8 * i)) & 0xFF; idx += 1
        out[idx] = r86[n - 1]; idx += 1
    else:
        for i in range(5):
            out[idx] = (tmp[144] >> (8 * i)) & 0xFF; idx += 1
        for i in range(5):
            out[idx] = (tmp[145] >> (8 * i)) & 0xFF; idx += 1
        out[idx] = ((((tmp[144] >> 40) & 31) << 3) | (tmp[146] & 7)) & 0xFF; idx += 1
        out[idx] = ((((tmp[145] >> 40) & 31) << 3) | ((tmp[146] >> 3) & 7)) & 0xFF; idx += 1
        out[idx] = (tmp[146] >> 6) & 0xFF; idx += 1
    # base-3 part
    tmp32 = [0] * (nw3 + 1)
    for w in range(nw3):
        tmp32[w] = sum(r3[17 * w + j] * PK_T3[j] for j in range(17))
    tmp32[nw3] = sum(r3[n - m3rem + i] * PK_T3[i] for i in range(m3rem))
    j = 0
    t = tmp32[nw3] if n == 1024 else 0
    for i in range(0, nw3, 6):
        for k in range(5):
            rk = (tmp32[i + k] | (((tmp32[i + 5] >> (5 * k)) & 31) << 27)) & 0xFFFFFFFF
            out[m3_off + 4 * (j + k):m3_off + 4 * (j + k) + 4] = rk.to_bytes(4, 'little')
        j += 5
        t = ((t << 2) | (tmp32[i + 5] >> 25)) & 0xFFFFFFFF
    if n == 512:
        out[512] = t & 0xFF
        out[513] = ((t >> 8) | (tmp32[30] << 2)) & 0xFF
    else:
        out[1023] = t & 0xFF; out[1024] = (t >> 8) & 0xFF; out[1025] = (t >> 16) & 0xFF; out[1026] = (t >> 24) & 0xFF
    return bytes(out)

def decode_pk(b, P):
    n = P.n
    r86 = [0] * n; r3 = [0] * n
    if n == 512:
        nlow, nrem, hb_off = 72, 1, 384
    else:
        nlow, nrem, hb_off = 144, 2, 768
    tmp = [0] * (nlow + 3)
    idx = hb_off
    for i in range(0, nlow, 24):
        for k in range(7):
            tmp[i + 3 * k] = (b[idx + k] >> 1) << 38
        tmp[i + 21] = sum(((b[idx + k] & 1) << (38 + k)) for k in range(7))
        idx += 7
    j = 0
    for i in range(0, nlow, 3):
        r0 = int.from_bytes(b[8 * j:8 * j + 8], 'little'); r1 = int.from_bytes(b[8 * j + 8:8 * j + 16], 'little')
        tmp[i + 0] += ((r1 & 0x7FFFF) << 19) | (r0 & 0x7FFFF)
        tmp[i + 1] = r0 >> 19
        tmp[i + 2] = r1 >> 19
        j += 2
    if n == 512:
        tmp[72] = int.from_bytes(b[405:411], 'little')
        r86[n - 1] = b[411]
        nfull = 73
    else:
        tmp[144] = int.from_bytes(b[810:815], 'little') | (((b[820] >> 3) & 31) << 40)
        tmp[145] = int.from_bytes(b[815:820], 'little') | (((b[821] >> 3) & 31) << 40)
        tmp[146] = (b[820] & 7) | ((b[821] & 7) << 3) | (b[822] << 6)
        nfull = 146
    for w in range(nfull):
        x = tmp[w]
        for k in range(7):
            r86[7 * w + k] = x % 86; x //= 86
    if n == 1024:
        x = tmp[146]
        r86[n - 2] = x % 86; r86[n - 1] = (x // 86) % 86
    # base-3 part
    if n == 512:
        m3_off, nw3, m3rem = 412, 30, 2
        tmp32 = [0] * 31
        tmp32[30] = b[513] >> 2
        t = ((b[513] & 3) << 8) | b[512]
    else:
        m3_off, nw3, m3rem = 823, 60, 4
        tmp32 = [0] * 61
        tail = int.from_bytes(b[1023:1027], 'little')
        t = tail
    j = 0
    for i in range(0, nw3, 6):
        rs = [int.from_bytes(b[m3_off + 4 * (j + k):m3_off + 4 * (j + k) + 4], 'little') for k in range(5)]
        for k in range(5):
            tmp32[i + k] = rs[k] & 0x07FFFFFF
        tmp32[i + 5] = sum(((rs[k] >> 27) << (5 * k)) for k in range(5))
        j += 5
    ngroups = nw3 // 6
    # the 2-bit top fragments were shifted in group order: last group in the lowest bits
    for gi in range(ngroups):
        i = 6 * gi
        tmp32[i + 5] += ((t >> (2 * (ngroups - 1 - gi))) & 3) << 25
    if n == 1024:
        tmp32[60] = t >> (2 * ngroups)
    for w in range(nw3):
        x = tmp32[w]
        for k in range(17):
            r3[17 * w + k] = x % 3; x //= 3
    x = tmp32[nw3]
    for k in range(m3rem):
        r3[n - m3rem + k] = x % 3; x //= 3
    return [(r86[i] * 87 + r3[i] * 172) % 258 for i in range(n)]

# ---------------------------------------------------------------------------
# encryption (DPKE) and the trapdoor decoder
# ---------------------------------------------------------------------------
def dpke_encrypt(h_ntt, s, P):
    """h_ntt: NTT-domain public key; s: binary list. returns (ct_coeffs in [0,de), x = h*s mod q)"""
    x = intt257(ntt257_mul(h_ntt, ntt257(s, P), P), P)   # h*s mod 257 in [0,257)
    return [compress(v, P) for v in x], x

def lift_center_parity(r, parity):
    """r in [0,257): representative mod 514 with given parity, centred to [-257, 256]"""
    r %= Q
    if (r & 1) != parity:
        r += Q
    if r > 256:
        r -= 2 * Q
    return r

def expanded_key(sk, P):
    f, g, F, G, w = sk
    n = P.n
    one = [1] * n
    ek = {
        'f': f, 'g': g, 'F': F, 'G': G, 'w': w,
        'g1': pmul(g, one, n), 'f1': pmul(f, one, n), 'G1': pmul(G, one, n), 'F1': pmul(F, one, n),
        'FQ': [QP * F[i] - x for i, x in enumerate(pmul(w, f, n))],
        'GQ': [QP * G[i] - x for i, x in enumerate(pmul(w, g, n))],
        'p1': (sum(g) + P.delta * sum(f)) & 1,
        'p2': (sum(G) + P.delta * sum(F)) & 1,
        'pw': (sum(w) & 1),
    }
    ek['pwc1'] = ek['p1'] & ek['pw']     # parity of every coefficient of w*c1
    return ek

def fcl(scores, m0):
    """indices of the m0 smallest scores (ties: smaller index first)"""
    return sorted(range(len(scores)), key=lambda i: (scores[i], i))[:m0]

def candidates(S, B0, B1):
    """candidate wrap polynomials as dicts {index: sigma}; zero candidate first"""
    from itertools import combinations, product
    out = [dict()]
    for k in range(1, B0 + 1):
        for I in combinations(S, k):
            for sig in product([v for v in range(-B1, B1 + 1) if v != 0], repeat=k):
                out.append(dict(zip(I, sig)))
    return out

def apply_sparse(poly, cand, n, scale=1):
    """scale * sum_i sigma_i x^i * poly"""
    r = [0] * n
    for i, sg in cand.items():
        t = xpow_mul(poly, i, n)
        for k in range(n):
            r[k] += scale * sg * t[k]
    return r

def bnp_initial(ek, cprime, P):
    """Steps 1-4 of Alg.6 in the doubled domain. Returns dict with c1, c2, C2, sbar_p, mubar_p."""
    n = P.n
    f, g, F, G, w = ek['f'], ek['g'], ek['F'], ek['G'], ek['w']
    fc = pmul(f, cprime, n); Fc = pmul(F, cprime, n)
    c1 = [lift_center_parity(2 * fc[i] - ek['g1'][i] - P.delta * ek['f1'][i], ek['p1']) for i in range(n)]
    c2 = [lift_center_parity(2 * Fc[i] - ek['G1'][i] - P.delta * ek['F1'][i], ek['p2']) for i in range(n)]
    u = pmul(w, c1, n)                                   # exact w*c1
    assert all((x & 1) == ek['pwc1'] for x in u), "parity structure of w*c1 violated"
    C2 = [center_even(QP * c2[i] - u[i], TWO_QQP) for i in range(n)]
    num_s = pmul(c1, ek['FQ'], n); t = pmul(C2, f, n)
    num_s = [num_s[i] - t[i] for i in range(n)]
    num_mu = pmul(C2, g, n); t = pmul(c1, ek['GQ'], n)
    num_mu = [num_mu[i] - t[i] for i in range(n)]
    assert all(x % QQP == 0 for x in num_s) and all(x % QQP == 0 for x in num_mu), "qQ does not divide (V1 violated)"
    return {'c1': c1, 'c2': c2, 'C2': C2, 'u': u,
            'sbar_p': [x // QQP for x in num_s], 'mubar_p': [x // QQP for x in num_mu]}

def decode_trap(ek, cprime, P, full_score=False, diag=None):
    """Alg.6 (END.DecodeTrap) with the corrected signs; returns (s, sbar, mubar) or None on failure."""
    n = P.n
    f, g, F, G, w = ek['f'], ek['g'], ek['F'], ek['G'], ek['w']
    st = bnp_initial(ek, cprime, P)
    c1, C2, sbar_p, mubar_p = st['c1'], st['C2'], st['sbar_p'], st['mubar_p']
    # ---- gf stage --------------------------------------------------------
    d1 = [Q - abs(x) for x in c1]
    S1 = fcl(d1, P.m0_gf)
    best = None
    for cand in candidates(S1, P.B_gf[0], P.B_gf[1]):
        we = apply_sparse(w, cand, n)                     # w * e~
        r = [round_div(C2[i] - 2 * Q * we[i], TWO_QQP) for i in range(n)]
        rf = pmul(r, f, n)
        eF = apply_sparse(F, cand, n)
        s2 = [sbar_p[i] + 2 * eF[i] + 2 * rf[i] for i in range(n)]
        score = sum(x * x for x in s2)
        if full_score:
            rg = pmul(r, g, n); eG = apply_sparse(G, cand, n)
            mu2 = [mubar_p[i] - 2 * eG[i] - 2 * rg[i] for i in range(n)]
            score += sum(x * x for x in mu2)
        if best is None or score < best[0]:
            C2t = [C2[i] - 2 * Q * we[i] - TWO_QQP * r[i] for i in range(n)]
            best = (score, cand, s2, C2t)
    score, cand_gf, sbar_pp, C2t = best
    if diag is not None:
        diag['S_gf'] = S1; diag['cand_gf'] = cand_gf; diag['score_gf'] = score
        diag['c1'] = c1; diag['C2'] = C2; diag['C2t'] = C2t
    # ---- GF stage --------------------------------------------------------
    d2 = [QQP - abs(x) for x in C2t]
    S2 = fcl(d2, P.m0_GF)
    for cand in candidates(S2, P.B_GF[0], P.B_GF[1]):
        ef = apply_sparse(f, cand, n)
        sbar = [sbar_pp[i] - 2 * ef[i] for i in range(n)]
        if all(x == 1 or x == -1 for x in sbar):
            if diag is not None:
                diag['S_GF'] = S2; diag['cand_GF'] = cand
            s = [(x + 1) // 2 for x in sbar]
            return s, sbar
    if diag is not None:
        diag['S_GF'] = S2; diag['cand_GF'] = None
    return None

# ---------------------------------------------------------------------------
# PKE / KEM wrappers
# ---------------------------------------------------------------------------
def pke_encrypt(h_ntt, m, r, P):
    mprime = xor(m, H1(r, P))
    s = s_from_mr(mprime, r, P)
    c, _ = dpke_encrypt(h_ntt, s, P)
    return encode_ct(c, P), s

def pke_decrypt(ek, ct, P, diag=None):
    c = decode_ct(ct, P)
    cprime = [decompress(v, P) for v in c]
    res = decode_trap(ek, cprime, P, diag=diag)
    if res is None:
        return None
    s, _ = res
    mprime, r = mr_from_s(s, P)
    return xor(mprime, H1(r, P)), r

def kem_encaps(h_ntt, m, P):
    r = H1(m, P)                     # H2 = H1 in the reference
    ct, _ = pke_encrypt(h_ntt, m, r, P)
    K = hashlib.sha3_256(m + ct).digest()
    return ct, K

def kem_decaps(ek, h_ntt, k, ct, P, diag=None):
    res = pke_decrypt(ek, ct, P, diag=diag)
    ok = False
    if res is not None:
        m, _ = res
        ct2, _ = pke_encrypt(h_ntt, m, H1(m, P), P)
        ok = (ct2 == ct)
    if not ok:
        m = k
    return hashlib.sha3_256(m + ct).digest(), ok

# ---------------------------------------------------------------------------
# KAT utilities
# ---------------------------------------------------------------------------
def parse_kat(path):
    entries = []; cur = {}
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line:
                if cur:
                    entries.append(cur); cur = {}
                continue
            k, v = [x.strip() for x in line.split('=', 1)]
            cur[k] = v
    if cur:
        entries.append(cur)
    return entries

def i8(hexs):
    b = bytes.fromhex(hexs)
    return [x - 256 if x > 127 else x for x in b]

def i16(hexs):
    b = bytes.fromhex(hexs)
    return [int.from_bytes(b[2 * i:2 * i + 2], 'little', signed=True) for i in range(len(b) // 2)]

def sk_from_entry(e):
    return (i8(e['f']), i8(e['g']), i8(e['F']), i8(e['G']), i16(e['w']))

def check_entry(e, P, verbose=False):
    n = P.n
    f, g, F, G, w = sk_from_entry(e)
    pk = bytes.fromhex(e['pk']); ct = bytes.fromhex(e['ct']); m = bytes.fromhex(e['m'])
    ss = bytes.fromhex(e['ss']); k = bytes.fromhex(e['k'])
    # NTRU equation
    gF = pmul(g, F, n); fG = pmul(f, G, n)
    eq = [gF[i] - fG[i] for i in range(n)]
    assert eq[0] == Q and all(x == 0 for x in eq[1:]), "gF - fG != q"
    # public key
    h_ntt = decode_pk(pk, P)
    assert encode_pk(h_ntt, P) == pk, "pk codec round trip"
    fi = ntt257_inv(ntt257(f, P), P)
    assert fi is not None
    h2 = ntt257_mul(ntt257(g, P), fi, P)
    assert h2 == h_ntt, "h != g/f"
    # encaps
    ct2, K = kem_encaps(h_ntt, m, P)
    assert ct2 == ct, "ciphertext mismatch"
    assert K == ss, "shared secret mismatch"
    # decaps
    ek = expanded_key((f, g, F, G, w), P)
    diag = {}
    K2, ok = kem_decaps(ek, h_ntt, k, ct, P, diag=diag)
    assert ok and K2 == ss, "decaps mismatch"
    assert diag['cand_gf'] == {} and diag['cand_GF'] == {}, "unexpected wrap in honest KAT"
    return ek

def parity_stats(entries, P):
    """how often would the 'even representative' CRT of the reference be wrong?"""
    cnt = {'p1': 0, 'pw': 0, 'odd_wc1': 0}
    for e in entries:
        ek = expanded_key(sk_from_entry(e), P)
        cnt['p1'] += ek['p1']; cnt['pw'] += ek['pw']; cnt['odd_wc1'] += ek['pwc1']
    return cnt

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('kat'); ap.add_argument('--n', type=int, default=512); ap.add_argument('--max', type=int, default=20)
    a = ap.parse_args()
    P = P512 if a.n == 512 else P1024
    ents = parse_kat(a.kat)[:a.max]
    for i, e in enumerate(ents):
        check_entry(e, P)
    print("checked %d KAT entries for END-%d: all OK (numpy=%s)" % (len(ents), a.n, HAVE_NP))
    print("parity stats:", parity_stats(ents, P))
