#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wrap-injection tests for the golden decoder (design brief V1-V5, S1, S6).

For a reference key (from the KAT files) we construct decompressed ciphertexts c' = h*s + mu' (mod q)
whose unreduced residues wrap in chosen coordinates, and check that
  * the V1/V5 identities hold exactly with the true wrap vectors e1, e2;
  * decode_trap recovers s (gf single/double wrap, GF single wrap, gf+GF);
  * the reference's gf-stage FCL score (128 - |c1|, ranked by |.|) does NOT localise the wrap,
    while the correct boundary score (257 - |c1|) does  (S1);
  * the reference's candidate correction sign convention (s' - r f - e~ F) is wrong (S6/V2).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from end_golden import *

def sign(x):
    return (x > 0) - (x < 0)

def make_instance(ek, h_ntt, P, rng, gf_wraps=(), gf_strength=1.0, GF_wrap=False, GF_target=1.6):
    """Return (cprime, s, mubar, info). gf_wraps: list of coordinates where c1 should wrap by +1."""
    n = P.n
    f, g, F, G = ek['f'], ek['g'], ek['F'], ek['G']
    s = [rng.randrange(2) for _ in range(n)]
    sbar = [2 * x - 1 for x in s]
    mubar_choice = {1: (-3, -1, 1, 3), 0: (-2, 0, 2, 4)}[P.delta]     # valid honest mu~ values
    mubar = [rng.choice(mubar_choice) for _ in range(n)]
    if GF_wrap:
        # align s~ with the sign of G_Q (coordinate 0 of G_Q s~) on the largest entries until target*qQ
        GQ, FQ = ek['GQ'], ek['FQ']
        order = sorted(range(n), key=lambda j: -abs(GQ[j]))
        acc = 0
        for j in order:
            if acc >= GF_target * QQP:
                break
            if j == 0:
                sbar[0] = sign(GQ[0]) or 1
                acc += abs(GQ[0])
            else:
                sbar[n - j] = -(sign(GQ[j]) or 1)
                acc += abs(GQ[j])
        s = [(x + 1) // 2 for x in sbar]
    fixed = set()
    for t in gf_wraps:
        # (f mu~)_t = f_0 mu~_t + sum_{j>=1} f_j * (x^j mu~)_t  ;  (x^j mu~)_t = mu~_{t-j} (t>=j) or -mu~_{t-j+n}
        target = 2 * Q * gf_strength / 2.0 + Q * 0.25          # aim at 1.25 q .. keeps e1 = +1
        acc = 0.0
        for j in range(n):
            if f[j] == 0:
                continue
            if acc >= target:
                break
            i = t - j
            sgn = 1
            if i < 0:
                i += n; sgn = -1
            if i in fixed:
                continue
            want = sign(f[j]) * sgn            # want sgn * f_j * mu~_i > 0
            mubar[i] = max(mubar_choice) if want > 0 else min(mubar_choice)
            fixed.add(i)
            acc += abs(f[j]) * abs(mubar[i])
    # mu' from mu~ ; c' = (h s mod q) + mu' , folded into [0,256] by changing mu' where needed
    _, x = dpke_encrypt(h_ntt, s, P)
    muprime = [(m + P.delta) // 2 for m in mubar]
    cprime = []
    for i in range(n):
        v = x[i] + muprime[i]
        if v < 0 or v > 256:
            muprime[i] = 0; v = x[i]
        cprime.append(v)
    mubar = [2 * m - P.delta for m in muprime]
    return cprime, s, mubar

def true_wraps(ek, s, mubar, P):
    n = P.n
    sbar = [2 * x - 1 for x in s]
    v1 = pmul(ek['g'], sbar, n); t = pmul(ek['f'], mubar, n)
    v1 = [v1[i] + t[i] for i in range(n)]                      # g s~ + f mu~
    e1 = [round_div(v, 2 * Q) for v in v1]
    v2 = pmul(ek['GQ'], sbar, n); t = pmul(ek['FQ'], mubar, n); u = pmul(ek['w'], e1, n)
    v2 = [v2[i] + t[i] + 2 * Q * u[i] for i in range(n)]       # G_Q s~ + F_Q mu~ + 2 q w e1
    e2 = [round_div(v, TWO_QQP) for v in v2]
    e2p = [round_div(v - 2 * Q * u[i], TWO_QQP) for i, v in enumerate(v2)]   # wrap after removing e1
    return e1, e2, e2p

def check_identities(ek, cprime, s, mubar, P):
    n = P.n
    st = bnp_initial(ek, cprime, P)
    sbar = [2 * x - 1 for x in s]
    e1, e2, e2p = true_wraps(ek, s, mubar, P)
    # V5: (mu~', -s~') = (mu~, -s~) + 2 e1 (G,F) - 2 e2 (g,f)
    t1 = pmul(e1, ek['F'], n); t2 = pmul(e2, ek['f'], n)
    exp_s = [sbar[i] - 2 * t1[i] + 2 * t2[i] for i in range(n)]
    t3 = pmul(e1, ek['G'], n); t4 = pmul(e2, ek['g'], n)
    exp_mu = [mubar[i] + 2 * t3[i] - 2 * t4[i] for i in range(n)]
    assert st['sbar_p'] == exp_s, "V5 (s) violated"
    assert st['mubar_p'] == exp_mu, "V5 (mu) violated"
    return e1, e2, e2p

def ref_fcl_gf(c1, m0):
    """reference-style gf FCL (pke.c:298-316): score 128-|c1|, ranked by |score|"""
    sc = [abs(128 - abs(x)) for x in c1]
    return sorted(range(len(c1)), key=lambda i: (sc[i], i))[:m0]

def ref_sign_candidate(ek, st, cand, P):
    """reference-style gf correction in the {0,1} domain: s' - r f - e~ F  (pke.c:361-363)"""
    n = P.n
    sp = [(x + 1) // 2 for x in st['sbar_p']]
    we = apply_sparse(ek['w'], cand, n)
    r = [round_div(st['C2'][i] - 2 * Q * we[i], TWO_QQP) for i in range(n)]
    rf = pmul(r, ek['f'], n); eF = apply_sparse(ek['F'], cand, n)
    return [sp[i] - rf[i] - eF[i] for i in range(n)]

DUMP = []

def run(P, katfile, seed=1):
    rng = random.Random(seed)
    e = parse_kat(katfile)[0]
    ek = expanded_key(sk_from_entry(e), P)
    h_ntt = decode_pk(bytes.fromhex(e['pk']), P)
    n = P.n
    results = []

    # --- 1. honest instance (no wrap) -----------------------------------
    cprime, s, mubar = make_instance(ek, h_ntt, P, rng)
    e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
    assert not any(e1) and not any(e2), "unexpected wrap in honest instance"
    diag = {}
    res = decode_trap(ek, cprime, P, diag=diag)
    results.append(("honest, no wrap", res is not None and res[0] == s and diag['cand_gf'] == {} and diag['cand_GF'] == {}))
    DUMP.append(("honest no wrap", cprime, s, 1))

    # --- 2. single gf wrap at coordinate 0 -------------------------------
    for t in (0, 17, n - 1):
        cprime, s, mubar = make_instance(ek, h_ntt, P, rng, gf_wraps=[t])
        e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
        wraps = [i for i in range(n) if e1[i] != 0]
        diag = {}
        res = decode_trap(ek, cprime, P, diag=diag)
        ok = (wraps == [t]) and res is not None and res[0] == s and diag['cand_gf'] == {t: e1[t]}
        results.append(("single gf wrap at %d (e1=%+d, |e2'|=%d)" % (t, e1[t], sum(abs(x) for x in e2p)), ok))
        DUMP.append(("single gf wrap at %d" % t, cprime, s, 1 if ok else 0))
        if t == 0:
            # S1: reference FCL boundary 128 vs. correct boundary 257
            S_ref = ref_fcl_gf(diag['c1'], P.m0_gf)
            results.append(("S1: correct FCL contains wrap index", t in diag['S_gf']))
            results.append(("S1: reference 128-boundary FCL misses wrap index", t not in S_ref))
            # S6/V2: reference sign convention does not yield s for the true candidate
            st = bnp_initial(ek, cprime, P)
            s_ref = ref_sign_candidate(ek, st, {t: e1[t]}, P)
            s_true = [(x + 1) // 2 for x in res[1]]
            results.append(("S6: reference-sign correction (s'-rf-eF) != s", s_ref != s))
            results.append(("V2: corrected-sign correction == s (via decoder)", s_true == s))

    # --- 3. double gf wrap ------------------------------------------------
    cprime, s, mubar = make_instance(ek, h_ntt, P, rng, gf_wraps=[3, 200])
    e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
    wraps = sorted(i for i in range(n) if e1[i] != 0)
    diag = {}
    res = decode_trap(ek, cprime, P, diag=diag)
    ok = (wraps == [3, 200]) and res is not None and res[0] == s and set(diag['cand_gf']) == {3, 200}
    results.append(("double gf wrap at {3,200} (actual e1 support %s, |e2'| support %s, S_gf=%s, cand=%s, ok=%s)" % (
        wraps, [i for i in range(n) if e2p[i]], diag.get('S_gf'), diag.get('cand_gf'), res is not None and res[0] == s), ok))
    DUMP.append(("double gf wrap at 3,200", cprime, s, 1 if ok else 0))

    # --- 4. single GF wrap (no gf wrap) ------------------------------------
    for target in (1.3, 1.15):
        cprime, s, mubar = make_instance(ek, h_ntt, P, rng, GF_wrap=True, GF_target=target)
        e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
        wraps2 = [i for i in range(n) if e2[i] != 0]
        diag = {}
        res = decode_trap(ek, cprime, P, diag=diag)
        ok = (not any(e1)) and wraps2 == [0] and res is not None and res[0] == s and diag['cand_GF'] == {0: e2[0]}
        results.append(("single GF wrap at 0 (target %.2f qQ, e2 support %s, e1 support %s, S_GF=%s, cand_GF=%s)" % (
            target, [(i, e2[i]) for i in wraps2], [i for i in range(n) if e1[i]], diag.get('S_GF'), diag.get('cand_GF')), ok))
        DUMP.append(("single GF wrap at 0 target %.2f" % target, cprime, s, 1 if ok else 0))

    # --- 5. gf wrap + GF wrap -------------------------------------------
    cprime, s, mubar = make_instance(ek, h_ntt, P, rng, gf_wraps=[5], GF_wrap=True, GF_target=1.5)
    e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
    diag = {}
    res = decode_trap(ek, cprime, P, diag=diag)
    ok = res is not None and res[0] == s and diag['cand_gf'] == {5: e1[5]} and diag['cand_GF'] != {}
    results.append(("gf wrap at 5 + GF wrap (e2' support %s)" % [i for i in range(n) if e2p[i]], ok))
    DUMP.append(("gf wrap at 5 + GF wrap", cprime, s, 1 if ok else 0))

    # --- 6. FCL miss -> decryption failure is detected (3 gf wraps > B0) ----
    cprime, s, mubar = make_instance(ek, h_ntt, P, rng, gf_wraps=[1, 150, 300])
    e1, e2, e2p = check_identities(ek, cprime, s, mubar, P)
    nw = sum(1 for i in range(n) if e1[i])
    res = decode_trap(ek, cprime, P)
    results.append(("%d gf wraps (beyond B0=2): decoder cannot return the right s (returns %s)" % (
        nw, "None" if res is None else ("s" if res[0] == s else "wrong s")), nw > 2 and (res is None or res[0] != s)))
    DUMP.append(("%d gf wraps beyond B0" % nw, cprime, s, 0))
    return results

def dump(path, P):
    with open(path, 'w') as fp:
        for name, cprime, s, ok in DUMP:
            fp.write('name = %s' % name + chr(10))
            fp.write('cprime = %s' % b''.join(int(x).to_bytes(2, 'little') for x in cprime).hex() + chr(10))
            fp.write('s = %s' % bytes(s).hex() + chr(10))
            fp.write('expect_ok = %d' % ok + chr(10) + chr(10))

if __name__ == '__main__':
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'kat')
    do_dump = '--dump' in sys.argv
    allok = True
    for P, kat in ((P512, 'kat_end512_ref.txt'), (P1024, 'kat_end1024_ref.txt')):
        print("== END-%d ==" % P.n)
        del DUMP[:]
        for name, ok in run(P, os.path.join(base, kat)):
            print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
            allok &= ok
        if do_dump:
            dump(os.path.join(base, 'wraps_end%d.txt' % P.n), P)
            print("  dumped %d instances" % len(DUMP))
    print("ALL PASS" if allok else "SOME FAILURES")
    sys.exit(0 if allok else 1)
