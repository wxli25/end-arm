#!/usr/bin/env python3
"""Turn run_m4model.sh logs into Cortex-M4 comparison tables.

    python3 pqm4/m4model/m4model_report.py <log-dir> [--board FILE[=LABEL] ...]
                                           [--csv out.csv] [--profile scheme/impl ...] [--top 12]

Per implementation (logs <name>_{none,kg10,keypair,enc,dec}.log):
    keygen  = (kg10 - none) / 10          (average over ten seeds; falls back to keypair - none)
    op2     = enc - keypair,   op3 = dec - enc
for executed instructions and modelled cycles.  For KEMs op2/op3 are encaps/decaps, for signature
schemes (pqm4/icount_sign) sign/verify; sign is a single run of a variable-time algorithm.

--board takes pqm4 speed tables with measured board cycles (results/speed/<family>.md, or the
legacy benchmarks.md); when several files list an implementation, the first one wins.  Rows with
a measurement get the model error next to them.  The hashing share counts the symmetric
primitives of pqm4's mupq/common (Keccak/FIPS 202, SHA-2, AES, seed expanders), the functions
pqm4's hashing benchmark times.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

OPS = ("keygen", "op2", "op3")
LABELS = {"kem": ("KeyGen", "Encaps", "Decaps"), "sign": ("KeyGen", "Sign", "Verify")}
HASH_RE = re.compile(r"keccak|shake|sha3|sha2|sha256|sha512|hashblocks|aes|seedexpander|cshake|kmac|fips202", re.I)
RNG_RE = re.compile(r"^(randombytes|surf)$")
PREFIXES = {"mupq_pqclean_crypto_kem_": "kem", "mupq_crypto_kem_": "kem", "crypto_kem_": "kem",
            "mupq_pqclean_crypto_sign_": "sign", "mupq_crypto_sign_": "sign", "crypto_sign_": "sign"}


def parse_log(path):
    d = {"insns": 0, "cycles": 0, "stat": {}, "sym": {}}
    for line in path.read_text().splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] in ("insns", "cycles"):
            d[f[0]] = int(f[1])
        elif f[0] == "stat":
            d["stat"][f[1]] = int(f[2])
        elif f[0] == "sym":
            d["sym"][f[1]] = (int(f[2]), int(f[3]))
    return d


def split_name(name):
    for p, kind in PREFIXES.items():
        if name.startswith(p):
            scheme, impl = name[len(p):].rsplit("_", 1)
            return kind, scheme, impl
    raise ValueError(name)


def diff(a, b, scale=1.0):
    """(a - b) / scale for totals, stats and per-symbol profiles."""
    out = {"insns": (a["insns"] - b["insns"]) / scale, "cycles": (a["cycles"] - b["cycles"]) / scale,
           "stat": {k: (a["stat"].get(k, 0) - b["stat"].get(k, 0)) / scale for k in a["stat"]}, "sym": {}}
    for s in set(a["sym"]) | set(b["sym"]):
        ai, ac = a["sym"].get(s, (0, 0))
        bi, bc = b["sym"].get(s, (0, 0))
        if ai != bi or ac != bc:
            out["sym"][s] = ((ai - bi) / scale, (ac - bc) / scale)
    return out


def load(logdir):
    runs = defaultdict(dict)
    for p in logdir.glob("*.log"):
        m = re.match(r"(.+)_(none|kg10|keypair|enc|dec)\.log$", p.name)
        if m:
            runs[m.group(1)][m.group(2)] = parse_log(p)
    res = {}
    for name, r in runs.items():
        if not {"none", "keypair", "enc", "dec"} <= set(r):
            print("incomplete runs for", name, file=sys.stderr)
            continue
        kind, scheme, impl = split_name(name)
        kg = diff(r["kg10"], r["none"], 10.0) if "kg10" in r else diff(r["keypair"], r["none"])
        res[(scheme, impl)] = {"kind": kind, "keygen": kg, "op2": diff(r["enc"], r["keypair"]),
                               "op3": diff(r["dec"], r["enc"])}
    return res


def share(op, regex):
    tot = op["cycles"]
    part = sum(c for s, (_, c) in op["sym"].items() if regex.search(s))
    return part / tot if tot else 0.0


def parse_board(path, label, out):
    """Speed tables of a pqm4 results page: {(scheme, impl): (label, (kg, op2, op3))}."""
    top = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("# "):
            top = line[2:].strip()
        if not (top or "").startswith("Speed Evaluation") or not line.startswith("| ") or "AVG" not in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        scheme = cells[0].split(" (")[0]
        avg = tuple(int(re.search(r"AVG: ([\d,]+)", c).group(1).replace(",", "")) for c in cells[2:5])
        out.setdefault((scheme, cells[1]), (label, avg))


def k(x):
    return f"{x / 1000:,.0f}k" if x < 10e6 else f"{x / 1e6:,.2f}M"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logdir", type=Path)
    ap.add_argument("--board", action="append", default=[], help="pqm4 speed table FILE[=LABEL]")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--profile", nargs="*", default=[], help="scheme/impl to print a per-function profile for")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    res = load(args.logdir)
    board = {}
    for spec in args.board:
        path, _, label = spec.partition("=")
        parse_board(path, label or Path(path).name, board)

    rows = []
    for kind in ("kem", "sign"):
        items = sorted((key, r) for key, r in res.items() if r["kind"] == kind)
        if not items:
            continue
        L = LABELS[kind]
        print(f"## Modelled Cortex-M4 cycles and executed instructions ({'KEM' if kind == 'kem' else 'signatures'})\n")
        print(f"| scheme | impl | {L[0]} cyc | {L[1]} cyc | {L[2]} cyc | {L[0]} insn | {L[1]} insn | {L[2]} insn "
              f"| CPI ({L[2]}) | hash % ({'/'.join(L)}) |")
        print("|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for (scheme, impl), r in items:
            cpi = r["op3"]["cycles"] / r["op3"]["insns"]
            hs = "/".join(f"{100 * share(r[o], HASH_RE):.1f}" for o in OPS)
            print(f"| {scheme} | {impl} | " + " | ".join(f"{r[o]['cycles']:,.0f}" for o in OPS) + " | "
                  + " | ".join(f"{r[o]['insns']:,.0f}" for o in OPS) + f" | {cpi:.3f} | {hs} |")
            b = board.get((scheme, impl))
            row = {"kind": kind, "scheme": scheme, "impl": impl, "board_source": b[0] if b else ""}
            for o, lab in zip(OPS, L):
                lab = lab.lower()
                row[lab + "_cycles"] = round(r[o]["cycles"])
                row[lab + "_insns"] = round(r[o]["insns"])
                row[lab + "_hash_pct"] = round(100 * share(r[o], HASH_RE), 2)
                row[lab + "_rng_pct"] = round(100 * share(r[o], RNG_RE), 3)
                row[lab + "_board"] = b[1][OPS.index(o)] if b else ""
            rows.append(row)
        print()

    if board:
        print("## Model vs. measured board cycles\n")
        print("| scheme | impl | source | op 1 model / board | err | op 2 model / board | err | op 3 model / board | err |")
        print("|---|---|---|---|---:|---|---:|---|---:|")
        errs = defaultdict(list)
        for (scheme, impl), r in sorted(res.items()):
            b = board.get((scheme, impl))
            if not b:
                continue
            cells = []
            for i, o in enumerate(OPS):
                if r["kind"] == "sign" and o == "op2":     # signing time varies with rejections
                    cells.append(f"{k(r[o]['cycles'])} / {k(b[1][i])} (1 run) | n/a")
                    continue
                e = r[o]["cycles"] / b[1][i] - 1
                errs[b[0]].append(e)
                cells.append(f"{k(r[o]['cycles'])} / {k(b[1][i])} | {100 * e:+.1f}%")
            print(f"| {scheme} | {impl} | {b[0]} | " + " | ".join(cells) + " |")
        print()
        for src, es in errs.items():
            mae = sum(abs(e) for e in es) / len(es)
            print(f"- {src}: {len(es)} points, mean |error| {100 * mae:.1f}%, "
                  f"range {100 * min(es):+.1f}% .. {100 * max(es):+.1f}%")

    for want in args.profile:
        scheme, impl = want.split("/")
        r = res.get((scheme, impl))
        if not r:
            print("no data for", want, file=sys.stderr)
            continue
        for o, lab in zip(OPS, LABELS[r["kind"]]):
            tot = r[o]["cycles"]
            print(f"\n### Profile {scheme}/{impl} {lab}: {tot:,.0f} cycles, {r[o]['insns']:,.0f} insns\n")
            print("| function | cycles | share | insns | CPI |")
            print("|---|---:|---:|---:|---:|")
            items = sorted(r[o]["sym"].items(), key=lambda kv: -kv[1][1])
            for s, (ins, cyc) in items[:args.top]:
                cpi = cyc / ins if ins else 0
                print(f"| `{s}` | {cyc:,.0f} | {100 * cyc / tot:.1f}% | {ins:,.0f} | {cpi:.2f} |")
            st = r[o]["stat"]
            print("\nclass counts: " + ", ".join(f"{kk} {v:,.0f}" for kk, v in st.items()))

    if args.csv and rows:
        fields = sorted({f for row in rows for f in row}, key=lambda f: (f not in ("kind", "scheme", "impl"), f))
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
