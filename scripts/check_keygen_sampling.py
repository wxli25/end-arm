#!/usr/bin/env python3
"""Compare END-512 sampling paths using the 100 public KAT seeds."""

import argparse
import csv
import ctypes as C
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SAMPLE_BYTES = 384
PAIRS = [(i, i + 136) for i in range(200, 248)]


class KeccakState(C.Structure):
    _fields_ = [("s", C.c_uint64 * 25), ("pos", C.c_uint)]


def load_keccak(compiler, directory):
    suffix = ".dll" if sys.platform == "win32" else ".so"
    library = Path(directory) / ("keccak" + suffix)
    subprocess.run(
        [compiler, "-std=c99", "-O2", "-shared", "-fPIC",
         str(ROOT / "src/common/fips202.c"), "-o", str(library)],
        check=True,
    )
    keccak = C.CDLL(str(library))
    for name in ("shake256_absorb_once", "shake256_squeeze",
                 "shake128_squeezeblocks", "shake128_squeeze"):
        fn = getattr(keccak, name)
        if name.endswith("absorb_once"):
            fn.argtypes = [C.POINTER(KeccakState), C.c_void_p, C.c_size_t]
        else:
            fn.argtypes = [C.c_void_p, C.c_size_t, C.POINTER(KeccakState)]
        fn.restype = None
    return keccak


def sample_bytes(keccak, seed, legacy):
    state = KeccakState()
    stream = C.create_string_buffer(SAMPLE_BYTES)
    k = C.create_string_buffer(32)
    keccak.shake256_absorb_once(C.byref(state), seed, len(seed))
    keccak.shake256_squeeze(k, len(k), C.byref(state))
    if legacy:
        keccak.shake128_squeezeblocks(stream, 2, C.byref(state))
        keccak.shake128_squeeze(C.byref(stream, 336), 48, C.byref(state))
    else:
        keccak.shake256_squeeze(stream, len(stream), C.byref(state))
    return k.raw, stream.raw


def ternary(buf):
    def coefficient(x):
        return (x > 54) - (x < 9)

    result = []
    for offset in range(0, len(buf), 3):
        x, y, z = buf[offset:offset + 3]
        high_bits = ((x >> 6) << 4) | ((y >> 6) << 2) | (z >> 6)
        result.extend(coefficient(i) for i in (x & 63, y & 63, z & 63, high_bits))
    return result


def coefficient_index(byte):
    group, lane = divmod(byte, 3)
    return 4 * group + lane


def equal_coefficients(poly):
    return sum(poly[coefficient_index(a)] == poly[coefficient_index(b)] for a, b in PAIRS)


def read_kat(name):
    path = ROOT / "kat" / name
    entries = {}
    for block in path.read_text(encoding="ascii").strip().split("\n\n"):
        fields = dict(line.split(" = ", 1) for line in block.splitlines())
        entries[fields["keygen_seed"]] = fields
    return entries


def measure(keccak):
    old_kat = read_kat("kat_end512_ref.txt")
    new_kat = read_kat("kat_end512_shake256.txt")
    if old_kat.keys() != new_kat.keys():
        raise ValueError("The KAT seed sets differ")

    for seed_hex, old_entry in old_kat.items():
        seed = bytes.fromhex(seed_hex)
        new_entry = new_kat[seed_hex]
        old_k, old = sample_bytes(keccak, seed, legacy=True)
        new_k, new = sample_bytes(keccak, seed, legacy=False)
        standard = hashlib.shake_256(seed).digest(32 + SAMPLE_BYTES)
        if new_k + new != standard or old_k.hex() != old_entry["k"] or new_k.hex() != new_entry["k"]:
            raise ValueError(f"SHAKE or KAT mismatch for seed {seed_hex}")
        old_f = bytes.fromhex(old_entry["f"])
        new_f = bytes.fromhex(new_entry["f"])
        yield {
            "kat": old_entry["count"],
            "seed_hex": seed_hex,
            "old_first_equal_bytes": sum(old[a] == old[b] for a, b in PAIRS),
            "new_first_equal_bytes": sum(new[a] == new[b] for a, b in PAIRS),
            "old_first_equal_coefficients": equal_coefficients(ternary(old)),
            "new_first_equal_coefficients": equal_coefficients(ternary(new)),
            "old_final_f_equal_coefficients": equal_coefficients(old_f),
            "new_final_f_equal_coefficients": equal_coefficients(new_f),
            "old_first_left_8_hex": old[200:208].hex(),
            "old_first_right_8_hex": old[336:344].hex(),
            "new_first_left_8_hex": new[200:208].hex(),
            "new_first_right_8_hex": new[336:344].hex(),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", default="gcc", help="C compiler for fips202.c")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        keccak = load_keccak(args.cc, directory)
        rows = list(measure(keccak))
        if sys.platform == "win32":
            # Release the DLL before TemporaryDirectory removes it on Windows.
            C.windll.kernel32.FreeLibrary.argtypes = [C.c_void_p]
            C.windll.kernel32.FreeLibrary(keccak._handle)
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
