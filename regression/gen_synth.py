#!/usr/bin/env python
"""Generate well-formed synthetic cgMLST/SNP profiles for regression testing.

Produces profiles with controlled properties:
  - clonal expansion (a founder mutated stepwise -> star + chains)
  - duplicate rows (to exercise `nonredundant` embeded groups)
  - missing-data gradient (0 / N / - markers)
  - integer allele identifiers (MLST) or nucleotide chars (SNP), TSV format

Deterministic: fixed seed per dataset name. Output = GrapeTree profile TSV.
"""
import random
import sys
import zlib


def gen(name, n_samples, n_loci, missing_rate=0.0, n_founders=3,
        mut_rate=0.05, dup_frac=0.0, alphabet="int", missing_char="0",
        seed=0):
    # Deterministic seed: Python's builtin hash() is salted per-process
    # (PYTHONHASHSEED), which would make datasets differ every run. crc32 of
    # the name is stable across processes and machines.
    rng = random.Random((zlib.crc32(name.encode()) ^ seed) & 0xFFFFFFFF)

    def new_allele():
        if alphabet == "int":
            return str(rng.randint(1, 40))
        else:  # nucleotide SNP-style
            return rng.choice("ACGT")

    # founders
    founders = [[new_allele() for _ in range(n_loci)] for _ in range(n_founders)]
    rows = []
    for i in range(n_samples):
        base = list(rng.choice(founders))
        # stepwise mutations
        n_mut = 0
        for k in range(n_loci):
            if rng.random() < mut_rate:
                base[k] = new_allele()
                n_mut += 1
        rows.append(base)

    # inject exact duplicates
    n_dup = int(dup_frac * n_samples)
    for _ in range(n_dup):
        src = rng.randrange(len(rows))
        rows.append(list(rows[src]))

    # inject missing data
    if missing_rate > 0:
        for r in rows:
            for k in range(n_loci):
                if rng.random() < missing_rate:
                    r[k] = missing_char

    # write TSV with a '#'-prefixed header
    header = ["#Strain"] + [f"L{k+1}" for k in range(n_loci)]
    out = ["\t".join(header)]
    for i, r in enumerate(rows):
        out.append("\t".join([f"S{i}"] + r))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    # gen.py <name> <n_samples> <n_loci> [missing_rate] [dup_frac] [alphabet] [missing_char]
    name = sys.argv[1]
    n_samples = int(sys.argv[2])
    n_loci = int(sys.argv[3])
    missing_rate = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    dup_frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    alphabet = sys.argv[6] if len(sys.argv) > 6 else "int"
    missing_char = sys.argv[7] if len(sys.argv) > 7 else "0"
    gen(name, n_samples, n_loci, missing_rate=missing_rate,
        dup_frac=dup_frac, alphabet=alphabet, missing_char=missing_char)
