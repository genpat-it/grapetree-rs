#!/usr/bin/env python
"""Standalone runner for GrapeTree's backend, bypassing the Flask package __init__.

Usage: oracle_run.py <MSTrees.py path> <profile> <method> [key=value ...]
Prints the backend() output (newick or PHYLIP matrix) to stdout.
"""
import sys, importlib.util, os

def load_backend(mstrees_path):
    spec = importlib.util.spec_from_file_location("gt_mstrees", mstrees_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gt_mstrees"] = mod
    spec.loader.exec_module(mod)
    return mod.backend

def main():
    mstrees_path = sys.argv[1]
    profile = sys.argv[2]
    method = sys.argv[3]
    kwargs = {}
    for kv in sys.argv[4:]:
        k, v = kv.split("=", 1)
        if v in ("True", "False"):
            v = (v == "True")
        elif v.lstrip("-").isdigit():
            v = int(v)
        else:
            try:
                v = float(v)
            except ValueError:
                pass
        kwargs[k] = v
    backend = load_backend(mstrees_path)
    out = backend(profile=profile, method=method, **kwargs)
    sys.stdout.write(out)
    if not out.endswith("\n"):
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
