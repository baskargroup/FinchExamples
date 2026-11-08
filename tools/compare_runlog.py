#!/usr/bin/env python3
"""Compare the summary numbers of two Finch run logs.

    python3 tools/compare_runlog.py REF NEW [--rtol 1e-8] [--atol 1e-12]

Every line of the form `name = v`, `name: a -> b, ...`, `dt = v, steps = n ...`
is a key (the text before the first `=` or `:`) with the floats found in it.
Exit 0 when every key of REF is in NEW with all its numbers within rtol of the
largest magnitude on that line (a balance error next to O(1) totals is compared
on the totals' scale), or within atol absolutely (round-off-level drifts);
exit 1 with one line per mismatch otherwise. Lines starting with `run time`
are ignored. Lines without a `=` or `:` (e.g. "61 snapshots written") use the
whole line, without its numbers, as the key.
"""
import argparse
import re
import sys

NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def parse(path):
    keys = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("run time"):
            continue
        m = re.search(r"[=:]", line)
        # the numbers are those after the key: "L1 error = ..." has one, not two
        key = line[: m.start()].rstrip() if m else NUM.sub("#", line).strip()
        nums = [float(x) for x in NUM.findall(line[m.end():] if m else line)]
        if m is None and not nums:
            continue
        keys[key] = nums
    return keys


def close(a, b, rtol, atol, scale):
    return abs(a - b) <= rtol * max(abs(a), abs(b), scale) or abs(a - b) <= atol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--atol", type=float, default=1e-12)
    a = ap.parse_args()
    ref, new = parse(a.ref), parse(a.new)
    bad, count = [], 0
    for key, rv in ref.items():
        if key not in new:
            bad.append(f"{key!r}  missing in NEW")
            continue
        nv = new[key]
        if len(nv) != len(rv):
            bad.append(f"{key!r}  {len(rv)} numbers in REF, {len(nv)} in NEW")
            continue
        scale = max([abs(v) for v in rv] + [0.0])
        for i, (x, y) in enumerate(zip(rv, nv)):
            count += 1
            if not close(x, y, a.rtol, a.atol, scale):
                rel = abs(x - y) / max(abs(x), abs(y), 1e-300)
                bad.append(f"{key!r}  [{i}]  ref {x!r}  new {y!r}  rel {rel:.3e}")
    if bad:
        print("\n".join(bad))
        return 1
    print(f"OK ({count} numbers, {len(ref)} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
