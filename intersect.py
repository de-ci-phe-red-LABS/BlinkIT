"""Intersect per-column candidate lists from several same-key blink_attack logs."""
import re, sys
sets = {}
for f in sys.argv[1:]:
    txt = open(f).read()
    true = dict(re.findall(r"column (\d): +\d+ candidate\(s\) of 2\^16   true ([0-9a-f]{4})", txt))
    for c, lst in re.findall(r"CANDS (\d):((?: [0-9a-f]{4})*)", txt):
        s = set(lst.split()); sets.setdefault(c, []).append(s)
    print(f"{f}: " + "  ".join(f"col{c}={len(s[-1])}" for c, s in sorted(sets.items())))
prod = 1
for c in sorted(sets):
    inter = set.intersection(*sets[c]); prod *= len(inter)
    print(f"column {c}: {' ∩ '.join(str(len(s)) for s in sets[c])} -> {len(inter)}  true {true[c]} {'IN' if true[c] in inter else 'MISSING'}")
import math; print(f"product of intersections = {prod} = 2^{math.log2(prod):.2f}")
