#!/bin/bash
# Tests for tools/compare_runlog.py. Run: bash tools/test_compare_runlog.sh
cd "$(dirname "$0")/.."
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cat > "$T/ref.txt" <<'REF'
dt = 0.0002499974259809647, steps = 25133 (one revolution, T = 2*pi)
max error = 0.6398676897668724
L1 error = 0.003385368754585961
min phi = 9.909731958264485e-14, max phi = 0.9999987571495362
relative mass change = -5.052801217160516e-13
61 snapshots written
REF
fail=0
check() { # name expected_exit grep_word_or_empty
  local name=$1 want=$2 word=$3
  out=$(python3 tools/compare_runlog.py "$T/ref.txt" "$T/new.txt" 2>&1); got=$?
  if [ "$got" != "$want" ]; then echo "FAIL $name: exit $got, wanted $want"; echo "$out"; fail=1; return; fi
  if [ -n "$word" ] && ! echo "$out" | grep -q "$word"; then echo "FAIL $name: output lacks '$word'"; echo "$out"; fail=1; return; fi
  echo "pass $name"
}
cp "$T/ref.txt" "$T/new.txt";                                             check identical 0 OK
sed 's/0.003385368754585961/0.0033853687/' "$T/ref.txt" > "$T/new.txt";   check "changed L1" 1 L1
{ cat "$T/ref.txt"; echo "run time = 1 s"; } > "$T/new.txt";              check "extra run time line" 0 OK
grep -v "max phi" "$T/ref.txt" > "$T/new.txt";                            check "missing line" 1 "min phi"
{ cat "$T/ref.txt"; echo "run time = 1 s"; } > "$T/new.txt"
# a run time line in REF is ignored too
cp "$T/new.txt" "$T/ref.txt"; echo "run time = 99 s" >> "$T/new.txt";     check "run time ignored" 0 OK
# tolerance scales with the largest number on the line: a near-zero balance error
# next to O(1) volumes may differ in its 7th digit, a lone small number may not
printf 'gas volume: 0.25 -> 0.41079902933452633, created 0.16079863239518402, balance error 3.969393423397527e-7\nmax error = 1.0e-7\n' > "$T/ref.txt"
printf 'gas volume: 0.25 -> 0.41079902933452633, created 0.16079863239518402, balance error 3.969388275293362e-7\nmax error = 1.0e-7\n' > "$T/new.txt"; check "scaled tolerance" 0 OK
printf 'gas volume: 0.25 -> 0.41079902933452633, created 0.16079863239518402, balance error 3.969393423397527e-7\nmax error = 1.00001e-7\n' > "$T/new.txt"; check "lone small number" 1 "max error"
# an absolute floor for round-off-level numbers (a mass drift of 1e-9 +- 1e-14)
printf 'relative mass change = -7.475342541264091e-9\n' > "$T/ref.txt"
printf 'relative mass change = -7.475355312345479e-9\n' > "$T/new.txt"; check "absolute floor" 0 OK
printf 'relative mass change = -7.4e-9\n' > "$T/new.txt";              check "above the floor" 1 "mass"
[ $fail = 0 ] && echo ALL PASS || { echo SOME FAILED; exit 1; }
