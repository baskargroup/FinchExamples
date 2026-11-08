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
[ $fail = 0 ] && echo ALL PASS || { echo SOME FAILED; exit 1; }
