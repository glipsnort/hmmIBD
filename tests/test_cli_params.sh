#!/usr/bin/env bash
# Tests for the hmmIBD command-line parameters: every option is parsed, reported at startup,
# and changes the computation in the expected way. Uses the bundled samp_data.
#
#   tests/test_cli_params.sh            # builds ./hmmIBD from hmmIBD.c, runs all tests
#   HMMIBD=./hmmIBD tests/test_cli_params.sh   # use an existing binary
#
# Requires: cc (unless HMMIBD is given), bash, awk, sort, diff. Exit status 0 iff all tests pass.
set -u
cd "$(dirname "$0")/.."
D=samp_data; G=$D/pf3k_Cambodia_13.txt; F=$D/freqs_pf3k_Cambodia_13.txt
G2=$D/pf3k_Ghana_13.txt; F2=$D/freqs_pf3k_Ghana_13.txt
T=$(mktemp -d "${TMPDIR:-/tmp}/hmmibd_test.XXXXXX"); trap 'rm -rf "$T"' EXIT
pass=0; fail=0
ok()   { pass=$((pass+1)); printf 'PASS  %s\n' "$1"; }
bad()  { fail=$((fail+1)); printf 'FAIL  %s\n      %s\n' "$1" "$2"; }
check(){ # check <name> <condition-exit-status> <detail>
  if [ "$2" -eq 0 ]; then ok "$1"; else bad "$1" "$3"; fi; }
run()  { # run <prefix> <args...>  -> stdout in $T/<prefix>.out, exit status in $T/<prefix>.rc
  local p=$1; shift; "$HMMIBD" "$@" -o "$T/$p" > "$T/$p.out" 2> "$T/$p.err"; echo $? > "$T/$p.rc"; }
rc()   { cat "$T/$1.rc"; }
frac() { sed 1d "$T/$1.hmm_fract.txt"; }
segs() { sed 1d "$T/$1.hmm.txt"; }
col()  { frac "$1" | cut -f"$2"; }          # frac column: 3 N_inform 4 discord 5 log_p 6 iter 7 k_rec 10 fract_IBD
same() { cmp -s "$1" "$2"; }
nl()   { wc -l < "$1" | tr -d ' '; }

# ---------------------------------------------------------------- build
if [ -z "${HMMIBD:-}" ]; then
  HMMIBD=$T/hmmIBD
  cc -O2 -Wall -Wextra -o "$HMMIBD" hmmIBD.c -lm > "$T/build.log" 2>&1
  check "build: compiles" $? "$(head -5 "$T/build.log")"
  nw=$(grep -c 'warning:' "$T/build.log")
  check "build: no compiler warnings with -Wall -Wextra" $([ "$nw" -eq 0 ]; echo $?) "$nw warnings: $(grep 'warning:' "$T/build.log" | head -3)"
fi

# ---------------------------------------------------------------- help
"$HMMIBD" -h > "$T/help.txt" 2>&1; check "help: -h exits 0" $? ""
missing=""; for o in i o I f F b g m n r N e R x c s k d D t T u h; do grep -q "^  -$o " "$T/help.txt" || missing="$missing -$o"; done
check "help: lists every option" $([ -z "$missing" ]; echo $?) "missing:$missing"
nodef=""; for o in m n r N e R x c s k d D t T u f F b g; do grep "^  -$o " "$T/help.txt" | grep -q '\[.*\]$' || nodef="$nodef -$o"; done
check "help: every optional parameter shows a default in brackets" $([ -z "$nodef" ]; echo $?) "no default:$nodef"
# defaults in help must equal the values the program reports when run with no options
run def -i "$G" -f "$F"
hd() { grep "^  -$1 " "$T/help.txt" | sed 's/.*\[\(.*\)\]$/\1/'; }
check "help: default for -m matches startup report" $([ "$(hd m)" = "$(grep -o 'iterations allowed: [0-9]*' "$T/def.out" | awk '{print $3}')" ]; echo $?) "help [$(hd m)]"
check "help: default for -s matches startup report" $([ "$(hd s)" = "$(grep -o 'spacing (bp): [0-9]*' "$T/def.out" | awk '{print $3}')" ]; echo $?) "help [$(hd s)]"
check "help: default for -k matches startup report" $([ "$(hd k)" = "$(grep -o 'informative markers: [0-9]*' "$T/def.out" | awk '{print $3}')" ]; echo $?) "help [$(hd k)]"
check "help: default for -R matches startup report" $([ "$(hd R)" = "$(grep -o 'generation): .*' "$T/def.out" | awk '{print $2}')" ]; echo $?) "help [$(hd R)]"
check "help: default for -x matches startup report" $([ "$(hd x)" = "$(grep -o 'allele index: [0-9]*' "$T/def.out" | awk '{print $3}')" ]; echo $?) "help [$(hd x)]"
check "help: default for -c matches startup report" $([ "$(hd c)" = "$(grep -o 'chromosomes: [0-9]*' "$T/def.out" | awk '{print $2}')" ]; echo $?) "help [$(hd c)]"

# ---------------------------------------------------------------- defaults
check "defaults: run succeeds" "$(rc def)" "$(cat "$T/def.err")"
run expl -i "$G" -f "$F" -m 5 -N 1 -e 0.001 -R 7.4e-7 -x 8 -c 14 -s 5 -k 10 -d 0 -D 1 -t 0.001 -T 0.01 -u 0.001
check "defaults: explicit defaults give identical output" $(same "$T/def.hmm.txt" "$T/expl.hmm.txt" && same "$T/def.hmm_fract.txt" "$T/expl.hmm_fract.txt"; echo $?) ""
if [ -f "$D/output_Cambodia.hmm.txt" ]; then
  check "defaults: identical to shipped reference output" $(same "$T/def.hmm.txt" "$D/output_Cambodia.hmm.txt" && same "$T/def.hmm_fract.txt" "$D/output_Cambodia.hmm_fract.txt"; echo $?) ""
fi

# ---------------------------------------------------------------- propagation: values reach the startup report
run rep -i "$G" -f "$F" -m 3 -s 7 -k 12 -d 0.1 -D 0.9 -c 3 -x 4 -R 1e-6 -N 2.5 -t 0.002 -T 0.02 -u 0.003 -e 0.02 -n 7 -r 0.3
check "report: run with every option succeeds" "$(rc rep)" "$(cat "$T/rep.err")"
for want in "Maximum fit iterations allowed: 3" "Minimum marker spacing (bp): 7" "Minimum informative markers: 12" \
            "Discordance accepted in range \[0.1, 0.9\]" "Number of chromosomes: 3" "Maximum allele index: 4" \
            "Recombination rate (per bp per generation): 1e-06" "Initial N generations: 2.5" \
            "Fit thresholds: dpi 0.002, dk 0.02, drelk 0.003" "Genotyping error rate: 2.00%" \
            "Number of generations capped at 7.00" "IBD fract fixed at 0.30"; do
  grep -q "$want" "$T/rep.out"; check "report: '$want'" $? "not found in startup output"
done

# ---------------------------------------------------------------- effects: each parameter changes the computation as specified
# -m
run m1 -i "$G" -f "$F" -m 1; run m2 -i "$G" -f "$F" -m 2
check "-m 1: every pair reports N_fit_iteration = 1" $([ "$(col m1 6 | sort -u)" = "1" ]; echo $?) "$(col m1 6 | sort -u | tr '\n' ' ')"
check "-m 2: no pair exceeds 2 iterations, output differs from default" $([ "$(col m2 6 | sort -n | tail -1)" -le 2 ] && ! same "$T/m2.hmm_fract.txt" "$T/def.hmm_fract.txt"; echo $?) ""
# -n
run n05 -i "$G" -f "$F" -n 0.5
check "-n 0.5: every N_generation <= 0.5" $([ "$(col n05 7 | awk '$1>0.5' | wc -l | tr -d ' ')" -eq 0 ]; echo $?) "max $(col n05 7 | sort -g | tail -1)"
check "-n 0.5: output differs from default" $(! same "$T/n05.hmm_fract.txt" "$T/def.hmm_fract.txt"; echo $?) ""
# -r
run r03 -i "$G" -f "$F" -r 0.3
check "-r 0.3: output differs from default" $(! same "$T/r03.hmm.txt" "$T/def.hmm.txt"; echo $?) ""
# -e
run e01 -i "$G" -f "$F" -e 0.01
check "-e 0.01: log_p differs from default for every pair" $([ "$(paste <(col def 5) <(col e01 5) | awk '$1==$2' | wc -l | tr -d ' ')" -eq 0 ]; echo $?) ""
# -N
run N3 -i "$G" -f "$F" -N 3 -m 1
check "-N 3 (with -m 1): N_generation differs from -m 1 alone" $(! same <(col N3 7) <(col m1 7); echo $?) ""
# -R / -N / -T exact scale invariance: ptrans = k * R * d, so doubling R, halving N0 and the absolute k threshold
# must leave the HMM unchanged and halve the fitted N_generation
run inv -i "$G" -f "$F" -R 1.48e-6 -N 0.5 -T 0.005
check "-R x2, -N /2, -T /2: segments identical to default" $(same "$T/inv.hmm.txt" "$T/def.hmm.txt"; echo $?) ""
check "-R x2, -N /2, -T /2: all frac columns except N_generation identical" $(same <(frac def | cut -f1-6,8-11) <(frac inv | cut -f1-6,8-11); echo $?) ""
nbad=$(paste <(col def 7) <(col inv 7) | awk '$1>0.01 && !($2/$1>0.498 && $2/$1<0.502)' | wc -l | tr -d ' ')
check "-R x2, -N /2, -T /2: N_generation exactly halved (pairs with k>0.01)" $([ "$nbad" -eq 0 ]; echo $?) "$nbad pairs off"
# -x
run x2 -i "$G" -f "$F" -x 2
check "-x 2: variants with too many alleles are dropped" $(grep -q 'Variants with too many alleles: [1-9]' "$T/x2.out" && [ "$(grep -o '^[0-9]* variants used' "$T/x2.out" | awk '{print $1}')" -lt "$(grep -o '^[0-9]* variants used' "$T/def.out" | awk '{print $1}')" ]; echo $?) "$(grep 'variants' "$T/x2.out" | tr '\n' ' ')"
# -c
run c2 -i "$G" -f "$F" -c 2
check "-c 2: segments only on chromosomes <= 2 (default has higher)" $([ "$(segs c2 | awk -F'\t' '$3>2' | wc -l | tr -d ' ')" -eq 0 ] && [ "$(segs def | awk -F'\t' '$3>2' | wc -l | tr -d ' ')" -gt 0 ]; echo $?) ""
# -s
run s0 -i "$G" -f "$F" -s 0; run s100 -i "$G" -f "$F" -s 100
vu() { grep -o '^[0-9]* variants used' "$T/$1.out" | awk '{print $1}'; }
check "-s 0: no variants skipped for spacing, more variants used than default" $(grep -q 'skipped for spacing: 0$' "$T/s0.out" && [ "$(vu s0)" -gt "$(vu def)" ]; echo $?) "$(vu s0) vs $(vu def)"
check "-s 100: fewer variants used than default" $([ "$(vu s100)" -lt "$(vu def)" ]; echo $?) "$(vu s100) vs $(vu def)"
# -k
run k1000 -i "$G" -f "$F" -k 1000
check "-k 1000: no pair analyzed, output files have header only" $(grep -q 'markers): 0$' "$T/k1000.out" && [ "$(nl "$T/k1000.hmm_fract.txt")" -eq 1 ] && [ "$(nl "$T/k1000.hmm.txt")" -eq 1 ]; echo $?) ""
head -41 "$G" > "$T/g40.txt"; head -40 "$F" > "$T/f40.txt"
run k0 -i "$T/g40.txt" -f "$T/f40.txt" -k 0; run k10 -i "$T/g40.txt" -f "$T/f40.txt"
check "-k 0 on 40 variants: all 45 pairs reported (default -k 10 reports none)" $([ "$(frac k0 | nl /dev/stdin)" -eq 45 ] && [ "$(frac k10 | nl /dev/stdin)" -eq 0 ]; echo $?) "$(frac k0 | nl /dev/stdin) / $(frac k10 | nl /dev/stdin)"
# -d / -D: pick a threshold strictly between two observed discordance values, then the kept set is exactly the default rows on one side of it
thr=$(col def 4 | sort -g | awk 'NR==20{a=$1} NR==21{b=$1} END{printf "%.5f", (a+b)/2}')
# (prefixes differ in more than letter case: output paths must be distinct on case-insensitive filesystems)
run dmin -i "$G" -f "$F" -d "$thr"; run dmax -i "$G" -f "$F" -D "$thr"
check "-d $thr: kept pairs == default pairs with discordance >= threshold" $(same <(frac def | awk -F'\t' -v t="$thr" '$4>=t') <(frac dmin); echo $?) "$(frac dmin | nl /dev/stdin) kept"
check "-D $thr: kept pairs == default pairs with discordance <= threshold" $(same <(frac def | awk -F'\t' -v t="$thr" '$4<=t') <(frac dmax); echo $?) "$(frac dmax | nl /dev/stdin) kept"
check "-d/-D: the two subsets partition the default pairs" $([ $(( $(frac dmin | nl /dev/stdin) + $(frac dmax | nl /dev/stdin) )) -eq "$(frac def | nl /dev/stdin)" ]; echo $?) ""
# -t / -T / -u
run never -i "$G" -f "$F" -t 0 -T 0 -u 0
check "-t 0 -T 0 -u 0: fit never converges, every pair runs all 5 iterations" $([ "$(col never 6 | sort -u)" = "5" ]; echo $?) "$(col never 6 | sort -u | tr '\n' ' ')"
# |delta pi| <= 1 always; |delta k| can be hundreds, so the k thresholds must be large to be "always met"
run conv_abs -i "$G" -f "$F" -t 1 -T 1e9
check "-t 1 -T 1e9: fit converges after the first update for every pair" $([ "$(col conv_abs 6 | sort -u)" = "1" ]; echo $?) "$(col conv_abs 6 | sort -u | tr '\n' ' ')"
run conv_rel -i "$G" -f "$F" -t 1 -u 1e9
check "-t 1 -u 1e9: relative criterion alone also suffices" $([ "$(col conv_rel 6 | sort -u)" = "1" ]; echo $?) "$(col conv_rel 6 | sort -u | tr '\n' ' ')"
run pi_only -i "$G" -f "$F" -t 1
check "-t 1 alone: N criteria still apply, some pairs need > 1 iteration" $([ "$(col pi_only 6 | sort -n | tail -1)" -gt 1 ]; echo $?) ""
run k_only -i "$G" -f "$F" -T 1e9 -u 1e9
check "-T 1e9 -u 1e9 alone: pi criterion still applies, some pairs need > 1 iteration" $([ "$(col k_only 6 | sort -n | tail -1)" -gt 1 ]; echo $?) ""

# ---------------------------------------------------------------- existing options still behave
run twopop -i "$G" -I "$G2" -f "$F" -F "$F2"; run twopop_e -i "$G" -I "$G2" -f "$F" -F "$F2" -e 0.001 -R 7.4e-7
check "two-population run unchanged by explicit defaults" $(same "$T/twopop.hmm.txt" "$T/twopop_e.hmm.txt" && same "$T/twopop.hmm_fract.txt" "$T/twopop_e.hmm_fract.txt"; echo $?) ""
if [ -f "$D/output_Cambodia_Ghana.hmm.txt" ]; then
  check "two-population run identical to shipped reference output" $(same "$T/twopop.hmm.txt" "$D/output_Cambodia_Ghana.hmm.txt" && same "$T/twopop.hmm_fract.txt" "$D/output_Cambodia_Ghana.hmm_fract.txt"; echo $?) ""
fi
s1=$(frac def | head -1 | cut -f1); s2=$(frac def | head -1 | cut -f2)
printf '%s\t%s\n' "$s2" "$s1" > "$T/good.txt"; run good -i "$G" -f "$F" -g "$T/good.txt"
check "-g: one listed pair (reverse order) gives exactly that default row" $(same <(frac def | head -1) <(frac good); echo $?) ""
printf '%s\n%s\n' "$s1" "$s2" > "$T/bad.txt"; run badf -i "$G" -f "$F" -b "$T/bad.txt"
check "-b: excluded samples absent, remaining rows identical to default" $([ "$(frac badf | grep -c "$s1\|$s2")" -eq 0 ] && same <(frac def | grep -v "$s1" | grep -v "$s2") <(frac badf); echo $?) ""

# ---------------------------------------------------------------- validation
badarg() { # badarg <name> <args...>: must exit non-zero with a message on stderr and write no output files
  local name=$1; shift; run v "$@" -i "$G" -f "$F"
  check "reject: $name" $([ "$(rc v)" -ne 0 ] && [ -s "$T/v.err" ]; echo $?) "rc=$(rc v) err='$(head -1 "$T/v.err")'"; rm -f "$T/v".hmm*; }
badarg "-e 1.5"            -e 1.5
badarg "-e abc"            -e abc
badarg "-m 0"              -m 0
badarg "-m 3x (trailing junk)" -m 3x
badarg "-x 0"              -x 0
badarg "-c 0"              -c 0
badarg "-s -1"             -s -1
badarg "-k -1"             -k -1
badarg "-d 0.5 -D 0.2"     -d 0.5 -D 0.2
badarg "-D 1.5"            -D 1.5
badarg "-R 0"              -R 0
badarg "-N 0"              -N 0
badarg "-n 0"              -n 0
badarg "-r 1.0"            -r 1.0
badarg "-t -1"             -t -1
badarg "unknown option -z" -z 1
badarg "-m without value (last arg)" -m
"$HMMIBD" -i "$G" > "$T/noo.err" 2>&1; check "reject: missing -o" $([ $? -ne 0 ]; echo $?) ""

echo; echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
