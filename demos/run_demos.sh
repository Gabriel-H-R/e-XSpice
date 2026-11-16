#!/usr/bin/env bash
# Runs every showcase netlist through the CLI and rebuilds the README figures.
#
#   ./demos/run_demos.sh [path-to-e-XSpice] [path-to-csv2svg]
#
# Writes CSVs and console logs to demos/out and regenerates the SVG plots in
# demos/plots with the csv2svg tool. Safe to re-run at any time.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EXE="${1:-$ROOT/build/e-XSpice}"
CSV2SVG="${2:-$ROOT/build/csv2svg}"
OUT="$ROOT/demos/out"
PLOTS="$ROOT/demos/plots"

mkdir -p "$OUT" "$PLOTS"

for d in rc_step ac_rlc_resonance noise_divider sens_divider hb_clipper clipper_tran bsim4_inverter; do
    "$EXE" file "$ROOT/demos/$d.cir" "$OUT/$d.csv" > "$OUT/$d.log" 2>&1
    echo "  $d ok"
done

"$CSV2SVG" "$OUT/rc_step.csv" "$PLOTS/rc_step.svg" \
    --title "RC step response - tau = 1 ms" --series 1,2 \
    --xlabel "Time (s)" --ylabel "Voltage (V)"
"$CSV2SVG" "$OUT/ac_rlc_resonance.csv" "$PLOTS/ac_rlc.svg" \
    --title "Series RLC band-pass - f0 = 15.9 kHz, Q = 10" --xlog --series 5 \
    --xlabel "Frequency (Hz)" --ylabel "|V(out)| (V/V)"
"$CSV2SVG" "$OUT/noise_divider.csv" "$PLOTS/noise_psd.svg" \
    --title "Output noise PSD - resistive divider" --xlog --series 1 \
    --xlabel "Frequency (Hz)" --ylabel "S_out (V^2/Hz)"
"$CSV2SVG" "$OUT/hb_clipper.csv" "$PLOTS/hb_spectrum.svg" \
    --type bar --title "Harmonic balance - diode clipper (f0 = 1 kHz)" \
    --series 2,4 --ylabel "|V| (V)"
"$CSV2SVG" "$OUT/clipper_tran.csv" "$PLOTS/clipper_tran.svg" \
    --title "Diode clipper - transient" --series 1,2 \
    --xlabel "Time (s)" --ylabel "Voltage (V)"
"$CSV2SVG" "$OUT/bsim4_inverter.csv" "$PLOTS/bsim4_inverter.svg" \
    --title "CMOS inverter (BSIM4) - transient" --series 2,3 \
    --xlabel "Time (s)" --ylabel "Voltage (V)"
"$CSV2SVG" "$OUT/sens_divider.csv" "$PLOTS/sens_bar.svg" \
    --type bar --title "DC sensitivity of V(out) - divider" --series 1 \
    --ylabel "dV(out)/dR (1/ohm)"

echo
echo "Demos e figuras atualizados: demos/out e demos/plots"
