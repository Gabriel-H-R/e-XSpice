<#
.SYNOPSIS
  Runs every showcase netlist through the CLI and rebuilds the README figures.

.DESCRIPTION
  Writes CSVs and console logs to demos/out and regenerates the SVG plots in
  demos/plots with the csv2svg tool. Safe to re-run at any time.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File demos/run_demos.ps1
#>

param(
    [string]$Exe = "$PSScriptRoot/../build/e-XSpice.exe",
    [string]$Csv2Svg = "$PSScriptRoot/../build/csv2svg.exe"
)

$ErrorActionPreference = 'Stop'
$out = "$PSScriptRoot/out"
$plots = "$PSScriptRoot/plots"
New-Item -ItemType Directory -Force $out, $plots | Out-Null

$demos = @(
    'rc_step',
    'ac_rlc_resonance',
    'noise_divider',
    'sens_divider',
    'hb_clipper',
    'clipper_tran',
    'bsim4_inverter'
)

foreach ($d in $demos) {
    & $Exe file "$PSScriptRoot/$d.cir" "$out\$d.csv" *> "$out\$d.log"
    if ($LASTEXITCODE -ne 0) { throw "demo '$d' failed (see $out\$d.log)" }
    Write-Host "  $d ok"
}

& $Csv2Svg "$out\rc_step.csv" "$plots\rc_step.svg" `
    --title "RC step response - tau = 1 ms" --series 1,2 `
    --xlabel "Time (s)" --ylabel "Voltage (V)"
& $Csv2Svg "$out\ac_rlc_resonance.csv" "$plots\ac_rlc.svg" `
    --title "Series RLC band-pass - f0 = 15.9 kHz, Q = 10" --xlog --series 5 `
    --xlabel "Frequency (Hz)" --ylabel "|V(out)| (V/V)"
& $Csv2Svg "$out\noise_divider.csv" "$plots\noise_psd.svg" `
    --title "Output noise PSD - resistive divider" --xlog --series 1 `
    --xlabel "Frequency (Hz)" --ylabel "S_out (V^2/Hz)"
& $Csv2Svg "$out\hb_clipper.csv" "$plots\hb_spectrum.svg" `
    --type bar --title "Harmonic balance - diode clipper (f0 = 1 kHz)" `
    --series 2,4 --ylabel "|V| (V)"
& $Csv2Svg "$out\clipper_tran.csv" "$plots\clipper_tran.svg" `
    --title "Diode clipper - transient" --series 1,2 `
    --xlabel "Time (s)" --ylabel "Voltage (V)"
& $Csv2Svg "$out\bsim4_inverter.csv" "$plots\bsim4_inverter.svg" `
    --title "CMOS inverter (BSIM4) - transient" --series 2,3 `
    --xlabel "Time (s)" --ylabel "Voltage (V)"
& $Csv2Svg "$out\sens_divider.csv" "$plots\sens_bar.svg" `
    --type bar --title "DC sensitivity of V(out) - divider" --series 1 `
    --ylabel "dV(out)/dR (1/ohm)"

Write-Host ''
Write-Host 'Demos e figuras atualizados: demos/out e demos/plots'
