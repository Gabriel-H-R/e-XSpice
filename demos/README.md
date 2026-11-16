# Demos & Showcase

Runnable circuits that exercise the features highlighted in the main README.
Every demo is also an end-to-end **CTest** case (`ctest -R Demo_`), so the
showcase is checked in CI, not just documented.

## Run everything

```bash
./demos/run_demos.sh            # Linux / macOS (Git Bash works too)
```

```powershell
powershell -ExecutionPolicy Bypass -File demos/run_demos.ps1   # Windows
```

Both scripts run every netlist through the CLI, save CSV + console log to
`demos/out/` and rebuild the SVG figures in `demos/plots/` with the
`csv2svg` tool (built by CMake next to the simulator).

## The circuits

| Netlist | Analysis | What it shows |
|---|---|---|
| `rc_step.cir` | `.TRAN` | RC step response against the closed form `5·(1−e^(−t/RC))`, τ = 1 ms |
| `ac_rlc_resonance.cir` | `.AC` | Series RLC band-pass: f₀ = 15.9 kHz, Q = 10, worst \|H\| error 1.5×10⁻¹⁵ over 101 frequencies |
| `noise_divider.cir` | `.NOISE` | Resistive divider: S_out = 8.288035947×10⁻¹⁷ V²/Hz, matching 4kT/R theory to machine precision, plus the per-generator budget |
| `sens_divider.cir` | `.SENS` | DC sensitivity dV(out)/dR = ±2.5×10⁻⁴ against the exact closed form |
| `hb_clipper.cir` | `.HB` | Harmonic balance of a diode clipper: fundamental 0.488 V, 3rd harmonic 9.2 mV (1.9 %) |
| `clipper_tran.cir` | `.TRAN` | Nonlinear diode clipper waveform (clips at ±0.574 V) |
| `bsim4_inverter.cir` | `.TRAN` | CMOS inverter with Berkeley BSIM4 MOSFETs, rail-to-rail switching |
| `ac_ladder_200_freqpar.cir` | `.AC` | 5001-frequency sweep on a 202-unknown ladder — the parallel scaling benchmark |

## Parallel benchmark

```powershell
powershell -ExecutionPolicy Bypass -File demos/benchmark_scaling.ps1
```

Runs the dense AC ladder with 1/2/4/8 OpenMP threads and 1/2/4 MPI ranks,
records the simulation time and SHA-256 of each CSV in `demos/out/`, prints
the speedups and rebuilds `demos/plots/speedup.svg`.

Representative result on an Intel i5-9300HF (4 cores), CSV written to disk:

| Config | Time | Speedup |
|---|---|---|
| 1 thread | 6.46 s | 1.00× |
| 2 threads | 5.92 s | 1.09× |
| 4 threads | 6.09 s | 1.06× |
| 8 threads | 6.37 s | 1.01× |
| MPI 2 ranks | 5.46 s | 1.18× |
| MPI 4 ranks | 5.08 s | 1.27× |

Timings move with the CPU turbo/thermal state; the SHA-256 of every CSV is
printed so you can check the main claim directly: the distributed runs are
**byte-identical** to the serial one. That is also what the
`AC_Netlist_MPI` / `AC_Netlist_FreqPar` CTest cases pin for the reference
netlists.

## Plotting your own results

`csv2svg` is dependency-free (standard library only) and reads any simulator
CSV:

```bash
./build/csv2svg results.csv figure.svg \
    --title "My sweep" --xlog --series 2,4 \
    --xlabel "Frequency (Hz)" --ylabel "Magnitude"
```

Options: `--type line|bar`, `--xlog`, `--ymin/--ymax`, `--title`, `--xlabel`,
`--ylabel`. Bar charts accept a text category column, so they also work for
sensitivity/spectrum data.
