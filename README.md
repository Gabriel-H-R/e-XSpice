# e-XSpice: High-Performance Parallel Circuit Simulator

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.15%2B-success.svg)
![MPI](https://img.shields.io/badge/MPI-Distributed_Memory-orange.svg)
![OpenMP](https://img.shields.io/badge/OpenMP-Shared_Memory-yellow.svg)
![License](https://img.shields.io/badge/License-GPL_3.0-lightgrey.svg)

**e-XSpice** is an educational, yet highly scalable, electronic circuit simulator built from scratch in modern C++. Designed with High-Performance Computing (HPC) principles in mind, it utilizes a hybrid **MPI + OpenMP** parallelization strategy to accelerate the simulation of massive circuit netlists. 

This project explores the intersection of numerical analysis, sparse matrix computation, and distributed systems, providing a foundation for modern Electronic Design Automation (EDA) tools.

## Table of Contents
- [Core Features](#-core-features)
- [Mathematical Foundation (MNA)](#-mathematical-foundation-mna)
- [Software Architecture](#-software-architecture)
- [Building and Compilation](#-building-and-compilation)
- [Testing](#-testing)
- [Usage Examples](#-usage-examples)
- [Parallelization Strategy](#-parallelization-strategy)
- [Development Roadmap](#-development-roadmap)

---

## Core Features

### Numerical & Simulation Capabilities
- **SPICE Netlist Parser:** Native support for parsing `.cir` files — R/C/L, independent V/I sources (DC, PULSE, SIN, PWL and the `AC <mag> [phase]` small-signal spec), diodes with `.MODEL` cards (parameter expressions like `1k`, `100n`, `1u`), `IC=` initial conditions and `+` continuation lines. A source carries both specs independently, as in SPICE: `V1 in 0 DC 5 AC 1m 90` biases the node at 5 V and excites a 1 mV phasor at 90°, while `V1 in 0 AC 1` is a pure small-signal source (DC 0).
- **DAE formulation (Xyce-style):** the circuit is solved as the nonlinear DAE
  `C(x)·ẋ + F(x) = b(t)`. Every device contributes:
  - `F(x)` — static residual (Kirchhoff laws + constitutive laws),
  - `Q(x)` — charges/fluxes (capacitor `Q = C·V`, inductor `Φ = L·i_L`, junction charge),
  - `dF/dx` and `dQ/dx` — the two Jacobians; the **time integrator** (not the device)
    forms the effective system `J = dFdx + s·dQdx` with `s = 2/dt` (trapezoidal) or
    `1/dt` (backward Euler). No device contains companion-model or integration code.
- **Analysis types:**
  - **DC operating point (.OP)** with Newton-Raphson, GMIN homotopy stepping and
    per-junction GMIN for nonlinear devices;
  - **Transient (.TRAN)** with trapezoidal/backward-Euler integration,
    **LTE-controlled adaptive time stepping**, breakpoint clamping, initial-condition
    operating points (`IC=`), and local step rejection with state rollback;
  - **DC sweep (.DC)** in the Xyce subset — `[LIN|DEC|OCT]` over a source, a
    resistor or `TEMP`, with an optional second (inner) sweep specification:
    - one `DaeSystem` setup for the whole sweep (pattern, CSR slots, states and
      MPI plan built once) plus a warm-started continuation chain and a linear
      predictor, instead of one operating point per analysis;
    - source-step refinement (bisection) when a point fails, and
      `.OPTIONS NOFAIL` to keep sweeping past a non-convergence;
    - streamed CSV (full precision) with a decimated in-memory sample, so RAM
      stays O(nodes) on long sweeps;
    - `.OPTIONS SWEEPPAR` splits the sweep **points** across ranks (each rank
      keeps its own chain, points are exchanged in chunks). It pays off when
      the per-point solve is expensive (measured 1.5x on 4 ranks for a
      400-node diode ladder); for a handful of devices it loses to the
      sequential chain, which stays the default and bit-reproducible.
  - **Small-signal AC (.AC)** — `LIN|DEC|OCT` frequency sweep around the DC
    operating point, with per-source `AC <mag> [phase]` phasors (SPICE's dual
    spec, so one netlist carries both the bias and the excitation):
    - **one linearization for the whole sweep.** The frequency system is
      `(dF/dx + jω·dQ/dx)·X = B_ac`, and `dF/dx` / `dQ/dx` are exactly the two
      Jacobians the transient path already builds at the operating point, so
      devices need **no new physics hook** — only the excitation vector
      (`Device::loadACExcitation`). The device evaluation, the sparsity pattern
      and the symbolic factorization happen once; each frequency only refills
      the values (O(nnz)) and factorizes.
    - **native complex sparse LU** (`src/solver/ComplexSparseLU.{h,cpp}`)
      instead of Xyce's 2x2 real block form: half the factorization work and
      half the memory (4/3·N³ instead of 8/3·N³ real flops), split real/imag
      values (SoA, no `std::complex` temporaries), and a static fill-in pattern
      reused at every frequency.
    - **MNA-correct ordering:** the system is pivot-safe by construction — the
      ground node is eliminated (it is not an unknown), the node `gmin` is
      excluded (it is a DC convergence aid, not physics), and a static row
      permutation computed once from the first frequency's matrix gives every
      column a usable diagonal (MNA KCL rows of nodes wired only to
      sources/inductors have no self-term). An ideal source therefore pins its
      node **exactly** (measured: 2.000000000000000, 0°).
    - **frequencies run in parallel** — OpenMP across the sweep by default,
      and `.OPTIONS ACFREQPAR` splits them across MPI ranks in batched
      exchanges. Both are deterministic by construction: frequencies are
      independent and every rank factorizes the same replicated matrix, so a
      parallel run is **byte-identical** to the serial one (ctest checks it).
      Measured on a 401-unknown RC ladder with 501 frequencies: 4.24 s serial
      → 2.35 s with 4 OpenMP threads (1.8x) → 1.78 s with 4 ranks (2.4x,
      sub-linear because the operating point and the setup are replicated).
      On small circuits the sweep is not the bottleneck at all — for a
      401-unknown ladder the whole sweep is ~0.5 ms per frequency, so the
      operating point dominates and neither lever helps.
- **Circuit temperature (`.TEMP <celsius>`, also `.OPTIONS TEMP=`):** re-scales
  every model card that did not pin its own `TEMP=` (SPICE precedence:
  instance > card > circuit) and refreshes the devices that precompute
  temperature scalars, e.g. resistor TC1/TC2. `TEMP` is also a valid `.DC` sweep
  variable, in which case the shared cards are re-scaled per point while the
  same DAE system keeps solving.
- **Robust nonlinear devices:** PN junction diode with the full SPICE3/Xyce level-1 model
  (Shockley + high-injection knee + reverse third-power law + breakdown + transit-time
  diffusion charge + C⁰-continuous junction capacitance + PNJLIM limiting + GMIN).
- **Fixed-sparsity-pattern assembly:** devices declare their Jacobian stamp once and
  write into pre-resolved CSR slots — no `std::map`/sorting per iteration, with
  device-parallel assembly (OpenMP) above a size threshold.
- **Solvers:** Jacobi-preconditioned BiCGSTAB (scalable) + dense-LU direct fallback
  (robustness), with automatic switch for small/medium systems.
- **HPC:** hybrid MPI + OpenMP architecture with parallel-size thresholds, so small
  circuits never pay parallel-region overhead.

---

## Mathematical Foundation (MNA)

e-XSpice relies on Modified Nodal Analysis (MNA) to formulate the circuit equations. The simulator constructs and solves the following sparse linear system:

$$\mathbf{G} \mathbf{x} = \mathbf{b}$$

Where:
* $$\mathbf{G}$$ is the Conductance/Incidence matrix.
* $$\mathbf{x}$$ is the vector of unknown node voltages and branch currents of voltage-defined devices.
* $$\mathbf{b}$$ is the vector of known independent source values.

For Transient Analysis, the companion model approach is used, transforming differential equations of dynamic components ($C$, $L$) into equivalent discrete-time resistive circuits for each time step.

---

## Software Architecture

```text
e-XSpice/
├── src/
│   ├── core/           # Circuit topology, Node, Device (DAE contract), parser
│   ├── devices/        # Device instances (R, C, L, V/I sources, diode, BJT, waveforms)
│   ├── models/
│   │   ├── core/       # Model SYSTEM: DeviceModel, ModelCard, ModelRegistry
│   │   └── core_models/# Model MATH: Diode, Bjt, passives, shared Junction physics
│   ├── analysis/       # DaeSystem (DAE/NR core), DCAnalysis, DCSweep,
│   │                   #   TransientAnalysis, ACAnalysis (small-signal)
│   ├── solver/         # DaeMatrix (fixed pattern), StateStore, LinearSolver,
│   │                   #   ComplexSparseLU (AC factorization)
│   ├── parallel/       # MPI communicator wrappers and OpenMP pragmas
│   └── main.cpp        # CLI entry point
├── tests/              # CTest suite, incl. analytical physics verification
├── examples/           # Benchmark circuits and .cir netlists
└── CMakeLists.txt      # Main CMake configuration
---
```

### The DAE device contract

```
                 ┌──────────────────────────────┐
   netlist ─────▶│  Parser → Circuit → Devices  │
                 └──────────────┬───────────────┘
                                │  setup(): pattern + indices + states
                 ┌──────────────▼───────────────┐
                 │        DaeSystem             │
                 │  F, Q, B, dFdx, dQdx load    │
                 │  J = dFdx + s·dQdx           │
                 │  Newton + gmin homotopy      │
                 │  LTE step control (transient)│
                 └──────────────────────────────┘
```

Every device implements five setup methods (`numExtraUnknowns`,
`registerJacobianPattern`, `setJacobianSlots`, `registerStates`, `setTopology`) and
six evaluation methods (`loadF`, `loadQ`, `loadB`, `loadJacobianF`, `loadJacobianQ`,
`storeState`). Adding a device (JFET, IGBT, …) means implementing those — no
analysis code needs to change.

### The AC path: one linearization, many frequencies

```text
  .AC DEC 10 10 1Meg
      │
      ▼
  DaeSystem.solveNonlinear()      ← ONE operating point
      │
      ▼
  refreshLinearization()          ← dF/dx and dQ/dx AT that point (not the
      │                            previous Newton iterate)
      ▼
  (G, C) = (dF/dx, dQ/dx)         ← one snapshot, both on the same CSR pattern
      │
      ├── ground row/column eliminated  (x(gnd) = 0 is not an unknown)
      ├── node gmin excluded            (DC homotopy aid, not physics)
      ├── pivot-safe row permutation    (computed once, from the first ω)
      ▼
  ComplexSparseLU:  symbolic once ──▶ per frequency:  factorize(G + jωC) + solve
      │
      ▼
  |X|, arg(X) per node  +  |I|, arg(I) per branch unknown
```

Because `dF/dx` and `dQ/dx` are the *transient* Jacobians, a device supports AC
the moment it supports transient: the only addition to the device contract is
`loadACExcitation` (the phasor on the right-hand side, split into real and
imaginary vectors so the device layer stays free of `<complex>`). A resistor,
a MOSFET and a BSIM4 device need no AC code at all.

### The model layer (`src/models/`)

Following Xyce's `DeviceModelPKG` — but with the temperature math hoisted out
of the per-instance path — every device is agnostic to its model contents:

```
.MODEL DMOD D (IS=1e-14 N=1.0 ...)     netlist card
        │  generic ModelCard (no physics in the parser)
        ▼
ModelRegistry ──factory──▶ DiodeModel ──shared──▶ D1, D2, ... (RCP)
   (D/R/C/L/NPN/PNP/      (validates once, scales to TEMP once;
    NMOS/PMOS built-in;     instances read via pointer, zero loop cost)
    new type = 1 line)
```

- `DeviceModel`: abstract card (name, type letter, `processParams`, `updateTemperature`) **plus the declarative parameter table** (`declareParams` + `applyCardParams`, Xyce-`addPar`-style): each subclass lists its keys once — `makeParam` (canonical), `makeAlias` (SPICE spelling, canonical wins on conflict), `makeCelsiusParam` (`TNOM`/`TEMP` arrive °C, stored K), `makeIntParam` (`LEVEL`) — with one-line help each.
- `ModelRegistry`: owns cards by name; factories only construct defaults (one-line lambdas, never parse); the registry fills card values **generically** from the table, rejects unknown keys with the valid list (`... unknown parameter(s) 'BTEA' (valid parameters: IS N RS ...)`), and lists registered types on unknown `.MODEL` types. Typed checked lookup (`findAs<T>`).
- Adding a model = subclass + table + `processParams()` + one factory line. No parser change, no `fromCard` loop, no forgotten alias, no silent typo.
- Instance lines use the SAME machinery: each `*InstanceParams` struct
  declares its own option table (`AREA`/`M`/`W`/`L`/`TEMP`… with the
  `TEMPGiven` hook), and one generic `parseInstanceTail` fills it —
  unknown options fail with the valid list (a misspelled `AREA=` used to
  be silently ignored). New instance options never touch parser code.
- `DiodeModel`: full level-1 physics + one `DiodeThermal` snapshot per card (bandgap `VJ(T)`/`CJO(T)`, `IS(T)`, `tVcrit`, breakdown fit, `F2/F3/tF1` continuity constants). Instances keep only `AREA`/`M`/`TEMP` overrides — an instance-level `TEMP` gets its own snapshot, otherwise it reads the card's. `BV`/`VB` set `BVGiven` via the `onParamGiven` hook.
- `BjtModel`: **level 1 (Ebers-Moll transport) + level 2 (Gummel-Poon), NPN *and* PNP from one polarity-agnostic card** — `LEVEL`, `IS/BF/BR/NF/NR`, base-width modulation (`VAF`/`VAR`), transit time `TF`/`TR`, per-junction depletion capacitance (`CJE/VJE/MJE`, `CJC/VJC/MJC`), `IS(T)` per junction, `gmin`; level 2 adds the normalized base charge `qb` (high injection via `IKF/IKR`), non-ideal base currents (`ISE/NE`, `ISC/NC`), ohmic `RB/RC/RE` (internal nodes, zero-cost when absent), and the substrate diode + `CJS/VJS/MJS` (optional 4th terminal). With all GP extras off, level 2 reduces to level 1 (`qb → q1 = 1/K_E`). See `Bjt.h` for the complete derivation and stamping.
- `MosfetModel`: **level 1 (Shichman-Hodges), NMOS *and* PMOS from one polarity-agnostic card** — `VTO/KP/LAMBDA/PHI/GAMMA` (VTO follows the SPICE sign, negative for PMOS), bulk `IS/PB/MJ/FC` + `CBD/CBS`, overlap `CGSO/CGDO/CGBO`, `TOX/LD` (effective `BETA = KP·W/Leff`, `Cox = εox/TOX·W·Leff`), `IS(T)/PHI(T)/PB(T)`. Instances carry `W/L/M/AD/AS/TEMP`.
- `PassiveModels`: `ResistorModel` (`R` multiplier + `TC1`/`TC2`, wired: `R1 1 2 RMOD 1k`), `CapacitorModel`/`InductorModel` ready for the same path.
- `Junction.h`: the PN-junction math **shared by diode, BJT and MOSFET bulk** — `Q(V)`/`C(V)` C⁰-continuous at `FC·VJ`, the bandgap temperature scaling and `IS(T)`. Header-only and `inline`: every device inlines it into `loadF`/`loadQ`, one source of truth, zero call overhead.
- Parser dispatch is data-driven: device letters go through a creator map, model types through the registry — a new device needs no parser-core changes.

### The BJTs (`src/devices/Bjt.{h,cpp}`)

```
Q<name>  <collector> <base> <emitter> [substrate]  <NPN|PNP card>
         [AREA=<a>] [TEMP=<degC>]
```

* Level 1: three terminals, no internal nodes; one code path for both
  polarities (the device remaps the junction voltages, the card stays
  polarity-agnostic, and **every intrinsic Jacobian entry is identical for
  NPN and PNP** — verified by the mirror tests, which compare all terminal
  values bit-for-bit).
* Terminal currents: `i_b = I_BE + I_BC`, `i_c = I_CT − I_BC`, `i_e = −i_b − i_c`;
  charges: `Q_BE = dep(vbe) + τF·i_c`, `Q_BC = dep(vbc) − τR·i_c` — including
  the transit-charge cross derivatives, so Newton keeps quadratic convergence.
* Level 2 (`LEVEL=2`): Gummel-Poon — transport current divided by `qb`,
  non-ideal base currents, `RB/RC/RE` via internal unknowns (same
  extra-unknown range mechanism as the diode's `Pri` node; zero unknowns
  when the resistances are zero), optional substrate terminal with its own
  diode + charge + state. Frame rule for polarity: intrinsic laws are
  N-oriented (p² = 1, polarity-free Jacobian), ohmic drops stay physical
  (no `p` anywhere) — mixing the frames negates the PNP equation.
* Robustness: SPICE3f5 PNJLIM on **both** junctions, and the limiting
  write-back (`applySolutionLimit`) that restarts the Newton iteration at the
  limited point — without it a current-source-driven junction (no resistive
  path) never converges. The diode uses the same hook.

### The MOSFETs (`src/devices/Mosfet.{h,cpp}`)

```
M<name>  <drain> <gate> <source> <bulk>  <NMOS|PMOS card>
         [W=<w>] [L=<l>] [M=<m>] [AD=<a>] [AS=<a>] [TEMP=<degC>]
```

* Level 1 Shichman-Hodges, one polarity-agnostic code path (same `p`
  remapping trick as the BJT; `VTO` keeps the SPICE sign, `GAMMA/PHI`
  stay positive — reverse bulk bias is negative in the N-oriented frame
  for both). Source/drain swap internally when `vds < 0`, mapping
  currents, derivatives and Meyer caps back to the physical frame, so a
  single unified 4×4 stamp serves all quadrants.
* Channel `Id/gm/gds/gmbs` with exact region derivatives, always-on bulk
  diodes (`IS + gmin`), Meyer partition (`Cox/2` linear, `2Cox/3`
  saturation, `Cox` cutoff) + overlap + bulk depletion in one
  charge-conserving stamp (`Qs = −Qg−Qd−Qb`, three states).
* Limiting: PNJLIM on both bulk diodes + bounded-step clamp on the channel
  (quadratic needs no exponential limiter), gate-anchored write-back
  (never fights a voltage-driven gate).

### Newton limiting (write-back)

SPICE/Xyce semantics, now shared by every nonlinear device: an iteration that
a device rejects (`acceptSolution() == false`) is followed by the device
writing its **limited unknowns back into the solution vector**, so the next
iteration starts from the limited point instead of walking the exponential up
by `Vt` per iteration. Serial path only (the distributed path must keep `x`
identical on every rank). Two refinements over the textbook hook: the
write-back shifts external + internal terminals together (preserving every
ohmic drop — otherwise it fabricates kiloamp residuals into nearly singular
base blocks), and each device re-anchors its limiter history from raw values
(`resetIterationState`) whenever the solver restores a sane point, so a
poisoned history can never misfire `pnjlim` on a retry.

### Newton core hardening (`src/analysis/DaeSystem.{h,cpp}`)

* **Direct LU for small systems** (≤ 2000 rows): Krylov on a nearly singular
  Newton Jacobian can "converge" with arbitrary nullspace components while LU
  returns the deterministic answer — and is cheaper than up-to-2000 Krylov
  iterations per Newton step. Krylov stays the workhorse above that size.
* **Residual convergence test**: `|Δx|`-small no longer suffices — with a
  singular Jacobian it coexists with a big residual. Acceptance now also
  requires the freshly reassembled residual small against the equation scale
  (one extra assembly per solve), plus a re-check that no limiter fired.
* **Adaptive trust region** (node voltages only; branch currents are
  Lagrange multipliers and must stay free): grows ×2 on residual progress,
  shrinks ×0.5 on growth, between 1e-3 V and 100 V. Lets kilovolt circuits
  step freely while herding floating-node excursions down.
* **GMIN homotopy from the sane point**: a failed direct Newton restores the
  pre-attempt point (not its own diverged garbage) before stepping `gmin`
  down — homotopy then tracks the physical branch instead of a parasitic one.
* **DC warm start**: `DCAnalysis::run(use_initial_guess=true)` now honors
  its flag (previously ignored), seeding Newton with the previous solution
  for honest continuation sweeps past bistabilities.

## Building and Compilation 

The project uses a modern **CMake workflow**, ensuring a seamless and unified build process across **Windows, Linux, and macOS**.

### 1. Prerequisites

- C++17 compatible compiler (**MSVC, GCC, or MinGW**)
- **CMake 3.15+**
- **MPI Implementation:**
  - **Windows:** MS-MPI (install both `msmpisetup.exe` and `msmpisdk.msi`)
  - **Linux:** OpenMPI or MPICH
- **OpenMP** (usually bundled with modern C++ compilers)

---

### Cross-Platform Compilation

Open your terminal (or Developer Command Prompt / PowerShell on Windows) and run:

### Clone the repository
```text
git clone [https://github.com/Gabriel-H-R/e-XSpice.git](https://github.com/Gabriel-H-R/e-XSpice.git)
cd e-XSpice
```

### Compilation Details

 Windows (Using MinGW / Git Bash)

Ensure your MinGW compiler and MS-MPI paths are properly configured in your system's `PATH`

```text
# Generate build files for MinGW
cmake -S . -B build -G "MinGW Makefiles"

# Compile using multiple CPU cores
cmake --build build --config Release -j
```

 Linux (Ubuntu / Debian)

Install the required dependencies:


```text
sudo apt-get update
sudo apt-get install cmake build-essential libopenmpi-dev openmpi-bin
```

Compile the project:

```text
# Generate standard Unix Makefiles
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Compile using all available cores
cmake --build build -j $(nproc)
```
---

## Testing

e-XSpice includes a suite of automated unit tests to validate the physical accuracy of the solvers against expected SPICE results. To run the tests using **CTest**:

```text
cd build
ctest --output-on-failure
```

## Usage Examples

In the `examples/` folder, you will find some benchmark circuits and `.cir` netlists. After compilation, the executable will be located in the `build/` directory (`e-XSpice` on Linux/macOS, `e-XSpice.exe` on Windows).

### Mode 1: Simulating Custom Netlists (.cir)

You can run a custom SPICE netlist by passing the `file` argument:

```text
# Running an RC Step Response (Transient Analysis)
./build/e-XSpice file ./examples/rc_step.cir

# Running a Multi-loop Resistor Network (DC Analysis)
./build/e-XSpice file ./examples/multiloop.cir

# Small-signal frequency response of an RC low-pass (AC analysis)
./build/e-XSpice file ./tests/ac_rc_netlist.cir out.csv
```

An AC netlist only needs the excitation and the sweep:

```spice
* RC low-pass: |H| = 1/sqrt(1+(wRC)^2), phase = -atan(wRC)
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 1n
.AC DEC 10 10 1Meg
.END
```

```text
$ ./build/e-XSpice file ac_rc.cir ac_rc.csv
[ACAnalysis] DEC sweep 51 frequencies, 10 Hz .. 1e+06 Hz
[ACAnalysis] operating point in 1 Newton iterations
AC sweep: 51/51 frequencies (last 1e+06 Hz)
AC analysis completed: 51 frequencies (DEC), 3 unknowns
  OP iterations: 1 | factorizations: 51 | LU flops: 485316 | ranks: 1 | time: 0.0006 s

Magnitude / phase:
  f=10 Hz: |Vin|=1V@0deg |Vout|=1V@-0.0036deg |I|=6.28e-08A@-90deg
  ...
  f=100000 Hz: |Vin|=1V@0deg |Vout|=0.846733V@-32.1419deg |I|=0.000532A@-122.14deg
```

The CSV carries `frequency`, then `<node>_mag` / `<node>_phase` per node and
`<device>_imag_mag` / `<device>_imag_phase` per branch unknown (`I(V1)`, `I(L1)`).

### Mode 2: Programmatic Examples & MPI

The simulator also contains hardcoded examples for scalability testing.

**Serial Execution**

```text
./build/e-XSpice 1
```

**Distributed MPI Execution (setting 4 processes and 1k resistors)**

```text
mpiexec -np 4 ./build/e-XSpice 3 1000
```

Note:

On Linux, you may use `mpirun` instead of `mpiexec`
On Windows (MS-MPI), always use `mpiexec`

### OpenMP Thread Control

You can control the shared-memory parallelism by setting the `OMP_NUM_THREADS` environment variable before execution.

Windows (PowerShell)

```text
$env:OMP_NUM_THREADS="4"
mpiexec -np 2 ./build/e-xspice 3 10000
```

Linux (Bash)

```text
export OMP_NUM_THREADS=4
mpiexec -np 2 ./build/e-xspice 3 10000
```

---

## Parallelization Strategy

The simulator runs in **SPMD mode**: every rank builds the same circuit and
holds row-replicated assembly buffers plus the fixed sparsity pattern. Work
— not data — is partitioned (see `src/parallel/Distribution.h`):

### MPI (Distributed Memory)

- **Device-partitioned assembly:** each rank evaluates a disjoint device
  subset (greedy longest-processing-time assignment weighted by device cost:
  diodes/modal devices cost more than passives) into local `F/Q/B/J`
  buffers, fused by one SUM `Allreduce` per buffer per Newton iteration
  (`F/Q/B` travel packed in a single round trip). Active above 64 devices
  by default (configurable via `DaeSystem::setDistributionThresholds`).
- **Row-distributed Krylov solve:** each rank owns a contiguous row block;
  the solution stays replicated (Newton updates are identical everywhere),
  so a matvec needs only one `Allgatherv` of the locally computed block
  plus `Allreduce`d dots — no ghost exchange, no distributed matrix.
  Active above 20000 unknowns by default (configurable); below that the
  replicated Krylov + direct-LU fallback runs redundantly (cheaper than
  the collectives).
- **Identical control flow:** Newton/LTE/dt decisions derive only from
  globally synchronized quantities (limiter votes use one `AND`), so every
  rank takes the same branches with no extra coordination. Verified by the
  `MPI_Ladder_*` equivalence tests (serial vs 2/4-rank bit-level agreement).

### OpenMP (Shared Memory)

- Parallelizes device assembly and SpMV within each MPI rank (serial fast
  paths below size thresholds, so small circuits never pay region overhead).

### Tuning

```text
# Keep total threads ≈ core count (ranks × threads ≤ cores):
mpiexec -n 4 ./build/e-XSpice file examples/clipper_2000.cir   # + OMP_NUM_THREADS=2

# Force distribution on/off regardless of thresholds (tests, debugging):
EXSPICE_DISTRIBUTE=1  # force on     EXSPICE_DISTRIBUTE=0  # force serial
```

Rule of thumb: MPI distribution pays off for device-heavy circuits (costly
model evaluations amortize the ~5 collectives per Newton iteration) and for
very large systems via the distributed solve; on a single shared-memory box
with small/medium linear circuits, serial or pure OpenMP is often faster —
the thresholds encode exactly that trade-off. Demo circuits:
`examples/ladder_200.cir` (linear), `examples/clipper_500.cir` and
`examples/clipper_2000.cir` (nonlinear diode ladders).

For `.AC` the two levers are different in kind, and they multiply into
oversubscription if you use both blindly:

```text
# Default: every rank solves the whole sweep, OMP threads split the FREQUENCIES
OMP_NUM_THREADS=4 ./build/e-XSpice file ac.cir          # 1.8x on 4 threads

# `.OPTIONS ACFREQPAR`: ranks split the frequencies; keep OMP=1 to avoid
# oversubscription (measured 4 ranks x 4 threads = slower than either alone)
mpiexec -n 4 ./build/e-XSpice file ac_freqpar.cir       # 2.4x on 4 ranks
```

Both are exact: the CSV of a parallel run is byte-identical to the serial one,
which the CTest suite checks (`AC_Netlist_MPI`, `AC_Netlist_FreqPar`).


## Next Steps (Development Roadmap)
 - [x] Non-linear devices (Diodes + Newton-Raphson): PNJ + GMIN + PNJLIM + junction charge
 - [x] DeviceModel package (registry + factories, shared D/R/C/L cards, instance overrides)
 - [x] **Transistors NPN + PNP level 1** (Ebers-Moll transport: Early + transit time + per-junction capacitance, one polarity-agnostic code path)
 - [x] **Transistors NPN + PNP level 2** (Gummel-Poon: qb high injection, ISE/ISC non-ideal base, RB/RC/RE internal nodes, substrate terminal)
 - [x] **MOSFETs NMOS + PMOS level 1** (Shichman-Hodges + Meyer charges, symmetric vds<0 handling, one polarity-agnostic code path)
 - [x] SPICE-style limiting write-back (`applySolutionLimit`) for robust current-source-driven junctions
 - [x] Newton core hardening (direct LU ≤2000 rows, residual convergence test, adaptive trust region, sane-point homotopy restart, limiter-history reset, DC warm start)
 - [ ] Diode level-2 / sidewall junction (JSW, CJSW, VJSW)
 - [ ] BJT level-3 extras (base-resistance modulation RBM/IRB, excess phase PTF for AC, avalanche BV/IBV)
 - [ ] MOSFET level-2/3 (short-channel VMAX/THETA, subthreshold NFS, RD/RS/RG ohmic)
  - [x] **AC analysis (frequency domain) in the DAE framework** - one operating point
        and one linearization per sweep, native complex sparse LU with a static
        fill-in pattern, MNA-correct pivot-safe ordering, frequencies parallel
        over OpenMP/MPI with byte-identical results
  - [ ] AC extras (`.PRINT`/`.PLOT` column selection, `.NOISE`)
  - [x] MOSFET BSIM4 model (level 14, Berkeley-derived level-2 core)
  - [ ] Preconditioning (ILU)
  - [ ] Graph partitioning (METIS) - and, for the AC path, a row-distributed
        complex LU with coordinated pivoting for million-unknown systems
 - [x] Distributed (MPI) assembly + solve (device-partitioned Allreduce fusion + row-distributed BiCGSTAB)


### License

This project is licensed under the GPL-3.0 License.
See the ```LICENSE``` file for details.

## References

### Core Resources
- [SPICE (UC Berkeley)](https://people.eecs.berkeley.edu/~newton/spice/)
- [ngspice](https://ngspice.sourceforge.io)
- [Xyce Parallel Electronic Simulator](https://xyce.sandia.gov/)

### Books (Circuit Simulation & SPICE)
- Vladimirescu, A. — *The SPICE Book* (McGraw-Hill, 1994)  
- Najm, F. N. — *Circuit Simulation* (Wiley-IEEE Press, 2010) 
- Tuma, T., Burmen, Á. — *Circuit Simulation with SPICE OPUS: Theory and Practice* (Birkhäuser, 2009)

### Books (Numerical Methods & Advanced Simulation)
- Benner, P., Hinze, M., ter Maten, E. — *Model Reduction for Circuit Simulation* (Springer, 2011) 

### High Performance Computing (HPC) 
- Hager, G., Wellein, G. — *Introduction to High Performance Computing for Scientists and Engineers* (CRC Press, 2010)
