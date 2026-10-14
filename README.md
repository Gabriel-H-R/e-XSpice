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
- **SPICE Netlist Parser:** Native support for parsing `.cir` files — R/C/L, independent V/I sources (DC, PULSE, SIN, PWL, AC), diodes with `.MODEL` cards (parameter expressions like `1k`, `100n`, `1u`), `IC=` initial conditions and `+` continuation lines.
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
    operating points (`IC=`), and local step rejection with state rollback.
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
│   ├── devices/        # Physical models (R, C, L, V/I sources, diode, waveforms)
│   ├── analysis/       # DaeSystem (DAE/NR core), DCAnalysis, TransientAnalysis
│   ├── solver/         # DaeMatrix (fixed pattern), StateStore, LinearSolver
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
`storeState`). Adding a device (MOSFET, BJT, …) means implementing those — no
analysis code needs to change.

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
```

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

### MPI (Distributed Memory)

- Distributes circuit devices across processes
- Each process assembles a local portion of the global conductance matrix (G)
- Uses `MPI_Allreduce` for synchronization during iterative solving

### OpenMP (Shared Memory)

- Parallelizes matrix assembly within each MPI rank
- Accelerates Sparse Matrix-Vector multiplication (SpMV) in solvers (CG / BiCGSTAB)

### Scalability Profile

- For large circuits (>10k nodes), the hybrid parallelism footprint:

```text
4 MPI ranks × 4 OpenMP threads = 16 logical cores
```

achieves ~10–12× speedup vs serial execution


## Next Steps (Development Roadmap)
 - [x] Non-linear devices (Diodes + Newton-Raphson): PNJ + GMIN + PNJLIM + junction charge
 - [ ] Diode level-2 / sidewall junction (JSW, CJSW, VJSW)
 - [ ] MOSFET Basic Quadratic model
 - [ ] BJT Basic Exponencial model
 - [ ] AC Analysis (frequency domain) in the DAE framework
 - [ ] MOSFET (BSIM) models
 - [ ] Preconditioning (ILU)
 - [ ] Graph partitioning (METIS)
 - [ ] Distributed (MPI) assembly + solve


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