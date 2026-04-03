# ⚡ e-XSpice: High-Performance Parallel Circuit Simulator

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.15%2B-success.svg)
![MPI](https://img.shields.io/badge/MPI-Distributed_Memory-orange.svg)
![OpenMP](https://img.shields.io/badge/OpenMP-Shared_Memory-yellow.svg)
![License](https://img.shields.io/badge/License-GPL_3.0-lightgrey.svg)

**e-XSpice** is an educational, yet highly scalable, electronic circuit simulator built from scratch in modern C++. Designed with High-Performance Computing (HPC) principles in mind, it utilizes a hybrid **MPI + OpenMP** parallelization strategy to accelerate the simulation of massive circuit netlists. 

This project explores the intersection of numerical analysis, sparse matrix computation, and distributed systems, providing a foundation for modern Electronic Design Automation (EDA) tools.

## 📋 Table of Contents
- [Core Features](#-core-features)
- [Mathematical Foundation (MNA)](#-mathematical-foundation-mna)
- [Software Architecture](#-software-architecture)
- [Building and Compilation](#-building-and-compilation)
- [Testing](#-testing)
- [Usage Examples](#-usage-examples)
- [Parallelization Strategy](#-parallelization-strategy)
- [Development Roadmap](#-development-roadmap)

---

## ✨ Core Features

### 🧮 Numerical & Simulation Capabilities
- **SPICE Netlist Parser:** Native support for parsing `.cir` and standard SPICE netlist files.
- **Analysis Types:** - **DC Operating Point Analysis (.OP):** Steady-state resolution of linear and multi-loop networks.
  - **Transient Analysis (.TRAN):** Time-domain simulation using Numerical Integration (backward Euler method) for time-dependent devices (Capacitors and Inductors).
- **Formulation:** Modified Nodal Analysis (MNA) for robust linear system modeling.
- **Iterative Solvers:** Custom implementations of Conjugate Gradient (CG) and BiCGSTAB, optimized for sparse matrices.
- **Supported Devices:** Ideal Resistors, Capacitors, Inductors, Independent Voltage Sources, and Independent Current Sources.

### 🚀 High-Performance Computing (HPC)
- **Hybrid Architecture:** Distributed memory scaling across nodes via **MPI**, coupled with shared-memory thread parallelism within nodes via **OpenMP**.
- **Vectorization-Ready:** Inner loops structured to facilitate compiler auto-vectorization.

---

## 📐 Mathematical Foundation (MNA)

e-XSpice relies on Modified Nodal Analysis (MNA) to formulate the circuit equations. The simulator constructs and solves the following sparse linear system:

$$\mathbf{G} \mathbf{x} = \mathbf{b}$$

Where:
* $$\mathbf{G}$$ is the Conductance/Incidence matrix.
* $$\mathbf{x}$$ is the vector of unknown node voltages and branch currents of voltage-defined devices.
* $$\mathbf{b}$$ is the vector of known independent source values.

For Transient Analysis, the companion model approach is used, transforming differential equations of dynamic components ($C$, $L$) into equivalent discrete-time resistive circuits for each time step.

---

## 🏗️ Software Architecture

```text
e-XSpice/
├── src/
│   ├── core/           # Circuit topology, Node management, SPICE Parser
│   ├── devices/        # Physical models (Resistor, Capacitor, Inductor, etc.)
│   ├── analysis/       # Simulation engines (DCAnalysis, TransientAnalysis)
│   ├── solver/         # Numerical linear algebra (BiCGSTAB, CG)
│   ├── parallel/       # MPI communicator wrappers and OpenMP pragmas
│   └── main.cpp        # CLI entry point
├── tests/              # Independent CTest-based unit tests
├── examples/           # Benchmark circuits and .cir netlists
└── CMakeLists.txt      # Main CMake configuration
---
```

## 🛠️ Building and Compilation 

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

💻 Windows (Using MinGW / Git Bash)

Ensure your MinGW compiler and MS-MPI paths are properly configured in your system's `PATH`

```text
# Generate build files for MinGW
cmake -S . -B build -G "MinGW Makefiles"

# Compile using multiple CPU cores
cmake --build build --config Release -j
```

🐧 Linux (Ubuntu / Debian)

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

## 🧪 Testing

e-XSpice includes a suite of automated unit tests to validate the physical accuracy of the solvers against expected SPICE results. To run the tests using **CTest**:

```text
cd build
ctest --output-on-failure
```

## 💻 Usage Examples

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

## ⚡ Parallelization Strategy

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


## 🗺️ Next Steps (Development Roadmap)
 - [ ] Non-linear devices (Diodes + Newton-Raphson)
 - [ ] MOSFET Basic Quadratic model
 - [ ] BJT Basic Exponencial model
 - [ ] AC Analysis (frequency domain)
 - [ ] MOSFET (BSIM) models
 - [ ] Preconditioning (ILU)
 - [ ] Graph partitioning (METIS)


### 📜 License

This project is licensed under the GPL-3.0 License.
See the ```LICENSE``` file for details.

## 📚 References

### Core Resources
- [SPICE (UC Berkeley)](https://people.eecs.berkeley.edu/~newton/spice/)
- [ngspice] (https://ngspice.sourceforge.io)
- [Xyce Parallel Electronic Simulator](https://xyce.sandia.gov/)

### Books (Circuit Simulation & SPICE)
- Vladimirescu, A. — *The SPICE Book* (McGraw-Hill, 1994)  
- Najm, F. N. — *Circuit Simulation* (Wiley-IEEE Press, 2010) 
- Tuma, T., Burmen, Á. — *Circuit Simulation with SPICE OPUS: Theory and Practice* (Birkhäuser, 2009)

### Books (Numerical Methods & Advanced Simulation)
- Benner, P., Hinze, M., ter Maten, E. — *Model Reduction for Circuit Simulation* (Springer, 2011) 

### High Performance Computing (HPC) 
- Hager, G., Wellein, G. — *Introduction to High Performance Computing for Scientists and Engineers* (CRC Press, 2010)