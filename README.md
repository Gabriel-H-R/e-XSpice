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
- [Usage Examples](#-usage-examples)
- [Parallelization Strategy](#-parallelization-strategy)
- [Development Roadmap](#-development-roadmap)

---

## ✨ Core Features

### 🧮 Numerical & Simulation Capabilities
- **Analysis Types:** DC Operating Point Analysis (Transient and AC in development).
- **Formulation:** Modified Nodal Analysis (MNA) for robust linear system modeling.
- **Iterative Solvers:** Custom implementations of Conjugate Gradient (CG) and BiCGSTAB, optimized for sparse matrices.
- **Matrix Formats:** Compressed Sparse Row (CSR) and Coordinate List (COO) for memory-efficient matrix assembly.
- **Supported Devices:** Ideal Resistors, Independent Voltage Sources, Independent Current Sources.

### 🚀 High-Performance Computing (HPC)
- **Hybrid Architecture:** Distributed memory scaling across nodes via **MPI**, coupled with shared-memory thread parallelism within nodes via **OpenMP**.
- **Vectorization-Ready:** Inner loops structured to facilitate compiler auto-vectorization.

---

## 📐 Mathematical Foundation (MNA)

e-XSpice relies on Modified Nodal Analysis (MNA) to formulate the circuit equations. The simulator constructs and solves the following sparse linear system:

$$\mathbf{G} \mathbf{x} = \mathbf{b}$$

Where:
* $$\mathbf{G}$$ is the Conductance/Incidence matrix (sparse, symmetric for pure linear resistor circuits).
* $$\mathbf{x}$$ is the vector of unknown node voltages and branch currents of voltage-defined devices.
* $$\mathbf{b}$$ is the vector of known independent source values (currents and voltages).

### Stamping Examples

For a **Resistor** ($R$) with conductance $g = 1/R$ connected between nodes $i$ and $j$:
$$G_{i,i} \mathrel{+}= g \quad G_{i,j} \mathrel{-}= g$$
$$G_{j,i} \mathrel{-}= g \quad G_{j,j} \mathrel{+}= g$$

For a **Voltage Source** ($V$) between nodes $i$ and $j$, an auxiliary current variable $k$ is introduced:
$$G_{i,k} = 1 \quad G_{k,i} = 1$$
$$G_{j,k} = -1 \quad G_{k,j} = -1$$
$$b_k = V$$

---

## 🏗️ Software Architecture

```text
e-XSpice/
├── src/
│   ├── core/           # Circuit topology, Node management, Device base classes
│   ├── devices/        # Physical models (Resistor, VoltageSource, etc.)
│   ├── analysis/       # Simulation engines (DCAnalysis, ACAnalysis...)
│   ├── solver/         # Numerical linear algebra (BiCGSTAB, CG, CSR/COO structures)
│   ├── parallel/       # MPI communicator wrappers and OpenMP pragmas
│   └── main.cpp        # CLI entry point
├── tests/              # CTest-based unit tests
├── examples/           # Benchmark circuits and netlists
└── CMakeLists.txt      # CMake configuration
```
---

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

### 2. Environment Setup

#### Windows
- Ensure your C++ compiler (MinGW or MSVC via Visual Studio Build Tools) is in your `PATH`
- If using VS Code, the **CMake Tools** extension is highly recommended

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install cmake build-essential libopenmpi-dev openmpi-bin
```

### Cross-Platform Compilation

Open your terminal (or Developer Command Prompt / PowerShell on Windows) and run:

### Clone the repository
```text
git clone https://github.com/yourusername/e-XSpice.git
cd e-XSpice
```

### Configure the project
```text
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

### Build using all cores
```text
cmake --build build --config Release -j
```

After compilation, the executable will be located at:

```text
./build/e-xspice        (Linux/macOS)
./build/e-xspice.exe    (Windows)
```

---

## 💻 Usage Examples

In the examples folder you are going to find some examples circuits.

Running the CLI Engine

### Example 1: Serial execution

```text
./build/e-xspice 1
```

### Example 2: Distributed MPI execution (4 processes and 1000 resistors)

```text
mpiexec -np 4 ./build/e-xspice 3 1000
```

Note:

On Linux, you may use **mpirun** instead of **mpiexec**
On Windows (MS-MPI), always use **mpiexec**

### OpenMP Thread Control

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
Total computation footprint:

```text
2 MPI processes × 4 OpenMP threads = 8 logical cores
```
---

## ⚡ Parallelization Strategy

### MPI (Distributed Memory)

- Distributes circuit devices across processes
- Each process assembles a local portion of the global conductance matrix (G)
- Uses MPI_Allreduce for synchronization during iterative solving

### OpenMP (Shared Memory)

- Parallelizes matrix assembly within each MPI rank
- Accelerates Sparse Matrix-Vector multiplication (SpMV) in solvers (CG / BiCGSTAB)

### Scalability Profile

- For large circuits (>10k nodes), hybrid parallelism:

```text
4 MPI ranks × 4 OpenMP threads
```

achieves ~10–12× speedup vs serial execution


## 🗺️ Next Steps (Development Roadmap)
 - [ ] Non-linear devices (Diodes + Newton-Raphson)
 - [ ] Transient Analysis (time integration)
 - [ ] AC Analysis (frequency domain)
 - [ ] MOSFET (BSIM) models
 - [ ] SPICE netlist parser (.cir, .sp)
 - [ ] Preconditioning (ILU)
 - [ ] Graph partitioning (METIS)


### 📜 License

This project is licensed under the GPL-3.0 License.
See the ```LICENSE``` file for details.

## 📚 References

### Core Resources
- [SPICE (UC Berkeley)](https://people.eecs.berkeley.edu/~newton/spice/)
- [ngspice] ()
- [Xyce Parallel Electronic Simulator](https://xyce.sandia.gov/)

### Books (Circuit Simulation & SPICE)
- Vladimirescu, A. — *The SPICE Book* (McGraw-Hill, 1994)  
- Najm, F. N. — *Circuit Simulation* (Wiley-IEEE Press, 2010) 
- Tuma, T., Burmen, Á. — *Circuit Simulation with SPICE OPUS: Theory and Practice* (Birkhäuser, 2009)

### Books (Numerical Methods & Advanced Simulation)
- Benner, P., Hinze, M., ter Maten, E. — *Model Reduction for Circuit Simulation* (Springer, 2011) 

### High Performance Computing (HPC) 
- Hager, G., Wellein, G. — *Introduction to High Performance Computing for Scientists and Engineers* (CRC Press, 2010)