/**
 * @author Gabriel Henrique Silva
 * @file DCAnalysis.cpp
 * @brief Implementation of DC Operating Point Analysis
 */

#include "DCAnalysis.h"
#include "../devices/VoltageSource.h"
#include <iostream>
#include <chrono>
#include <algorithm>

namespace e_XSpice {

DCAnalysis::DCAnalysis(Circuit& circuit)
    : circuit_(circuit)
    , solver_type_(SolverType::BICGSTAB)
    , max_iterations_(10000)
    , tolerance_(1e-9)
    , verbose_(false)
{
}

void DCAnalysis::createNodeMap(std::vector<int>& node_map) {
    // Creates mapping: node_id -> matrix_index
    // Ground (id=0) does not enter the matrix
    
    int num_nodes = circuit_.getNumNodes();
    node_map.resize(num_nodes, -1);
    
    // Ground always maps to -1 (not in the matrix)
    node_map[0] = -1;
    
    // Other nodes map sequentially
    const auto& nodes = circuit_.getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_map[nodes[i]->getId()] = i;
    }
}

void DCAnalysis::assignVoltageSourceIndices(std::vector<int>& node_map) {
    // Voltage sources need extra variables for current
    // Their indices come after the node voltages
    
    int num_nodes = circuit_.getNodes().size();  // Excluding ground
    int vsource_idx = num_nodes;  // Index starts after the nodes
    
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::VOLTAGE_SOURCE) {
            VoltageSource* vs = dynamic_cast<VoltageSource*>(device);
            if (vs) {
                vs->setCurrentIndex(vsource_idx);
                vsource_idx++;
            }
        }
    }
}

void DCAnalysis::assembleMNA_Serial(SparseMatrix& G,
                                    std::vector<double>& b,
                                    std::vector<int>& node_map) {
    // Initializes vectors
    int matrix_size = circuit_.getNodes().size() + circuit_.countVoltageSources();
    b.assign(matrix_size, 0.0);
    
    // Each device stamps its contribution
    for (auto device : circuit_.getDevices()) {
        device->stamp(G, b, node_map);
    }
}

void DCAnalysis::assembleMNA_MPI(SparseMatrix& G,
                                 std::vector<double>& b,
                                 std::vector<int>& node_map) {
    MPIManager& mpi = MPIManager::getInstance();
    
    int matrix_size = circuit_.getNodes().size() + circuit_.countVoltageSources();
    b.assign(matrix_size, 0.0);
    
    // Each process works on a subset of devices
    const auto& devices = circuit_.getDevices();
    int num_devices = devices.size();
    
    int start_dev, end_dev;
    mpi.getLocalRange(num_devices, start_dev, end_dev);
    
    // Local matrix for this process
    SparseMatrix G_local(matrix_size, matrix_size);
    std::vector<double> b_local(matrix_size, 0.0);
    
    // OpenMP parallelization within the MPI range
    #pragma omp parallel
    {
        // Each thread has its own temporary matrix
        SparseMatrix G_thread(matrix_size, matrix_size);
        std::vector<double> b_thread(matrix_size, 0.0);
        
        #pragma omp for schedule(dynamic)
        for (int i = start_dev; i < end_dev; ++i) {
            devices[i]->stamp(G_thread, b_thread, node_map);
        }
        
        // Combines thread contributions (critical region)
        #pragma omp critical
        {
            for (int row = 0; row < matrix_size; ++row) {
                for (int col = 0; col < matrix_size; ++col) {
                    double val = G_thread.get(row, col);
                    if (std::abs(val) > 1e-15) {
                        G_local.add(row, col, val);
                    }
                }
                b_local[row] += b_thread[row];
            }
        }
    }
    
    // MPI reduction: combines matrices from all processes
    if (mpi.getSize() > 1) {
        // AllReduce of vector b
        std::vector<double> b_global(matrix_size);
        mpi.allReduceSum(b_local, b_global);
        b = b_global;
        
        // For the matrix, we need to transfer COO triplets
        // Simplification: master process collects everything
        if (mpi.isMaster()) {
            // Master already has G_local, now receives from others
            for (int proc = 1; proc < mpi.getSize(); ++proc) {
                // Receive matrix from process proc and add to G
                // (full implementation requires SparseMatrix serialization)
                // For now, we assume serial or OpenMP-only
            }
            G = G_local;
        } else {
            // Send G_local to master
            // (MPI implementation of sparse matrix)
        }
        
        // Broadcast of the final matrix
        mpi.barrier();
        
    } else {
        // Only one process
        G = G_local;
        b = b_local;
    }
}

void DCAnalysis::assembleMNA(SparseMatrix& G,
                             std::vector<double>& b,
                             std::vector<int>& node_map,
                             bool use_mpi) {
    createNodeMap(node_map);
    assignVoltageSourceIndices(node_map);
    
    int matrix_size = circuit_.getNodes().size() + circuit_.countVoltageSources();
    G = SparseMatrix(matrix_size, matrix_size);
    
    if (use_mpi && MPIManager::getInstance().isInitialized() && 
        MPIManager::getInstance().getSize() > 1) {
        assembleMNA_MPI(G, b, node_map);
    } else {
        assembleMNA_Serial(G, b, node_map);
    }
    
    // Converts to CSR for solver
    G.convertToCSR();
}

void DCAnalysis::extractSolution(const std::vector<double>& x,
                                 const std::vector<int>& node_map) {
    // Extracts node voltages
    const auto& nodes = circuit_.getNodes();
    for (auto node : nodes) {
        int idx = node_map[node->getId()];
        if (idx >= 0 && idx < static_cast<int>(x.size())) {
            node->setVoltage(x[idx]);
        }
    }
    
    // Extracts currents from voltage sources
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::VOLTAGE_SOURCE) {
            VoltageSource* vs = dynamic_cast<VoltageSource*>(device);
            if (vs) {
                int idx = vs->getCurrentIndex();
                if (idx >= 0 && idx < static_cast<int>(x.size())) {
                    vs->setCurrent(x[idx]);
                }
            }
        }
    }
    
    // Saves solution for next iteration
    previous_solution_ = x;
}

DCResult DCAnalysis::run(bool use_mpi, bool use_initial_guess) {
    DCResult result;
    result.success = false;
    
    // Validates circuit
    std::string error_msg;
    if (!circuit_.validate(error_msg)) {
        result.message = "Circuit validation failed: " + error_msg;
        return result;
    }
    
    if (verbose_) {
        MPIManager& mpi = MPIManager::getInstance();
        if (!mpi.isInitialized() || mpi.isMaster()) {
            std::cout << "\n===== DC Operating Point Analysis =====" << std::endl;
            std::cout << "Nodes: " << circuit_.getNumNodes() << std::endl;
            std::cout << "Devices: " << circuit_.getNumDevices() << std::endl;
            std::cout << "Voltage Sources: " << circuit_.countVoltageSources() << std::endl;
        }
    }
    
    // Assembly of the MNA system
    auto assembly_start = std::chrono::high_resolution_clock::now();
    
    SparseMatrix G(1, 1); // Temporary initialization
    std::vector<double> b;
    std::vector<int> node_map;
    
    assembleMNA(G, b, node_map, use_mpi);
    
    auto assembly_end = std::chrono::high_resolution_clock::now();
    result.assembly_time = std::chrono::duration<double>(assembly_end - assembly_start).count();
    
    if (verbose_) {
        MPIManager& mpi = MPIManager::getInstance();
        if (!mpi.isInitialized() || mpi.isMaster()) {
            std::cout << "Matrix assembly time: " << result.assembly_time << " s" << std::endl;
            G.printStats();
        }
    }
    
    // Resolution of the linear system
    LinearSolver solver(solver_type_);
    solver.setMaxIterations(max_iterations_);
    solver.setTolerance(tolerance_);
    solver.setVerbose(verbose_);
    
    std::vector<double> x;
    const std::vector<double>* initial_guess = nullptr;
    if (use_initial_guess && !previous_solution_.empty()) {
        initial_guess = &previous_solution_;
    }
    
    auto solve_start = std::chrono::high_resolution_clock::now();
    result.solver_result = solver.solve(G, b, x, initial_guess);
    auto solve_end = std::chrono::high_resolution_clock::now();
    
    result.solve_time = std::chrono::duration<double>(solve_end - solve_start).count();
    
    if (!result.solver_result.converged) {
        result.message = "Solver failed to converge: " + result.solver_result.message;
        return result;
    }
    
    // Extracts solution
    extractSolution(x, node_map);
    
    result.success = true;
    result.message = "DC analysis completed successfully";
    result.node_voltages = x;
    
    if (verbose_) {
        MPIManager& mpi = MPIManager::getInstance();
        if (!mpi.isInitialized() || mpi.isMaster()) {
            std::cout << "\nSolver: " << result.solver_result.message << std::endl;
            std::cout << "Iterations: " << result.solver_result.iterations << std::endl;
            std::cout << "Residual: " << result.solver_result.residual_norm << std::endl;
            std::cout << "Solve time: " << result.solve_time << " s" << std::endl;
            std::cout << "Total time: " << (result.assembly_time + result.solve_time) << " s" << std::endl;
            std::cout << "========================================\n" << std::endl;
        }
    }
    
    return result;
}

} // namespace e_XSpice
