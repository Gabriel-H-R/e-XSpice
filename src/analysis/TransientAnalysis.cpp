/**
 * @author Gabriel Henrique Silva
 * @file TransientAnalysis.cpp
 * @brief Implementation of transient analysis
 */

#include "TransientAnalysis.h"
#include "../devices/VoltageSource.h"
#include "../devices/Capacitor.h"
#include "../devices/Inductor.h"
#include <iostream>
#include <chrono>
#include <cmath>

namespace e_XSpice {

TransientAnalysis::TransientAnalysis(Circuit& circuit)
    : circuit_(circuit)
    , method_(IntegrationMethod::TRAPEZOIDAL)
    , solver_type_(SolverType::BICGSTAB)
    , max_iterations_(10000)
    , tolerance_(1e-9)
    , verbose_(false)
    , save_interval_(1)
{
}

void TransientAnalysis::createNodeMap(std::vector<int>& node_map) {
    int num_nodes = circuit_.getNumNodes();
    node_map.resize(num_nodes, -1);
    node_map[0] = -1;  // Ground
    
    const auto& nodes = circuit_.getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_map[nodes[i]->getId()] = i;
    }
    
    // Assigns indices for voltage sources and inductors
    int num_node_vars = nodes.size();
    int extra_idx = num_node_vars;
    
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::VOLTAGE_SOURCE) {
            VoltageSource* vs = dynamic_cast<VoltageSource*>(device);
            if (vs) {
                vs->setCurrentIndex(extra_idx++);
            }
        } else if (device->getType() == DeviceType::INDUCTOR) {
            Inductor* l = dynamic_cast<Inductor*>(device);
            if (l) {
                l->setCurrentIndex(extra_idx++);
            }
        }
    }
}

void TransientAnalysis::assembleTransientMNA(SparseMatrix& G,
                                             std::vector<double>& b,
                                             const std::vector<int>& node_map,
                                             double timestep,
                                             bool use_mpi) {
    // Counts extra variables (V sources and inductors)
    int num_vsources = circuit_.countVoltageSources();
    int num_inductors = 0;
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::INDUCTOR) {
            num_inductors++;
        }
    }
    
    int matrix_size = circuit_.getNodes().size() + num_vsources + num_inductors;
    b.assign(matrix_size, 0.0);
    
    // Stamp of each device
    const auto& devices = circuit_.getDevices();
    
    for (auto device : devices) {
        DeviceType type = device->getType();
        
        if (type == DeviceType::RESISTOR ||
            type == DeviceType::VOLTAGE_SOURCE ||
            type == DeviceType::CURRENT_SOURCE) {
            // Normal DC devices
            device->stamp(G, b, node_map);
            
        } else if (type == DeviceType::CAPACITOR) {
            // Capacitor uses companion model
            Capacitor* cap = dynamic_cast<Capacitor*>(device);
            if (cap) {
                cap->setTimeStep(timestep);
                cap->stampTransient(G, b, node_map, timestep);
            }
            
        } else if (type == DeviceType::INDUCTOR) {
            // Inductor uses companion model
            Inductor* ind = dynamic_cast<Inductor*>(device);
            if (ind) {
                ind->setTimeStep(timestep);
                ind->stampTransient(G, b, node_map, timestep);
            }
        }
    }
    
    (void)use_mpi;  // MPI for assembly will be implemented later
}

void TransientAnalysis::extractSolution(const std::vector<double>& x,
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
}

void TransientAnalysis::updateDeviceStates() {
    // Updates state of reactive devices for next step
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::CAPACITOR ||
            device->getType() == DeviceType::INDUCTOR) {
            device->updateState();
        }
    }
}

TransientDataPoint TransientAnalysis::collectDataPoint(double time,
                                                       const std::vector<int>& node_map) {
    TransientDataPoint point;
    point.time = time;
    
    // Collects node voltages
    const auto& nodes = circuit_.getNodes();
    point.node_voltages.resize(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        point.node_voltages[i] = nodes[i]->getVoltage();
    }
    
    // Collects currents from sources and inductors
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::VOLTAGE_SOURCE ||
            device->getType() == DeviceType::INDUCTOR) {
            point.branch_currents.push_back(device->getCurrent());
        }
    }
    
    (void)node_map;
    return point;
}

bool TransientAnalysis::timeStep(double t, double dt,
                                const std::vector<int>& node_map,
                                bool use_mpi) {
    // Assembles system for this time step
    int matrix_size = circuit_.getNodes().size() + 
                     circuit_.countVoltageSources();
    
    // Counts inductors
    for (auto device : circuit_.getDevices()) {
        if (device->getType() == DeviceType::INDUCTOR) {
            matrix_size++;
        }
    }
    
    SparseMatrix G(matrix_size, matrix_size);
    std::vector<double> b;
    
    assembleTransientMNA(G, b, node_map, dt, use_mpi);
    G.convertToCSR();
    
    // Solves linear system
    LinearSolver solver(solver_type_);
    solver.setMaxIterations(max_iterations_);
    solver.setTolerance(tolerance_);
    solver.setVerbose(false);  // Not verbose in each step
    
    std::vector<double> x;
    SolverResult result = solver.solve(G, b, x);
    
    if (!result.converged) {
        if (verbose_) {
            std::cerr << "Warning: Solver did not converge at t=" << t << "s" << std::endl;
        }
        return false;
    }
    
    // Extracts solution
    extractSolution(x, node_map);
    
    // Updates device states for next step
    updateDeviceStates();
    
    return true;
}

TransientResult TransientAnalysis::run(double tstart, double tstop, double tstep,
                                      bool use_mpi) {
    TransientResult result;
    result.success = false;
    
    // Validates parameters
    if (tstep <= 0) {
        result.message = "Time step must be positive";
        return result;
    }
    
    if (tstop <= tstart) {
        result.message = "Stop time must be greater than start time";
        return result;
    }
    
    // Validates circuit
    std::string error_msg;
    if (!circuit_.validate(error_msg)) {
        result.message = "Circuit validation failed: " + error_msg;
        return result;
    }
    
    MPIManager& mpi = MPIManager::getInstance();
    if (verbose_ && (!mpi.isInitialized() || mpi.isMaster())) {
        std::cout << "\n===== Transient Analysis =====" << std::endl;
        std::cout << "Time: " << tstart << "s to " << tstop << "s" << std::endl;
        std::cout << "Step: " << tstep << "s" << std::endl;
        std::cout << "Method: " << (method_ == IntegrationMethod::TRAPEZOIDAL ? 
                                    "Trapezoidal" : "Backward Euler") << std::endl;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Creates node mapping
    std::vector<int> node_map;
    createNodeMap(node_map);
    
    // Calculates number of steps
    int num_steps = static_cast<int>(std::ceil((tstop - tstart) / tstep));
    
    // Reserves space for data
    result.data.reserve(num_steps / save_interval_ + 1);
    
    // Temporal simulation
    double t = tstart;
    int step = 0;
    
    // Saves initial condition
    result.data.push_back(collectDataPoint(t, node_map));
    
    while (t < tstop) {
        double dt = std::min(tstep, tstop - t);
        
        if (!timeStep(t, dt, node_map, use_mpi)) {
            result.message = "Time step failed at t=" + std::to_string(t) + "s";
            return result;
        }
        
        t += dt;
        step++;
        
        // Saves data (if interval allows)
        if (step % save_interval_ == 0) {
            result.data.push_back(collectDataPoint(t, node_map));
        }
        
        // Progress callback
        if (progress_callback_) {
            progress_callback_(t, step, num_steps);
        }
        
        // Prints progress
        if (verbose_ && step % 100 == 0 && (!mpi.isInitialized() || mpi.isMaster())) {
            std::cout << "  Step " << step << "/" << num_steps 
                     << " (t=" << t << "s)" << std::endl;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double>(end_time - start_time).count();
    
    result.success = true;
    result.message = "Transient analysis completed successfully";
    
    if (verbose_ && (!mpi.isInitialized() || mpi.isMaster())) {
        std::cout << "\nCompleted " << step << " time steps" << std::endl;
        std::cout << "Total time: " << result.total_time << " s" << std::endl;
        std::cout << "Time per step: " << (result.total_time / step * 1000) << " ms" << std::endl;
        std::cout << "==============================\n" << std::endl;
    }
    
    return result;
}

} // namespace e-XSpice
