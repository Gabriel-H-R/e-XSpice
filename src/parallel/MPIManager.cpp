/**
 * @author Gabriel Henrique Silva
 * @file MPIManager.cpp
 * @brief MPIManager implementation
 */

#include "MPIManager.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace e_XSpice {

MPIManager::MPIManager()
    : initialized_(false)
    , rank_(0)
    , size_(1)
    , name_len_(0)
{
    processor_name_[0] = '\0';
}

MPIManager& MPIManager::getInstance() {
    static MPIManager instance;
    return instance;
}

void MPIManager::initialize(int* argc, char*** argv) {
    if (initialized_) {
        return;
    }
    
    int provided;
    MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    
    if (provided < MPI_THREAD_FUNNELED) {
        std::cerr << "Warning: MPI thread support level lower than requested" 
                  << std::endl;
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    MPI_Get_processor_name(processor_name_, &name_len_);
    
    initialized_ = true;
}

void MPIManager::finalize() {
    if (initialized_) {
        MPI_Finalize();
        initialized_ = false;
    }
}

void MPIManager::barrier() const {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPIManager::broadcast(std::vector<double>& data) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    
    // First broadcast of size
    int size = data.size();
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Resizes non-master processes
    if (!isMaster()) {
        data.resize(size);
    }
    
    // Broadcast of data
    if (size > 0) {
        MPI_Bcast(data.data(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void MPIManager::broadcast(int& value) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void MPIManager::gather(const std::vector<double>& local_data,
                        std::vector<double>& global_data) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    
    int local_size = local_data.size();
    
    // Collects sizes from each process
    std::vector<int> sizes(size_);
    MPI_Gather(&local_size, 1, MPI_INT, 
               sizes.data(), 1, MPI_INT, 
               0, MPI_COMM_WORLD);
    
    // Calculates offsets
    std::vector<int> displs(size_, 0);
    int total_size = 0;
    
    if (isMaster()) {
        for (int i = 0; i < size_; ++i) {
            displs[i] = total_size;
            total_size += sizes[i];
        }
        global_data.resize(total_size);
    }
    
    // Gather of data
    MPI_Gatherv(local_data.data(), local_size, MPI_DOUBLE,
                global_data.data(), sizes.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
}

void MPIManager::scatter(const std::vector<double>& global_data,
                         std::vector<double>& local_data,
                         int local_size) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    
    local_data.resize(local_size);
    
    // Calculates sizes and offsets
    std::vector<int> sizes(size_);
    std::vector<int> displs(size_, 0);
    
    if (isMaster()) {
        int offset = 0;
        for (int i = 0; i < size_; ++i) {
            sizes[i] = local_size;
            displs[i] = offset;
            offset += local_size;
        }
    }
    
    MPI_Scatterv(global_data.data(), sizes.data(), displs.data(), MPI_DOUBLE,
                 local_data.data(), local_size, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
}

void MPIManager::allReduceSum(const std::vector<double>& local_data,
                              std::vector<double>& global_data) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    
    global_data.resize(local_data.size());
    
    MPI_Allreduce(local_data.data(), global_data.data(), 
                  local_data.size(), MPI_DOUBLE, 
                  MPI_SUM, MPI_COMM_WORLD);
}

void MPIManager::getLocalRange(int total_items, int& start, int& end) const {
    if (total_items <= 0) {
        start = end = 0;
        return;
    }
    
    // Balanced division
    int base_chunk = total_items / size_;
    int remainder = total_items % size_;
    
    if (rank_ < remainder) {
        // First processes get an extra item
        start = rank_ * (base_chunk + 1);
        end = start + base_chunk + 1;
    } else {
        start = rank_ * base_chunk + remainder;
        end = start + base_chunk;
    }
}

void MPIManager::printInfo() const {
    if (!initialized_) {
        std::cout << "MPI not initialized" << std::endl;
        return;
    }
    
    // Each process prints its information sequentially
    for (int i = 0; i < size_; ++i) {
        if (rank_ == i) {
            std::cout << "Process " << rank_ << "/" << size_ 
                     << " on " << processor_name_;
            if (isMaster()) {
                std::cout << " [MASTER]";
            }
            std::cout << std::endl;
        }
        barrier(); // Synchronizes for ordered printing
    }
}

} // namespace e_XSpice
