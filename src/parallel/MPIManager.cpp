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
    double dummy = 0.0;
    const double* send_ptr = local_data.empty() ? &dummy : local_data.data();
    double* recv_ptr = global_data.empty() ? &dummy : global_data.data();

    MPI_Gatherv(send_ptr, local_size, MPI_DOUBLE,
                recv_ptr, sizes.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    //MPI_Gatherv(local_data.data(), local_size, MPI_DOUBLE,
    //            global_data.data(), sizes.data(), displs.data(), MPI_DOUBLE,
    //            0, MPI_COMM_WORLD);
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
    // Master collects the expected local_size from each process
    MPI_Gather(&local_size, 1, MPI_INT, 
               sizes.data(), 1, MPI_INT, 
               0, MPI_COMM_WORLD);

    std::vector<int> displs(size_, 0);
    
    if (isMaster()) {
        int offset = 0;
        for (int i = 0; i < size_; ++i) {
            //sizes[i] = local_size;
            displs[i] = offset;
            offset += sizes[i];
        }
    }
    
    // Protection to MS-MPI against acesses to empty .data()
    const double* send_ptr = global_data.empty() ? nullptr : global_data.data();
    double* recv_ptr = local_data.empty() ? nullptr : local_data.data();

    MPI_Scatterv(send_ptr, sizes.data(), displs.data(), MPI_DOUBLE,
                 recv_ptr, local_size, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
}

void MPIManager::allReduceSum(const std::vector<double>& local_data,
                              std::vector<double>& global_data) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }

    if (local_data.empty()) return;
    
    global_data.resize(local_data.size());
    
    MPI_Allreduce(local_data.data(), global_data.data(), 
                  local_data.size(), MPI_DOUBLE, 
                  MPI_SUM, MPI_COMM_WORLD);
}

void MPIManager::allReduceSumInPlace(std::vector<double>& data) {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    if (data.empty() || size_ <= 1) {
        return;   // nothing to combine (single rank: already global)
    }
    MPI_Allreduce(MPI_IN_PLACE, data.data(),
                  static_cast<int>(data.size()), MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
}

double MPIManager::allReduceSumDouble(double local) const {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    if (size_ <= 1) return local;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

void MPIManager::allGatherV(const std::vector<double>& local_segment,
                            std::vector<double>& global,
                            int totalSize) const {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    if (totalSize <= 0) {
        global.clear();
        return;
    }
    // Block partition is a pure function of (totalSize, rank): every rank
    // derives the same counts/displacements without communicating.
    // (Same layout as getLocalRange.)
    std::vector<int> counts(size_), displs(size_);
    for (int r = 0; r < size_; ++r) {
        int s = 0, e = 0;
        blockRangeForRank(totalSize, r, s, e);
        counts[r] = e - s;
        displs[r] = s;
    }
    global.assign(totalSize, 0.0);

    double dummy = 0.0;
    const double* send_ptr = local_segment.empty() ? &dummy : local_segment.data();
    MPI_Allgatherv(send_ptr, counts[rank_], MPI_DOUBLE,
                   global.data(), counts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
}

bool MPIManager::allReduceAnd(bool local) const {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    if (size_ <= 1) return local;
    int in = local ? 1 : 0, out = 0;
    MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return out != 0;
}

bool MPIManager::allReduceOr(bool local) const {
    if (!initialized_) {
        throw std::runtime_error("MPI not initialized");
    }
    if (size_ <= 1) return local;
    int in = local ? 1 : 0, out = 0;
    MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    return out != 0;
}

void MPIManager::getLocalRange(int total_items, int& start, int& end) const {
    blockRangeForRank(total_items, rank_, start, end);
}

void MPIManager::blockRangeForRank(int total_items, int rank,
                                   int& start, int& end) const {
    if (total_items <= 0) {
        start = end = 0;
        return;
    }

    // Balanced division
    int base_chunk = total_items / size_;
    int remainder = total_items % size_;

    if (rank < remainder) {
        // First processes get an extra item
        start = rank * (base_chunk + 1);
        end = start + base_chunk + 1;
    } else {
        start = rank * base_chunk + remainder;
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
