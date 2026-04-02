/**
 * @author Gabriel Henrique Silva
 * @file MPIManager.h
 * @brief MPI communication management for distributed parallelization
 * 
 * This class encapsulates all MPI operations of the simulator, including:
 * - Initialization and finalization
 * - Circuit distribution among processes
 * - Communication of vectors and matrices
 * - Result synchronization
 */

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

//#include <mpi.h>
#include "D:/Program Files (x86)/Microsoft SDKs/MPI/Include/mpi.h"
#include <vector>
#include <string>

namespace e_XSpice {

/**
 * @brief Singleton to manage MPI
 */
class MPIManager {
public:
    /**
     * @brief Gets unique instance (Singleton)
     */
    static MPIManager& getInstance();
    
    /**
     * @brief Initializes MPI
     * @param argc Argc from main
     * @param argv Argv from main
     */
    void initialize(int* argc, char*** argv);
    
    /**
     * @brief Finalizes MPI
     */
    void finalize();
    
    /**
     * @brief Checks if MPI was initialized
     */
    bool isInitialized() const { return initialized_; }
    
    /**
     * @brief Returns rank of current process
     */
    int getRank() const { return rank_; }
    
    /**
     * @brief Returns total number of processes
     */
    int getSize() const { return size_; }
    
    /**
     * @brief Checks if it is the master process (rank 0)
     */
    bool isMaster() const { return rank_ == 0; }
    
    /**
     * @brief Synchronization barrier
     */
    void barrier() const;
    
    /**
     * @brief Broadcast vector from master to all
     * @param data Vector to be transmitted
     */
    void broadcast(std::vector<double>& data);
    
    /**
     * @brief Broadcast of integer
     */
    void broadcast(int& value);
    
    /**
     * @brief Gather: collects vectors from all processes in master
     * @param local_data Local data
     * @param global_data Global data (only valid on master)
     */
    void gather(const std::vector<double>& local_data,
                std::vector<double>& global_data);
    
    /**
     * @brief Scatter: distributes vector from master to all
     * @param global_data Global data (on master)
     * @param local_data Local data (on all)
     * @param local_size Size of each local chunk
     */
    void scatter(const std::vector<double>& global_data,
                std::vector<double>& local_data,
                int local_size);
    
    /**
     * @brief AllReduce: sums vectors from all processes
     * @param local_data Local data
     * @param global_data Result of sum (on all)
     */
    void allReduceSum(const std::vector<double>& local_data,
                      std::vector<double>& global_data);
    
    /**
     * @brief Calculates range of indices for this process
     * 
     * Divides total_items between size_ processes and returns
     * [start, end) for this rank.
     * 
     * @param total_items Total items to divide
     * @param start Initial index (output)
     * @param end Final index (output, exclusive)
     */
    void getLocalRange(int total_items, int& start, int& end) const;
    
    /**
     * @brief Prints process information
     */
    void printInfo() const;

private:
    MPIManager(); // Private constructor (Singleton)
    ~MPIManager() = default;
    
    // Forbid copy
    MPIManager(const MPIManager&) = delete;
    MPIManager& operator=(const MPIManager&) = delete;
    
    bool initialized_;
    int rank_;        // Process ID
    int size_;        // Total number of processes
    
    char processor_name_[MPI_MAX_PROCESSOR_NAME];
    int name_len_;
};

} // namespace e_XSpice

#endif // MPI_MANAGER_H
