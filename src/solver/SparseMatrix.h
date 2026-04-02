/**
 * @author Gabriel Henrique Silva
 * @file SparseMatrix.h
 * @brief Sparse matrix implementation for circuit analysis
 * 
 * Supports two formats:
 * - COO (Coordinate): for efficient construction
 * - CSR (Compressed Sparse Row): for fast operations
 * 
 * OpenMP parallelization for matrix-vector operations.
 */

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <map>
#include <utility>
#include <omp.h>

namespace e_XSpice {

/**
 * @brief Matrix storage format
 */
enum class MatrixFormat {
    COO,  // Coordinate format (triplets)
    CSR   // Compressed Sparse Row
};

/**
 * @brief Optimized sparse matrix class
 */
class SparseMatrix {
public:
    /**
     * @brief Constructor
     * @param rows Number of rows
     * @param cols Number of columns
     */
    SparseMatrix(int rows, int cols);
    
    /**
     * @brief Destructor
     */
    ~SparseMatrix() = default;
    
    // Basic getters
    int getRows() const { return rows_; }
    int getCols() const { return cols_; }
    int getNonZeros() const;
    MatrixFormat getFormat() const { return format_; }
    
    /**
     * @brief Adds value to matrix (accumulates if already exists)
     * @param row Row (0-indexed)
     * @param col Column (0-indexed)
     * @param value Value to add
     */
    void add(int row, int col, double value);
    
    /**
     * @brief Sets value in matrix (replaces if already exists)
     */
    void set(int row, int col, double value);
    
    /**
     * @brief Gets matrix value
     */
    double get(int row, int col) const;
    
    /**
     * @brief Clears entire matrix
     */
    void clear();
    
    /**
     * @brief Converts to CSR format (necessary before solving)
     */
    void convertToCSR();
    
    /**
     * @brief Multiplies matrix by vector: result = A * x
     * Parallelized with OpenMP
     * 
     * @param x Input vector
     * @param result Output vector (must have size rows_)
     */
    void multiplyVector(const std::vector<double>& x, 
                       std::vector<double>& result) const;
    
    /**
     * @brief Prints matrix (for debug)
     * @param dense If true, prints as dense matrix
     */
    void print(bool dense = false) const;
    
    /**
     * @brief Matrix statistics
     */
    void printStats() const;
    
    // Direct access to CSR data (for external solvers)
    const std::vector<double>& getValues() const { return csr_values_; }
    const std::vector<int>& getColIndices() const { return csr_col_indices_; }
    const std::vector<int>& getRowPointers() const { return csr_row_ptr_; }

private:
    int rows_;                    // Number of rows
    int cols_;                    // Number of columns
    MatrixFormat format_;         // Current format
    
    // COO Format (Coordinate) - used during construction
    std::map<std::pair<int,int>, double> coo_data_;
    
    // CSR Format (Compressed Sparse Row) - used for calculations
    std::vector<double> csr_values_;      // Non-zero values
    std::vector<int> csr_col_indices_;    // Column indices
    std::vector<int> csr_row_ptr_;        // Row start pointers
    
    /**
     * @brief Checks if indices are within bounds
     */
    bool isValid(int row, int col) const;
};

} // namespace e_XSpice

#endif // SPARSE_MATRIX_H
