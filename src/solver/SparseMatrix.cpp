/**
 * @author Gabriel Henrique Silva
 * @file SparseMatrix.cpp
 * @brief SparseMatrix class implementation
 */

#include "SparseMatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace e_XSpice {

SparseMatrix::SparseMatrix(int rows, int cols)
    : rows_(rows)
    , cols_(cols)
    , format_(MatrixFormat::COO)
{
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
}

bool SparseMatrix::isValid(int row, int col) const {
    return row >= 0 && row < rows_ && col >= 0 && col < cols_;
}

/**
 * @brief Sum a value to the position(row, col). If the position already has a value, it will be accumulated.
 */
void SparseMatrix::add(int row, int col, double value) {
    if (!isValid(row, col)) {
        throw std::out_of_range("Matrix indices out of range");
    }
    
    if (format_ != MatrixFormat::COO) {
        throw std::runtime_error("Can only add to COO format. Convert back first.");
    }
    
    // Ignores very small values (numeric threshold)
    if (std::abs(value) < 1e-15) {
        return;
    }
    
    auto key = std::make_pair(row, col);
    coo_data_[key] += value; // Accumulates if already exists
}

void SparseMatrix::set(int row, int col, double value) {
    if (!isValid(row, col)) {
        throw std::out_of_range("Matrix indices out of range");
    }
    
    if (format_ != MatrixFormat::COO) {
        throw std::runtime_error("Can only set in COO format");
    }
    
    if (std::abs(value) < 1e-15) {
        // Removes if value is too small
        auto key = std::make_pair(row, col);
        coo_data_.erase(key);
        return;
    }
    
    auto key = std::make_pair(row, col);
    coo_data_[key] = value; // Replaces
}

double SparseMatrix::get(int row, int col) const {
    if (!isValid(row, col)) {
        throw std::out_of_range("Matrix indices out of range");
    }
    
    if (format_ == MatrixFormat::COO) {
        auto key = std::make_pair(row, col);
        auto it = coo_data_.find(key);
        return (it != coo_data_.end()) ? it->second : 0.0;
    } else {
        // CSR format
        int row_start = csr_row_ptr_[row];
        int row_end = csr_row_ptr_[row + 1];
        
        for (int i = row_start; i < row_end; ++i) {
            if (csr_col_indices_[i] == col) {
                return csr_values_[i];
            }
        }
        return 0.0;
    }
}

void SparseMatrix::clear() {
    coo_data_.clear();
    csr_values_.clear();
    csr_col_indices_.clear();
    csr_row_ptr_.clear();
    format_ = MatrixFormat::COO;
}

int SparseMatrix::getNonZeros() const {
    if (format_ == MatrixFormat::COO) {
        return coo_data_.size();
    } else {
        return csr_values_.size();
    }
}

void SparseMatrix::convertToCSR() {
    if (format_ == MatrixFormat::CSR) {
        return; // Already in CSR
    }
    
    //Clear previous CSR data/structure
    csr_values_.clear();
    csr_col_indices_.clear();
    csr_row_ptr_.clear();
    csr_row_ptr_.resize(rows_ + 1, 0);
    
    // Count elements per line to build row pointers
    std::vector<int> row_counts(rows_, 0);
    for (const auto& entry : coo_data_) {
        row_counts[entry.first.first]++;
    }
    
    // Calculate row pointers (prefixsum)
    csr_row_ptr_[0] = 0;
    for (int i = 0; i < rows_; ++i) {
        csr_row_ptr_[i + 1] = csr_row_ptr_[i] + row_counts[i];
    }
    
    // Allocate space
    int nnz = coo_data_.size();
    csr_values_.resize(nnz);
    csr_col_indices_.resize(nnz);
    
    //Filling CSR values and column indices
    std::vector<int> current_pos = csr_row_ptr_;
    current_pos.pop_back(); // Remove last element
    
    for (const auto& entry : coo_data_) {
        int row = entry.first.first;
        int col = entry.first.second;
        double val = entry.second;
        
        int pos = current_pos[row]++;
        csr_values_[pos] = val;
        csr_col_indices_[pos] = col;
    }
    
    // Order by column index inside each row for better cache performance in multiplication
    // Uses OpenMP to parallelize the sorting of rows
    #pragma omp parallel for
    for (int row = 0; row < rows_; ++row) {
        int row_start = csr_row_ptr_[row];
        int row_end = csr_row_ptr_[row + 1];
        
        if (row_end - row_start > 1) {
            // Create pairs (column, value) for sorting
            std::vector<std::pair<int, double>> row_data;
            for (int i = row_start; i < row_end; ++i) {
                row_data.push_back({csr_col_indices_[i], csr_values_[i]});
            }
            
            //Order by column
            std::sort(row_data.begin(), row_data.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
            
            // Writing back sorted data to CSR structure
            for (size_t i = 0; i < row_data.size(); ++i) {
                csr_col_indices_[row_start + i] = row_data[i].first;
                csr_values_[row_start + i] = row_data[i].second;
            }
        }
    }
    
    format_ = MatrixFormat::CSR;
}

void SparseMatrix::multiplyVector(const std::vector<double>& x, 
                                  std::vector<double>& result) const {
    if (format_ != MatrixFormat::CSR) {
        throw std::runtime_error("Matrix must be in CSR format for multiplication");
    }
    
    if (static_cast<int>(x.size()) != cols_) {
        throw std::invalid_argument("Vector size must match matrix columns");
    }
    
    result.resize(rows_, 0.0);
    
    // Parallel multiplication with OpenMP
    #pragma omp parallel for schedule(dynamic, 64)
    for (int row = 0; row < rows_; ++row) {
        double sum = 0.0;
        int row_start = csr_row_ptr_[row];
        int row_end = csr_row_ptr_[row + 1];
        
        // Scalar product of line with vector
        for (int i = row_start; i < row_end; ++i) {
            sum += csr_values_[i] * x[csr_col_indices_[i]];
        }
        
        result[row] = sum;
    }
}

void SparseMatrix::print(bool dense) const {
    std::cout << "SparseMatrix " << rows_ << "x" << cols_ 
              << " (" << getNonZeros() << " non-zeros)" << std::endl;
    
    if (dense) {
        // Prints as dense matrix
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                std::cout << std::setw(10) << std::setprecision(4) 
                         << get(i, j) << " ";
            }
            std::cout << std::endl;
        }
    } else {
        // Prints non-zero elements
        if (format_ == MatrixFormat::COO) {
            for (const auto& entry : coo_data_) {
                std::cout << "(" << entry.first.first << "," 
                         << entry.first.second << ") = " 
                         << entry.second << std::endl;
            }
        } else {
            for (int row = 0; row < rows_; ++row) {
                int row_start = csr_row_ptr_[row];
                int row_end = csr_row_ptr_[row + 1];
                for (int i = row_start; i < row_end; ++i) {
                    std::cout << "(" << row << "," << csr_col_indices_[i] 
                             << ") = " << csr_values_[i] << std::endl;
                }
            }
        }
    }
}

void SparseMatrix::printStats() const {
    std::cout << "=== Sparse Matrix Statistics ===" << std::endl;
    std::cout << "Dimensions: " << rows_ << " x " << cols_ << std::endl;
    std::cout << "Non-zeros: " << getNonZeros() << std::endl;
    std::cout << "Sparsity: " << (1.0 - (double)getNonZeros() / (rows_ * cols_)) * 100 
              << "%" << std::endl;
    std::cout << "Format: " << (format_ == MatrixFormat::COO ? "COO" : "CSR") 
              << std::endl;
    std::cout << "=================================" << std::endl;
}

} // namespace e_XSpice
