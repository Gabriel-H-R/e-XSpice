/**
 * @author Gabriel Henrique Silva
 * @file DaeMatrix.cpp
 * @brief Fixed-sparsity-pattern sparse matrix for DAE-formulation assembly.
 */

#include "DaeMatrix.h"

#include <omp.h>
#include <algorithm>
#include <stdexcept>
#include <cstdio>

namespace e_XSpice {

DaeMatrix::DaeMatrix(int n)
    : n_(n) {
    if (n <= 0) {
        // Placeholder matrix: pattern phase becomes usable after assignPattern.
        finalized_ = true;   // empty but consistent
    }
}

void DaeMatrix::declareEntry(int row, int col) {
    if (row < 0 || row >= n_ || col < 0 || col >= n_) {
        throw std::out_of_range("DaeMatrix::declareEntry out of range");
    }
    if (finalized_) {
        throw std::runtime_error("Cannot declare entries after finalizePattern()");
    }
    pattern_.emplace_back(row, col);
}

void DaeMatrix::declareEntries(const std::vector<int>& rows,
                               const std::vector<int>& cols) {
    if (rows.size() != cols.size()) {
        throw std::invalid_argument("DaeMatrix::declareEntries: size mismatch");
    }
    pattern_.reserve(pattern_.size() + rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        declareEntry(rows[i], cols[i]);
    }
}

void DaeMatrix::finalizePattern() {
    if (finalized_) return;

    // 1) Dedup declarations (a device may declare the same entry twice).
    std::sort(pattern_.begin(), pattern_.end());
    pattern_.erase(std::unique(pattern_.begin(), pattern_.end()), pattern_.end());

    // 2) Count non-zeros per row.
    nnz_ = static_cast<int>(pattern_.size());
    csr_row_ptr_.assign(n_ + 1, 0);

    #pragma omp parallel for
    for (int i = 0; i < n_; ++i) {
        const auto lo = std::lower_bound(pattern_.begin(), pattern_.end(),
                                          JacobianEntry(i, 0));
        const auto hi = std::lower_bound(pattern_.begin(), pattern_.end(),
                                          JacobianEntry(i + 1, 0));
        csr_row_ptr_[i]     = static_cast<int>(lo - pattern_.begin());
        csr_row_ptr_[i + 1] = static_cast<int>(hi - pattern_.begin());
    }

    // 3) Fill column indices (pattern_ is sorted by (row,col)) and values.
    csr_col_idx_.resize(nnz_);
    csr_values_.assign(nnz_, 0.0);

    #pragma omp parallel for
    for (int i = 0; i < n_; ++i) {
        for (int k = csr_row_ptr_[i]; k < csr_row_ptr_[i + 1]; ++k) {
            csr_col_idx_[k] = pattern_[k].second;
        }
    }

    // 4) Slot lookup for resolveSlots/add/get.
    slot_lookup_.clear();
    slot_lookup_.reserve(static_cast<size_t>(nnz_) * 2);
    for (int k = 0; k < nnz_; ++k) {
        slot_lookup_.emplace(packKey(pattern_[k].first, pattern_[k].second), k);
    }

    // Release the temporary pattern (11x fewer bytes resident than COO map).
    pattern_.clear();
    pattern_.shrink_to_fit();

    finalized_ = true;
}

std::vector<int> DaeMatrix::resolveSlots(const std::vector<JacobianEntry>& entries) const {
    if (!finalized_) {
        throw std::runtime_error("resolveSlots before finalizePattern()");
    }
    std::vector<int> slots(entries.size());
    for (size_t i = 0; i < entries.size(); ++i) {
        auto it = slot_lookup_.find(packKey(entries[i].first, entries[i].second));
        if (it == slot_lookup_.end()) {
            char buf[192];
            std::snprintf(buf, sizeof(buf),
                          "resolveSlots: entry (%d,%d) not declared in pattern (nnz=%d)",
                          entries[i].first, entries[i].second, nnz_);
            throw std::runtime_error(buf);
        }
        slots[i] = it->second;
    }
    return slots;
}

void DaeMatrix::zeroValues() {
    // O(nnz) fill; structure is untouched.
    std::fill(csr_values_.begin(), csr_values_.end(), 0.0);
}

void DaeMatrix::accumulateSlots(const std::vector<int>& slots,
                                const double* values,
                                int count) {
    for (int i = 0; i < count; ++i) {
        csr_values_[slots[i]] += values[i];
    }
}

void DaeMatrix::add(int row, int col, double value) {
    if (!finalized_) {
        throw std::runtime_error("DaeMatrix::add before finalizePattern()");
    }
    auto it = slot_lookup_.find(packKey(row, col));
    if (it == slot_lookup_.end()) {
        throw std::runtime_error("DaeMatrix::add on undeclared entry");
    }
    csr_values_[it->second] += value;
}

double DaeMatrix::get(int row, int col) const {
    if (row < 0 || row >= n_ || col < 0 || col >= n_) {
        throw std::out_of_range("DaeMatrix::get out of range");
    }
    if (!finalized_) return 0.0;
    auto it = slot_lookup_.find(packKey(row, col));
    if (it == slot_lookup_.end()) return 0.0;
    return csr_values_[it->second];
}

void DaeMatrix::multiplyVector(const std::vector<double>& x,
                               std::vector<double>& y) const {
    if (static_cast<int>(x.size()) != n_) {
        throw std::invalid_argument("DaeMatrix::multiplyVector: size mismatch");
    }
    y.resize(n_, 0.0);

    // Tiny systems: serial (OpenMP region overhead dominates the work).
    constexpr int kParallelMin = 512;
    if (n_ < kParallelMin) {
        for (int row = 0; row < n_; ++row) {
            double sum = 0.0;
            const int start = csr_row_ptr_[row];
            const int end   = csr_row_ptr_[row + 1];
            for (int k = start; k < end; ++k) {
                sum += csr_values_[k] * x[csr_col_idx_[k]];
            }
            y[row] = sum;
        }
        return;
    }

    #pragma omp parallel for schedule(static)
    for (int row = 0; row < n_; ++row) {
        double sum = 0.0;
        const int start = csr_row_ptr_[row];
        const int end   = csr_row_ptr_[row + 1];
        for (int k = start; k < end; ++k) {
            sum += csr_values_[k] * x[csr_col_idx_[k]];
        }
        y[row] = sum;
    }
}

void DaeMatrix::printStats() const {
    const int n = n_ * n_;
    std::printf("DaeMatrix %dx%d (nnz=%d, %.3f%% dense, %s)\n",
                n_, n_, nnz_, 100.0 * nnz_ / n,
                finalized_ ? "finalized" : "pattern-phase");
}

} // namespace e_XSpice
