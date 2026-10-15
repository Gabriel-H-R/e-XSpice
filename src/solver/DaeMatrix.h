/**
 * @author Gabriel Henrique Silva
 * @file DaeMatrix.h
 * @brief Fixed-sparsity-pattern sparse matrix for DAE-formulation assembly.
 *
 * Motivation (performance/scalability):
 * The legacy SparseMatrix accumulates into a std::map<std::pair<int,int>,double>
 * (O(log n) per insert with pointer chasing) and re-sorts every row when
 * converting COO -> CSR at every Newton iteration / time step. For a
 * transient simulation with NR iterations this dominates the assembly cost.
 *
 * DaeMatrix separates SETUP from SOLVE:
 *
 *  SETUP (once, after indices are assigned):
 *      DaeMatrix m(n);
 *      m.declareEntry(r, c);      // called by every device (Jacobian stamp)
 *      m.finalizePattern();       // -> CSR structure + O(1) slot map
 *      slots = m.resolveSlots(device_entries);   // device-owned nnz positions
 *
 *  SOLVE (hot path, per iteration):
 *      m.zeroValues();            // O(nnz), no tree, no sort
 *      for each device:
 *          for each entry i:  m.values()[slot[i]] += value[i];   // O(1), cache-friendly
 *
 * Because every nnz slot belongs to exactly one device registration, the
 * assembly loop is embarrassingly parallel (OpenMP) with zero synchronization:
 * this is where the legacy code was serial/map-bound.
 */

#ifndef DAE_MATRIX_H
#define DAE_MATRIX_H

#include <vector>
#include <cstdint>
#include <utility>
#include <unordered_map>
#include "MatrixOps.h"

namespace e_XSpice {

/**
 * @brief (row, col) coordinate pair used in the pattern phase.
 */
using JacobianEntry = std::pair<int, int>;

class DaeMatrix : public MatrixOps {
public:
    /**
     * @brief Constructor
     * @param n Number of rows/columns (square system, MNA). n ≤ 0 creates an
     *          empty placeholder (useful for default-constructed members).
     */
    DaeMatrix() = default;
    explicit DaeMatrix(int n);

    ~DaeMatrix() = default;

    /**
     * @brief Number of rows (= columns). Systems are always square.
     */
    int rows() const override { return n_; }
    int cols() const { return n_; }

    // ---------------------------------------------------------------
    // SETUP PHASE
    // ---------------------------------------------------------------

    /**
     * @brief Declares a structural non-zero (must be called before finalize).
     * Duplicated declarations are collapsed during finalizePattern().
     */
    void declareEntry(int row, int col);

    /**
     * @brief Declares all entries of a device stamp (row list is symmetric
     *        with col list).
     */
    void declareEntries(const std::vector<int>& rows,
                        const std::vector<int>& cols);

    /**
     * @brief Freezes the structure and builds CSR arrays + slot lookup.
     * Required before resolveSlots/accumulate.
     */
    void finalizePattern();

    /**
     * @brief Returns true if pattern is frozen.
     */
    bool isFinalized() const { return finalized_; }

    /**
     * @brief Maps a device's (row,col) entries to nnz slot positions inside
     *        the CSR value array. The device then owns those slots.
     *
     * @return slot positions, same order/length as entries
     */
    std::vector<int> resolveSlots(const std::vector<JacobianEntry>& entries) const;

    /**
     * @brief Total structural non-zeros.
     */
    int    nonZeros() const { return nnz_; }
    int    capacity() const { return static_cast<int>(csr_values_.capacity()); }

    // ---------------------------------------------------------------
    // SOLVE PHASE (hot path)
    // ---------------------------------------------------------------

    /**
     * @brief Resets all values to 0.0 (keeps structure). O(nnz), no sorting.
     */
    void zeroValues();

    /**
     * @brief Accumulates value into a slot obtained from resolveSlots.
     *        O(1). Safe: each slot belongs to a single device.
     */
    void accumulateSlot(int slot, double value) { csr_values_[slot] += value; }

    /**
     * @brief Accumulates a scalar-set of slots with per-slot values.
     *        Used by devices with a small fixed stamp (2x2 blocks).
     */
    void accumulateSlots(const std::vector<int>& slots,
                         const double* values,
                         int count);

    /**
     * @brief Convenience RSS: accumulate value at (row,col) via hash lookup.
     *        Slower than accumulateSlot; reserved for setup/debug.
     */
    void add(int row, int col, double value);

    /**
     * @brief Reads a value (debug/tests).
     */
    double get(int row, int col) const override;

    // ---------------------------------------------------------------
    // LINEAR ALGEBRA
    // ---------------------------------------------------------------

    /**
     * @brief y = A * x (OpenMP parallel over rows).
     */
    void multiplyVector(const std::vector<double>& x,
                        std::vector<double>& y) const override;

    /**
     * @brief Computes only rows [rowBegin, rowEnd) of y = A*x.
     *        Rows outside the block are left untouched. Used by the
     *        row-distributed Krylov solver: each rank evaluates its own
     *        block from the replicated matrix and the replicated x, then
     *        the blocks are fused with one Allgatherv.
     */
    void multiplyVectorBlock(const std::vector<double>& x,
                             std::vector<double>& y,
                             int rowBegin, int rowEnd) const;

    // MatrixOps interface
    bool isSymmetric() const override { return false; }

    /**
     * @brief Direct access to CSR arrays (for solvers / debugging).
     */
    const std::vector<double>& values() const { return csr_values_; }
    std::vector<double>&       values()       { return csr_values_; }
    const std::vector<int>&    colIndices() const { return csr_col_idx_; }
    const std::vector<int>&    rowPointers() const { return csr_row_ptr_; }

    /**
     * @brief Debug: prints non-zeros.
     */
    void printStats() const;

private:
    int n_;                    // dimension (square)
    bool finalized_ = false;
    int  nnz_ = 0;

    // Structure
    std::vector<JacobianEntry> pattern_;              // raw declarations (dedup at finalize)
    std::vector<int>            csr_row_ptr_;        // size n_+1
    std::vector<int>            csr_col_idx_;        // size nnz_
    std::vector<double>         csr_values_;         // size nnz_

    // (row,col) -> nnz slot (packed key) for resolveSlots/add/get
    std::unordered_map<std::uint64_t, int> slot_lookup_;

    static std::uint64_t packKey(int row, int col) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(row)) << 32)
             ^ static_cast<std::uint32_t>(col);
    }
};

} // namespace e_XSpice

#endif // DAE_MATRIX_H
