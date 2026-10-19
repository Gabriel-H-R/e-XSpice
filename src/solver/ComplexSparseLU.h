/**
 * @author Gabriel Henrique Silva
 * @file ComplexSparseLU.h
 * @brief Sparse COMPLEX LU for the small-signal AC system, with the
 *        factorization pattern computed once and reused at every frequency.
 *
 * Why this shape of solver (and how it improves on Xyce)
 * -----------------------------------------------------
 * Xyce (AnalysisPKG/N_ANP_AC.C) builds the AC system as a 2x2 REAL block
 * matrix
 *
 *     [  G   -w*C ] [Vr]   [Br]
 *     [  w*C   G  ] [Vi] = [Bi]
 *
 * which is the real/imaginary split of (G + j*w*C) X = B. That is 2N x 2N
 * REAL: it needs twice the memory of the complex system and, more
 * importantly, twice the factorization work (8/3 N^3 real flops instead of
 * 4/3 N^3 for a complex LU on N x N).
 *
 * This class factorizes the complex matrix DIRECTLY on the shared CSR
 * pattern of the DAE layer:
 *
 *   - G = dF/dx and C = dQ/dx evaluated ONCE at the DC operating point
 *     (exactly what the transient linearization already produces — the AC
 *     analysis therefore needs NO new device physics hooks for the matrix,
 *     only an excitation vector),
 *   - the SPARSITY PATTERN of A = G + j*w*C is frequency independent, so
 *     the symbolic factorization runs once and every frequency only refills
 *     the values: O(nnz) to build + O(nnz(L)) to factorize,
 *   - values are stored as two real arrays (structure-of-arrays), which
 *     vectorizes better than std::complex and halves the temporary traffic,
 *   - the triangular solves are independent per frequency, so the analysis
 *     can run FREQUENCIES in parallel (OpenMP / MPI) instead of a serial
 *     loop of factorizations.
 *
 * Method: left-looking column factorization A = L·U with U unit diagonal
 * (so the backward substitution needs no division), natural ordering, no
 * pivoting, over the combined L|U column pattern computed once.
 *
 * Stability: no pivoting keeps the pattern static (the property that makes
 * the whole reuse strategy work). Pivots are checked against a relative
 * threshold; when a pivot is numerically zero the diagonal is perturbed by
 * the smallest complex delta that fixes it, which is reported through
 * perturbedPivot(). Matrices from passive networks and device linearizations
 * (G diagonally dominant) never need it; a pathological case would need a
 * pivoted/distributed LU, which is the documented next step.
 */

#ifndef COMPLEX_SPARSE_LU_H
#define COMPLEX_SPARSE_LU_H

#include <algorithm>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Sparse complex LU over a fixed CSR pattern.
 *
 * Lifecycle:  setPattern() once -> factorize(omega, G, C) per frequency ->
 *             solve(...) per right-hand side. Everything between setPattern
 *             and the last solve is reused.
 */
class ComplexSparseLU {
public:
    ComplexSparseLU();

    // ---------------------------------------------------------------
    // SETUP (once per analysis)
    // ---------------------------------------------------------------

    /**
     * @brief Adopts a CSR pattern (row pointers + column indices).
     *
     * The same arrays serve the DC linearization: the LU keeps an internal
     * CSC (column) view plus the mapping back into the CSR value arrays, so
     * per frequency it only has to read G/C through that mapping.
     *
     * @param rowPtr CSR row pointers (size n+1)
     * @param colIdx CSR column indices (size nnz)
     */
    void setPattern(const std::vector<int>& rowPtr,
                    const std::vector<int>& colIdx);

    int  size() const { return n_; }
    int  nonZeros() const { return nnz_; }
    long long factorEntries() const { return ll_; }   // nnz(L) + nnz(U)

    /**
     * @brief True when the factor pattern (combined L|U of a column)
     *        contains (row, col). Diagnostics/tests: lets a caller verify
     *        that the symbolic fill-in covers what a dense reference LU
     *        touches.
     */
    bool hasFactorEntry(int row, int col) const {
        if (row < 0 || col < 0 || col >= n_) return false;
        return std::binary_search(lcol_idx_.begin() + lcol_ptr_[col],
                                  lcol_idx_.begin() + lcol_ptr_[col + 1],
                                  row);
    }

    /**
     * @brief Complex values of the stored factors at (row, col).
     * @return false when the entry is outside the factor pattern.
     * Diagnostics/tests only (lets a caller check L·U == A).
     */
    bool factorValues(int row, int col,
                      double& l_re, double& l_im,
                      double& u_re, double& u_im) const {
        const auto b = lcol_idx_.begin() + lcol_ptr_[col];
        const auto e = lcol_idx_.begin() + lcol_ptr_[col + 1];
        const auto it = std::lower_bound(b, e, row);
        if (it == e || *it != row) return false;
        const int p = static_cast<int>(it - lcol_idx_.begin());
        l_re = lx_re_[static_cast<size_t>(p)];
        l_im = lx_im_[static_cast<size_t>(p)];
        u_re = ux_re_[static_cast<size_t>(p)];
        u_im = ux_im_[static_cast<size_t>(p)];
        return true;
    }

    /**
     * @brief Relative pivot threshold (a pivot below |tol|·|a|max is treated
     *        as singular and perturbed). Default 1e-13.
     */
    void setPivotTolerance(double tol) { pivot_tol_ = tol; }

    // ---------------------------------------------------------------
    // SOLVE PATH
    // ---------------------------------------------------------------

    /**
     * @brief Factorizes A = G + j·omega·C.
     *
     * @param omega  2*pi*f
     * @param G      real conductance values, on the CSR pattern
     * @param C      real charge (capacitance) values, same pattern
     * @return false only when the matrix cannot be factorized even with
     *         diagonal perturbation
     */
    bool factorize(double omega, const double* G, const double* C);

    /**
     * @brief Solves A x = b with the stored factors.
     * @param xr,xi right-hand side on entry, solution on exit (length n)
     */
    void solve(std::vector<double>& xr, std::vector<double>& xi) const;

    /**
     * @brief Solves the ADJOINT system A^T x = b with the stored factors.
     *
     * This is what makes a noise sweep affordable: the transfer from EVERY
     * noise generator to the output is a row of A^-1, and A^-T·b gives all of
     * them from ONE solve (symmetry of A), instead of one solve per generator.
     *
     * With A = L·U, A^T = U^T·L^T:
     *   - U^T is upper triangular with diagonal U(j,j), so a forward sweep in
     *     column order works; the entries of U^T(i,j) = U(j,i) live in column
     *     i BELOW its diagonal — the compact-LU subdiagonal slot shared with L;
     *   - L^T is upper triangular with UNIT diagonal, so a backward sweep in
     *     reverse column order finishes with no division.
     * Cost: two passes over nnz(L) plus n complex divisions — no extra
     * factorization, no transposed copy of the factors.
     *
     * @param xr,xi right-hand side on entry, solution on exit (length n)
     */
    void solveTranspose(std::vector<double>& xr, std::vector<double>& xi) const;

    // ---------------------------------------------------------------
    // DIAGNOSTICS
    // ---------------------------------------------------------------

    bool perturbedPivot() const { return perturbed_; }
    double minPivotMagnitude() const { return min_pivot_; }
    long long lastFactorFlops() const { return last_flops_; }
    const std::string& message() const { return message_; }

private:
    void buildColumnView();       // CSR -> CSC + slot mapping
    void symbolicFactorization(); // fill pattern of L and U (once)
    void allocateFactors();

    // Pattern (CSR of the input, CSC view of the same entries)
    int n_ = 0;
    int nnz_ = 0;
    std::vector<int> row_ptr_;       // input CSR (kept for interface clarity)
    std::vector<int> col_idx_;       // input CSR
    std::vector<int> col_ptr_;       // CSC column pointers (n+1)
    std::vector<int> row_of_;        // CSC row indices (nnz)
    std::vector<int> value_slot_;    // CSC entry -> CSR value index

    // Combined L|U column pattern (sorted row indices per column)
    std::vector<int> lcol_ptr_;      // n+1
    std::vector<int> lcol_idx_;      // sum of column lengths

    // Numeric factors
    std::vector<double> ux_re_, ux_im_;   // U values (rows >= column)
    std::vector<double> lx_re_, lx_im_;   // L values (rows < column, unit diagonal)
    std::vector<int> udiag_;              // index of U(j,j) inside column j
    std::vector<int> lcol_split_;         // column j: first index with row >= j

    // Per-factorization work
    std::vector<double> acc_re_, acc_im_; // dense accumulator (n)
    std::vector<int> acc_pos_;            // row -> position in the column list
    std::vector<int> acc_gen_;            // row -> column that last claimed it
    std::vector<double> g_work_;          // mutable G copy (perturbation)
    std::vector<int> xi_, mark_, stack_;  // symbolic work arrays

    double pivot_tol_ = 1e-13;
    bool   perturbed_ = false;
    double min_pivot_ = 0.0;
    long long last_flops_ = 0;
    long long ll_ = 0;
    std::string message_;
};

/**
 * @brief Static row permutation that gives every column a structural diagonal.
 *
 * MNA matrices (the ones a small-signal system always is) contain rows with
 * NO self-term: a node wired only to voltage sources and/or inductors has a
 * KCL row whose unknowns are branch currents, e.g.
 *
 *     row(node in):  i(V1) + i(L1) = 0        ->  diagonal entry is zero
 *
 * A factorization without pivoting has nothing to pivot on there. Adding a
 * tiny conductance would work numerically and change the circuit; perturbing
 * the diagonal leaks into every other unknown (measured: 5e-7 error on a node
 * that an ideal source pins exactly). The correct fix is an ORDERING.
 *
 * build() finds a permutation that assigns to every column a row it has a
 * USABLE entry in — "usable" being judged on the values of the matrix (see
 * below), not merely on the pattern. Row permutation creates NO fill-in (the
 * pattern is just relabelled), so the symbolic factorization still runs once
 * and every frequency reuses it.
 *
 * Why the values and not only the pattern: the DAE layer declares a diagonal
 * slot for EVERY row (the node-gmin / LU-robustness convention), so the
 * pattern claims a diagonal exists even where no device stamps one. The
 * decision therefore looks at the numbers: a column whose own diagonal is
 * negligible next to its largest entry is one that cannot carry its own pivot,
 * and it is given a partner row instead.
 *
 * The ordering is computed ONCE from the first frequency's matrix and then
 * frozen, which is what keeps the static pattern valid for the whole sweep:
 * in a circuit matrix, which entries exist and which are negligible does not
 * change with w (a device either stamps a self-conductance or it does not).
 * A later frequency whose pivot turns out to be numerically unusable is
 * repaired by the LU's diagonal perturbation, which is reported.
 *
 * The greedy pass prefers the natural diagonal whenever it is usable (a
 * conductance/gm diagonal is the well-conditioned pivot in a circuit matrix)
 * and, for the columns that have none, the highest-priority free row of the
 * pattern. The caller supplies the priority: for MNA, BRANCH rows are
 * preferred, because their couplings to node rows are the ±1 of a branch
 * equation — a pivot of magnitude one, always well scaled.
 */
struct StaticRowPermutation {
    std::vector<int> row_of_pos;   // row_of_pos[k] = source row placed at position k
    bool ok = false;               // false: the pattern admits no such ordering
    int  swaps = 0;                // how many columns are NOT on their own row
    int  weak_diagonals = 0;       // columns whose own diagonal is negligible
    std::string reason;            // filled when ok == false

    /**
     * @param re,im values of the reference matrix on the CSR pattern
     *        (the first frequency of the sweep).
     */
    static StaticRowPermutation build(const std::vector<int>& rowPtr,
                                      const std::vector<int>& colIdx,
                                      const std::vector<double>& re,
                                      const std::vector<double>& im,
                                      const std::vector<int>& priority);
};

} // namespace e_XSpice

#endif // COMPLEX_SPARSE_LU_H